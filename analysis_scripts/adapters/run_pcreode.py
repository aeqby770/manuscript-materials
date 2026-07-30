"""pCreode adapter: run pCreode and map its best graph to all cells.

Usage: python run_pcreode.py <canonical_h5ad> <output_csv>
"""
import importlib
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from itertools import combinations
from pathlib import Path

CONFIG = {
    "paths": {
        "vendor_root": Path(
            os.environ.get(
                "LLR_VENDOR_ROOT",
                Path(__file__).resolve().parents[3]
                / "LineaLinkR_benchmark_assets"
                / "02_vendor_algorithms",
            )
        ),
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
    "memory_budget_gb": 55.0,
    "native_num_runs": 100,
    "graph_worker_limit": 4,
    "worker_memory_fraction": 0.75,
    "worker_allocator_headroom": 1.20,
    "vendor_live_float64_square_arrays": 6,
    "vendor_live_intp_square_arrays": 1,
    "vendor_symmetrized_float64_square_arrays": 1,
    "vendor_csr_upper_triangle_fraction": 0.5,
    "worker_full_data_copies": 3,
    "worker_downsample_data_copies": 2,
    "worker_density_vector_copies": 4,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])

import matplotlib.pyplot as plt
import numpy as np
if not hasattr(np, 'int'):
    np.int = int
    np.float = float
import pandas as pd
import psutil
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import MinMaxScaler
from threadpoolctl import threadpool_limits

np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])

import igraph

_orig_weighted_adj = igraph.Graph.Weighted_Adjacency


def _sym_weighted_adj(matrix, *args, **kwargs):
    m = np.asarray(matrix)
    del matrix
    m = m + m.T
    m *= 0.5
    return _orig_weighted_adj(m, *args, **kwargs)
igraph.Graph.Weighted_Adjacency = staticmethod(_sym_weighted_adj)

sys.path.insert(0, str(CONFIG["paths"]["vendor_root"] / "pCreode"))
import pcreode

# Multiple pCreode centroids can map to one observed cell. Normalize endpoint
# indices and restore repeated distance columns for current igraph releases.
_PCREODE_FUNCTIONS = importlib.import_module("pcreode.functions")
_PCREODE_CORE = importlib.import_module("pcreode.pcreode")
_ORIGINAL_GRAPH_DISTANCE = _PCREODE_FUNCTIONS.get_graph_distance
_ORIGINAL_HIERARCHICAL_PLACEMENT = _PCREODE_FUNCTIONS.hierarchical_placement


def _graph_distance_with_repeated_targets(from_ind, to_ind, graph):
    target_array = np.asarray(to_ind, dtype=int)
    unique_targets, inverse = np.unique(target_array, return_inverse=True)
    unique_distances = _ORIGINAL_GRAPH_DISTANCE(
        from_ind,
        unique_targets,
        graph,
    )
    return unique_distances[:, inverse]


def _hierarchical_placement_with_unique_endpoints(graph, end_ind):
    endpoint_array = np.asarray(end_ind, dtype=int)
    endpoint_array = endpoint_array[
        (endpoint_array >= 0) & (endpoint_array < graph.vcount())
    ]
    unique_endpoints = np.unique(endpoint_array)
    return _ORIGINAL_HIERARCHICAL_PLACEMENT(graph, unique_endpoints)


_PCREODE_FUNCTIONS.get_graph_distance = _graph_distance_with_repeated_targets
_PCREODE_CORE.get_graph_distance = _graph_distance_with_repeated_targets
pcreode.get_graph_distance = _graph_distance_with_repeated_targets
_PCREODE_FUNCTIONS.hierarchical_placement = _hierarchical_placement_with_unique_endpoints
_PCREODE_CORE.hierarchical_placement = _hierarchical_placement_with_unique_endpoints
pcreode.hierarchical_placement = _hierarchical_placement_with_unique_endpoints

import _mem_guard
_mem_guard.cap_self(CONFIG["memory_budget_gb"])
_PROCESS_BASE_MEMORY = psutil.Process().memory_info()
_PROCESS_BASE_PRIVATE_BYTES = int(
    getattr(
        _PROCESS_BASE_MEMORY,
        "private",
        _PROCESS_BASE_MEMORY.rss,
    )
)

_GRAPH_DATA = None
_GRAPH_DENSITY = None
_GRAPH_NOISE = None
_GRAPH_TARGET = None
_GRAPH_TMPDIR = None
_GRAPH_POTENTIAL_CLUSTERS = None
_THREADPOOL_LIMITER = None
_SCORE_GRAPHS = None


def _prepare_pca(h5ad_path):
    import scanpy as sc

    adata = sc.read(h5ad_path)
    adata.X = (
        adata.layers["counts"].copy()
        if "counts" in adata.layers
        else adata.X.copy()
    )
    adata.uns.pop("log1p", None)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=min(2000, adata.n_vars),
    )
    matrix = adata[:, adata.var.highly_variable].X
    matrix = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
    matrix = np.nan_to_num(matrix, nan=0, posinf=0, neginf=0).astype(
        np.float64,
        copy=False,
    )
    n_components = min(20, matrix.shape[1], matrix.shape[0] - 1)
    data_pca = PCA(
        n_components=n_components,
        random_state=CONFIG["seed"],
    ).fit_transform(matrix).astype(np.float64, copy=False)
    return adata, data_pca


def _exact_density(data_pca):
    neighbors = NearestNeighbors(
        n_neighbors=min(10, data_pca.shape[0]),
        n_jobs=CONFIG["threads"],
    ).fit(data_pca)
    radius_guess = np.percentile(
        neighbors.kneighbors(data_pca, return_distance=True)[0][:, -1],
        50,
    )
    radius_graph = neighbors.radius_neighbors_graph(
        data_pca,
        radius=radius_guess,
        mode="connectivity",
    )
    density = np.asarray(radius_graph.getnnz(axis=1), dtype=int) - 1
    return density, radius_guess


def _downsample_probabilities(density, noise, target):
    outlier_density = np.percentile(density, noise)
    target_density = np.percentile(density, target)
    return np.select(
        (
            density < outlier_density,
            density <= target_density,
            density > target_density,
        ),
        (
            0.0,
            1.0,
            np.divide(
                target_density,
                density,
                out=np.zeros(density.shape, dtype=float),
                where=density != 0,
            ),
        ),
    )


def _exact_downsample_sizes(density, noise, target, run_seeds):
    probabilities = _downsample_probabilities(density, noise, target)
    return np.fromiter(
        map(
            lambda seed: np.count_nonzero(
                probabilities
                >= np.random.RandomState(int(seed)).random_sample(
                    probabilities.size
                )
            ),
            run_seeds,
        ),
        dtype=np.int64,
        count=len(run_seeds),
    )


def _graph_worker_count(
    data_pca,
    density,
    downsample_sizes,
    num_runs,
    system_available_bytes=None,
    parent_private_bytes=None,
):
    actual_downsample_sizes = np.asarray(downsample_sizes, dtype=np.int64)
    downsample_cells = int(np.max(actual_downsample_sizes))
    square_elements = downsample_cells * downsample_cells
    float64_bytes = np.dtype(np.float64).itemsize
    intp_bytes = np.dtype(np.intp).itemsize
    int32_bytes = np.dtype(np.int32).itemsize
    pointer_bytes = np.dtype(np.intp).itemsize
    python_float_bytes = sys.getsizeof(0.0)
    python_list_header_bytes = sys.getsizeof([])
    # Peak vendor state while igraph converts a symmetrized dense adjacency.
    live_square_bytes = square_elements * (
        CONFIG["vendor_live_float64_square_arrays"] * float64_bytes
        + CONFIG["vendor_live_intp_square_arrays"] * intp_bytes
        + CONFIG["vendor_symmetrized_float64_square_arrays"] * float64_bytes
        + pointer_bytes
        + python_float_bytes
    )
    csr_upper_bytes = int(
        square_elements
        * CONFIG["vendor_csr_upper_triangle_fraction"]
        * (float64_bytes + int32_bytes)
        + (downsample_cells + 1) * int32_bytes
    )
    nested_list_headers_bytes = (
        (downsample_cells + 1) * python_list_header_bytes
    )
    linear_data_bytes = int(
        CONFIG["worker_full_data_copies"] * data_pca.nbytes
        + CONFIG["worker_downsample_data_copies"]
        * downsample_cells
        * data_pca.shape[1]
        * float64_bytes
        + CONFIG["worker_density_vector_copies"] * density.nbytes
    )
    per_worker_bytes = int(
        np.ceil(
            (
                live_square_bytes
                + csr_upper_bytes
                + nested_list_headers_bytes
                + linear_data_bytes
                + _PROCESS_BASE_PRIVATE_BYTES
            )
            * CONFIG["worker_allocator_headroom"]
        )
    )
    process_memory = psutil.Process().memory_info()
    current_parent_private_bytes = int(
        getattr(process_memory, "private", process_memory.rss)
        if parent_private_bytes is None
        else parent_private_bytes
    )
    current_system_available_bytes = int(
        psutil.virtual_memory().available
        if system_available_bytes is None
        else system_available_bytes
    )
    job_remaining_bytes = max(
        0,
        int(CONFIG["memory_budget_gb"] * (1024 ** 3))
        - current_parent_private_bytes,
    )
    pool_capacity_bytes = min(
        current_system_available_bytes,
        job_remaining_bytes,
    )
    if per_worker_bytes > pool_capacity_bytes:
        gib = float(1024 ** 3)
        raise MemoryError(
            "pCreode graph preflight: "
            f"downsample_cells={downsample_cells}, "
            f"required_per_worker={per_worker_bytes / gib:.2f} GiB, "
            f"capacity={pool_capacity_bytes / gib:.2f} GiB"
        )
    pool_budget_bytes = int(
        pool_capacity_bytes * CONFIG["worker_memory_fraction"]
    )
    memory_workers = max(1, pool_budget_bytes // per_worker_bytes)
    workers = min(
        CONFIG["threads"],
        CONFIG["graph_worker_limit"],
        num_runs,
        memory_workers,
    )
    gib = float(1024 ** 3)
    return {
        "workers": workers,
        "downsample_min": int(np.min(actual_downsample_sizes)),
        "downsample_median": float(np.median(actual_downsample_sizes)),
        "downsample_max": downsample_cells,
        "per_worker_gb": per_worker_bytes / gib,
        "system_available_gb": current_system_available_bytes / gib,
        "parent_private_gb": current_parent_private_bytes / gib,
        "job_remaining_gb": job_remaining_bytes / gib,
        "pool_budget_gb": pool_budget_bytes / gib,
        "estimated_pool_peak_gb": workers * per_worker_bytes / gib,
        "square_bytes_per_cell_pair": (
            live_square_bytes + csr_upper_bytes
        ) / square_elements,
    }


def _init_graph_worker(
    data_pca,
    density,
    noise,
    target,
    tmpdir,
    potential_clusters,
):
    global _GRAPH_DATA
    global _GRAPH_DENSITY
    global _GRAPH_NOISE
    global _GRAPH_TARGET
    global _GRAPH_TMPDIR
    global _GRAPH_POTENTIAL_CLUSTERS
    global _THREADPOOL_LIMITER
    _GRAPH_DATA = data_pca
    _GRAPH_DENSITY = density
    _GRAPH_NOISE = noise
    _GRAPH_TARGET = target
    _GRAPH_TMPDIR = tmpdir
    _GRAPH_POTENTIAL_CLUSTERS = potential_clusters
    _THREADPOOL_LIMITER = threadpool_limits(limits=1)


def _run_graph(run_spec):
    run_id, run_seed = run_spec
    np.random.seed(run_seed)
    pcreode.pCreode(
        _GRAPH_DATA,
        _GRAPH_DENSITY,
        _GRAPH_NOISE,
        _GRAPH_TARGET,
        file_path=_GRAPH_TMPDIR + os.sep,
        num_runs=1,
        potential_clusters=_GRAPH_POTENTIAL_CLUSTERS,
        cls_thresh=0.0,
        start_id=run_id,
        mute=True,
    )
    return run_id


def _branch_path_matrix(graph):
    node_count = graph.vcount()
    degrees = np.asarray(graph.degree(), dtype=float)
    matrix = np.zeros((node_count, node_count), dtype=float)
    upper = np.triu_indices(node_count, 1)
    paths = map(
        lambda pair: graph.get_all_shortest_paths(pair[0], pair[1])[0],
        zip(upper[0], upper[1]),
    )
    values = np.fromiter(
        map(
            lambda path: np.sum(degrees[path][degrees[path] > 2.0])
            - 2.0 * np.sum(degrees[path] > 2.0),
            paths,
        ),
        dtype=float,
        count=upper[0].size,
    )
    matrix[upper] = values
    matrix[(upper[1], upper[0])] = values
    return matrix, upper


def _read_scoring_graph(data_pca, tmpdir, graph_id):
    indices = np.atleast_1d(
        np.genfromtxt(
            os.path.join(tmpdir, f"ind_{graph_id}.csv"),
            delimiter=",",
        )
    ).astype(int)
    adjacency = pd.read_csv(
        os.path.join(tmpdir, f"adj_{graph_id}.txt"),
        sep=" ",
        header=None,
    ).values
    coordinates = data_pca[indices]
    euclidean = pairwise_distances(coordinates, coordinates, n_jobs=1, metric="l2")
    weighted = np.multiply(euclidean, adjacency)
    graph = igraph.Graph.Weighted_Adjacency(
        weighted.tolist(),
        mode=igraph.ADJ_UNDIRECTED,
    )
    graph_distances = np.asarray(graph.distances(weights="weight"), dtype=float)
    branches, upper = _branch_path_matrix(graph)
    return {
        "indices": indices,
        "adjacency": adjacency,
        "coordinates": coordinates,
        "distances": graph_distances,
        "branches": branches,
        "upper": upper,
    }


def _directional_graph_difference(source, target, cross_distances):
    mapping = np.argmin(cross_distances, axis=1)
    nearest = cross_distances[np.arange(mapping.size), mapping]
    row, col = source["upper"]
    mapped_row = mapping[row]
    mapped_col = mapping[col]
    branch_difference = np.mean(
        np.abs(
            source["branches"][row, col]
            - target["branches"][mapped_row, mapped_col]
        )
    )
    distance_difference = np.mean(
        np.abs(
            target["distances"][mapped_row, mapped_col]
            - source["distances"][row, col]
            + nearest[row]
            + nearest[col]
        )
    )
    return branch_difference, distance_difference


def _score_graph_pair(pair):
    left_id, right_id = pair
    left = _SCORE_GRAPHS[left_id]
    right = _SCORE_GRAPHS[right_id]
    cross = pairwise_distances(
        left["coordinates"],
        right["coordinates"],
        n_jobs=1,
        metric="l2",
    )
    left_score = _directional_graph_difference(left, right, cross)
    right_score = _directional_graph_difference(right, left, cross.T)
    return (
        left_id,
        right_id,
        max(left_score[0], right_score[0]),
        max(left_score[1], right_score[1]),
    )


def _score_graphs(data_pca, tmpdir, num_runs):
    global _SCORE_GRAPHS
    _SCORE_GRAPHS = list(
        map(
            lambda graph_id: _read_scoring_graph(data_pca, tmpdir, graph_id),
            range(num_runs),
        )
    )
    graph_pairs = list(combinations(range(num_runs), 2))
    score_workers = min(CONFIG["threads"], len(graph_pairs))
    with ThreadPoolExecutor(max_workers=score_workers) as executor:
        scores = list(executor.map(_score_graph_pair, graph_pairs))
    scores = np.asarray(scores, dtype=float)
    rows = scores[:, 0].astype(int)
    columns = scores[:, 1].astype(int)
    branch_difference = np.zeros((num_runs, num_runs), dtype=float)
    distance_difference = np.zeros((num_runs, num_runs), dtype=float)
    branch_difference[rows, columns] = scores[:, 2]
    distance_difference[rows, columns] = scores[:, 3]
    branch_normalized = MinMaxScaler(feature_range=(0, 1)).fit_transform(
        branch_difference
    )
    distance_normalized = MinMaxScaler(feature_range=(0, 1)).fit_transform(
        distance_difference
    )
    combined = (
        branch_normalized
        + branch_normalized.T
        + distance_normalized
        + distance_normalized.T
    )
    np.savetxt(
        os.path.join(tmpdir, "branch_diff.csv"),
        branch_difference + branch_difference.T,
        delimiter=",",
    )
    np.savetxt(
        os.path.join(tmpdir, "graph_dist_diff.csv"),
        distance_difference + distance_difference.T,
        delimiter=",",
    )
    np.savetxt(
        os.path.join(tmpdir, "combined_norm_diff.csv"),
        combined,
        delimiter=",",
    )
    ranks = np.argsort(np.mean(combined, axis=0))
    return ranks, _SCORE_GRAPHS

def run_pcreode(h5ad_path, output_csv, num_runs=None, potential_clusters=10):
    num_runs = CONFIG["native_num_runs"] if num_runs is None else int(num_runs)
    adata, data_pca = _prepare_pca(h5ad_path)
    density, radius_guess = _exact_density(data_pca)
    noise, target = pcreode.get_thresholds(data_pca)
    tmpdir = str(Path(output_csv).parent / "pcreode_tmp")
    shutil.rmtree(tmpdir, ignore_errors=True)
    os.makedirs(tmpdir, exist_ok=True)
    run_seeds = np.random.SeedSequence(CONFIG["seed"]).generate_state(
        num_runs,
        dtype=np.uint32,
    )
    downsample_sizes = _exact_downsample_sizes(
        density,
        noise,
        target,
        run_seeds,
    )
    graph_schedule = _graph_worker_count(
        data_pca,
        density,
        downsample_sizes,
        num_runs,
    )
    run_specs = zip(range(num_runs), map(int, run_seeds))
    print(
        f"pCreode graph generation: runs={num_runs}, workers={graph_schedule['workers']}, "
        f"downsample_min={graph_schedule['downsample_min']}, "
        f"downsample_median={graph_schedule['downsample_median']:.1f}, "
        f"downsample_max={graph_schedule['downsample_max']}, "
        f"per_worker_estimate={graph_schedule['per_worker_gb']:.2f} GB, "
        f"pool_peak_estimate={graph_schedule['estimated_pool_peak_gb']:.2f} GB, "
        f"pool_budget={graph_schedule['pool_budget_gb']:.2f} GB, "
        f"system_available={graph_schedule['system_available_gb']:.2f} GB, "
        f"parent_private={graph_schedule['parent_private_gb']:.2f} GB, "
        f"density_radius={radius_guess:.6g}",
        flush=True,
    )
    with ProcessPoolExecutor(
        max_workers=graph_schedule["workers"],
        initializer=_init_graph_worker,
        initargs=(
            data_pca,
            density,
            noise,
            target,
            tmpdir,
            potential_clusters,
        ),
    ) as executor:
        completed = list(executor.map(_run_graph, run_specs))
    if completed != list(range(num_runs)):
        raise RuntimeError("pCreode graph run IDs are incomplete or out of order")
    ranks, scoring_graphs = _score_graphs(data_pca, tmpdir, num_runs)
    best_graph_id = int(ranks[0])
    best_graph = scoring_graphs[best_graph_id]
    indices = best_graph["indices"]
    adj = best_graph["adjacency"]
    distances = dijkstra(
        csr_matrix(adj > 0), directed=False, indices=0, unweighted=True
    )
    max_d = distances[distances < np.inf].max()
    if max_d > 0:
        distances = distances / max_d
    nearest_node = (
        NearestNeighbors(n_neighbors=1)
        .fit(data_pca[indices])
        .kneighbors(data_pca, return_distance=False)
        .ravel()
    )
    pseudotime = distances[nearest_node]
    out = pd.DataFrame({
        "cell_id": adata.obs_names.astype(str),
        "pseudotime": pseudotime,
    })
    out.to_csv(output_csv, index=False)
    shutil.rmtree(tmpdir, ignore_errors=True)
    n_valid = np.sum(np.isfinite(pseudotime))
    print(
        f"pCreode OK: {adata.n_obs} cells, {n_valid} assigned by nearest graph "
        f"node, best_graph={best_graph_id}"
    )
    return pseudotime

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_pcreode.py <canonical_h5ad> <output_csv>")
        sys.exit(1)
    run_pcreode(sys.argv[1], sys.argv[2])
