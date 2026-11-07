#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

struct StatePatch {
    std::vector<int> cells;
    std::vector<int> owned;
    int state = 0;
    int offset = 0;
    int rank = 0;
    double variance_fraction = 0;
    double gap_fraction = 0;
};

// Components are determined by observed cell edges before splitting chart labels.
struct Components {
    std::vector<int> parent;
    explicit Components(int size) : parent(size) {
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int cell) {
        int root = cell;
        while (parent[root] != root) root = parent[root];
        while (parent[cell] != cell) {
            int next = parent[cell];
            parent[cell] = root;
            cell = next;
        }
        return root;
    }
    void join(int left, int right) {
        left = find(left);
        right = find(right);
        if (left != right) parent[std::max(left, right)] = std::min(left, right);
    }
};

constexpr int chart_info_columns = 8;

} // namespace

// Calibration charts are local expression models, not biological lineage nodes.
// The reference forest is supplied by the existing rooted geometry. Its edges
// need not order state means. Calibration cannot worsen their anchor reversals;
// every cell keeps one shared time, including trunk cells.
extern "C" int llr_state_clock(
    int n_cells, int n_dims, const double* expression,
    const int* outer, const int* inner, const double* edge_weight,
    const int* input_state, const double* anchor_order,
    int n_states, const int* reference_parent, int threads,
    double* output_time, int* output_chart, int* output_component,
    double* output_info, double* output_diagnostics, int* output_ncharts,
    char* error, int error_size) {
    try {
        Eigen::setNbThreads(threads);
        using Sparse = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
        using Triplet = Eigen::Triplet<double>;
        Eigen::Map<const Eigen::MatrixXd> observed(expression, n_cells, n_dims);
        Eigen::Map<const Eigen::VectorXd> anchor(anchor_order, n_cells);
        if (!observed.allFinite() || !anchor.allFinite()) {
            throw std::runtime_error("State-clock expression or global coordinate is not finite");
        }
        Components components(n_cells);
        for (int cell = 0; cell < n_cells; ++cell) {
            for (int edge = outer[cell]; edge < outer[cell + 1]; ++edge) {
                if (edge_weight[edge] > 0) components.join(cell, inner[edge]);
            }
        }
        std::map<int, int> component_ids;
        std::map<std::pair<int, int>, int> chart_ids;
        std::vector<StatePatch> patches;
        std::vector<int> chart_component;
        std::vector<int> state_size(n_states, 0);
        for (int cell = 0; cell < n_cells; ++cell) {
            const int state = input_state[cell] - 1;
            if (state < 0 || state >= n_states) throw std::runtime_error("Invalid reference state");
            ++state_size[state];
            auto component = component_ids.emplace(components.find(cell), component_ids.size()).first->second;
            auto key = std::make_pair(component, state);
            auto found = chart_ids.find(key);
            if (found == chart_ids.end()) {
                found = chart_ids.emplace(key, patches.size()).first;
                patches.emplace_back();
                patches.back().state = state;
                chart_component.push_back(component);
            }
            const int chart = found->second;
            output_chart[cell] = chart + 1;
            patches[chart].owned.push_back(cell);
            patches[chart].cells.push_back(cell);
        }
        const int n_charts = patches.size();
        *output_ncharts = n_charts;
        for (int cell = 0; cell < n_cells; ++cell) {
            for (int edge = outer[cell]; edge < outer[cell + 1]; ++edge) {
                if (edge_weight[edge] <= 0) continue;
                const int neighbor = inner[edge];
                patches[output_chart[cell] - 1].cells.push_back(neighbor);
                patches[output_chart[neighbor] - 1].cells.push_back(cell);
            }
        }
        std::vector<int> incidence_count(n_cells, 0);
        for (auto& patch : patches) {
            std::sort(patch.cells.begin(), patch.cells.end());
            patch.cells.erase(std::unique(patch.cells.begin(), patch.cells.end()), patch.cells.end());
            for (int cell : patch.cells) ++incidence_count[cell];
        }
        // w_ci=1/m_i gives every cell total mass one. Q_c is an orthonormal
        // basis of sqrt(W_c)[1,X_c]. U_c=R_c' sqrt(W_c) Q_c, so 0 <= UU' <= I.
        // All numerically identifiable feature directions are retained; PC1 is
        // recorded only as a diagnostic and does not determine the clock.
        std::vector<Triplet> basis_entries;
        int n_basis = 0;
        for (int chart = 0; chart < n_charts; ++chart) {
            auto& patch = patches[chart];
            Eigen::MatrixXd local(patch.cells.size(), n_dims);
            Eigen::VectorXd weight(patch.cells.size());
            for (int row = 0; row < local.rows(); ++row) {
                local.row(row) = observed.row(patch.cells[row]);
                weight[row] = 1.0 / incidence_count[patch.cells[row]];
            }
            const double mass = weight.sum();
            const Eigen::RowVectorXd center = weight.transpose() * local / mass;
            local.rowwise() -= center;
            const Eigen::MatrixXd scatter = local.transpose() * weight.asDiagonal() * local;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(scatter);
            if (eigen.info() != Eigen::Success) throw std::runtime_error("Local expression decomposition failed");
            const double leading = eigen.eigenvalues()[n_dims - 1];
            const double tolerance = std::numeric_limits<double>::epsilon() * n_dims * std::max(leading, 0.0);
            patch.offset = n_basis;
            patch.rank = 1;
            for (int row = 0; row < local.rows(); ++row) {
                basis_entries.emplace_back(patch.cells[row], n_basis, weight[row] / std::sqrt(mass));
            }
            ++n_basis;
            for (int dim = 0; dim < n_dims; ++dim) {
                const double value = eigen.eigenvalues()[dim];
                if (value <= tolerance) continue;
                const Eigen::VectorXd coordinate = local * eigen.eigenvectors().col(dim) / std::sqrt(value);
                for (int row = 0; row < local.rows(); ++row) {
                    basis_entries.emplace_back(patch.cells[row], n_basis, weight[row] * coordinate[row]);
                }
                ++n_basis;
                ++patch.rank;
            }
            patch.variance_fraction = leading > 0 ? leading / eigen.eigenvalues().sum() : 0;
            patch.gap_fraction = leading > 0 && n_dims > 1 ?
                (leading - eigen.eigenvalues()[n_dims - 2]) / leading : 0;
            output_component[chart] = chart_component[chart] + 1;
        }
        Sparse basis(n_cells, n_basis);
        basis.setFromTriplets(basis_entries.begin(), basis_entries.end());
        basis.makeCompressed();
        basis_entries.clear();
        basis_entries.shrink_to_fit();
        const Eigen::MatrixXd gram = Eigen::MatrixXd(basis.transpose() * basis);
        const Eigen::MatrixXd reduced = 2 * Eigen::MatrixXd::Identity(n_basis, n_basis) - gram;
        Eigen::LDLT<Eigen::MatrixXd> factor(reduced);
        if (factor.info() != Eigen::Success || !factor.isPositive()) {
            throw std::runtime_error("Anchored expression system is not positive definite");
        }
        // B t <= max(B g,0): preserve supported mean order, but do not flatten
        // an existing reversal merely because a geometric edge points backward.
        // This is the smallest nonnegative bound making the anchor feasible.
        // The input forest and its separate roots are retained.
        std::vector<Triplet> order_entries;
        std::vector<int> edge_index(n_states, -1);
        int n_edges = 0;
        for (int state = 0; state < n_states; ++state) {
            if (state_size[state] == 0) throw std::runtime_error("Reference state has no cells");
            int node = state;
            for (int depth = 0; reference_parent[node] > 0; ++depth) {
                node = reference_parent[node] - 1;
                if (node < 0 || node >= n_states || depth >= n_states) {
                    throw std::runtime_error("Reference ancestry is not a forest");
                }
            }
            if (reference_parent[state] > 0) edge_index[state] = n_edges++;
        }
        for (int cell = 0; cell < n_cells; ++cell) {
            const int state = input_state[cell] - 1;
            if (edge_index[state] >= 0) {
                order_entries.emplace_back(edge_index[state], cell, -1.0 / state_size[state]);
            }
            for (int child = 0; child < n_states; ++child) {
                if (reference_parent[child] == state + 1) {
                    order_entries.emplace_back(edge_index[child], cell, 1.0 / state_size[state]);
                }
            }
        }
        Sparse order(n_edges, n_cells);
        order.setFromTriplets(order_entries.begin(), order_entries.end());
        const Eigen::VectorXd order_limit = (order * anchor).cwiseMax(0.0);
        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(n_cells, n_edges + 1);
        rhs.col(0) = anchor;
        rhs.rightCols(n_edges) = Eigen::MatrixXd(order.transpose());
        // Woodbury: H^{-1}Y = (Y + U (2I-U'U)^{-1} U'Y)/2, H=2I-UU'.
        // One reduced factorization serves the coordinate and all order constraints.
        const Eigen::MatrixXd fitted = 0.5 * (rhs + basis * factor.solve(basis.transpose() * rhs));
        Eigen::VectorXd time = fitted.col(0);
        Eigen::VectorXd multiplier = Eigen::VectorXd::Zero(n_edges);
        if (n_edges > 0) {
            const Eigen::MatrixXd dual = order * fitted.rightCols(n_edges);
            const Eigen::VectorXd violation = order * time - order_limit;
            std::vector<bool> active(n_edges, false);
            int iterations = 0;
            while (true) {
                const Eigen::VectorXd gradient = violation - dual * multiplier;
                const double tolerance = 64 * std::numeric_limits<double>::epsilon() *
                    (1 + violation.norm() + dual.norm() * multiplier.norm());
                int entering = -1;
                double largest = tolerance;
                for (int edge = 0; edge < n_edges; ++edge) {
                    if (!active[edge] && gradient[edge] > largest) {
                        entering = edge;
                        largest = gradient[edge];
                    }
                }
                if (entering < 0) break;
                active[entering] = true;
                while (true) {
                    if (++iterations > 10 * n_edges * n_edges + 1) {
                        throw std::runtime_error("State-order active-set solve did not converge");
                    }
                    std::vector<int> indices;
                    for (int edge = 0; edge < n_edges; ++edge) if (active[edge]) indices.push_back(edge);
                    Eigen::MatrixXd local(indices.size(), indices.size());
                    Eigen::VectorXd target(indices.size());
                    for (int row = 0; row < local.rows(); ++row) {
                        target[row] = violation[indices[row]];
                        for (int col = 0; col < local.cols(); ++col) local(row, col) = dual(indices[row], indices[col]);
                    }
                    const Eigen::VectorXd solution = local.ldlt().solve(target);
                    Eigen::VectorXd candidate = Eigen::VectorXd::Zero(n_edges);
                    double step = 1;
                    int leaving = -1;
                    for (int row = 0; row < solution.size(); ++row) {
                        const int edge = indices[row];
                        candidate[edge] = solution[row];
                        if (candidate[edge] <= 0) {
                            const double ratio = multiplier[edge] / (multiplier[edge] - candidate[edge]);
                            if (ratio <= step) { step = ratio; leaving = edge; }
                        }
                    }
                    if (leaving < 0) { multiplier = candidate; break; }
                    multiplier += step * (candidate - multiplier);
                    multiplier[leaving] = 0;
                    active[leaving] = false;
                }
            }
            time -= fitted.rightCols(n_edges) * multiplier;
        }
        const Eigen::VectorXd coordinates = basis.transpose() * time;
        const Eigen::VectorXd stationarity = 2 * time - basis * coordinates - anchor + order.transpose() * multiplier;
        const Eigen::VectorXd order_difference = order * time - order_limit;
        const double order_error = n_edges > 0 ? std::max(0.0, order_difference.maxCoeff()) : 0;
        const double complementarity = n_edges > 0 ? (multiplier.array() * order_difference.array()).abs().maxCoeff() : 0;
        const double normal_error = stationarity.norm() / (1 + anchor.norm() + time.norm());
        if (!time.allFinite() || normal_error > 1e-10 || order_error > 1e-10 || complementarity > 1e-8) {
            throw std::runtime_error("Anchored clock did not satisfy its numerical optimality conditions");
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(reduced, Eigen::EigenvaluesOnly);
        output_diagnostics[0] = (time - anchor).squaredNorm();
        output_diagnostics[1] = std::max(0.0, time.squaredNorm() - coordinates.squaredNorm());
        output_diagnostics[2] = normal_error;
        output_diagnostics[3] = order_error;
        output_diagnostics[4] = complementarity;
        output_diagnostics[5] = spectrum.eigenvalues()[0];
        output_diagnostics[6] = spectrum.eigenvalues()[n_basis - 1];
        output_diagnostics[7] = time.minCoeff();
        output_diagnostics[8] = time.maxCoeff();
        output_diagnostics[9] = (multiplier.array() > 0).count();
        const double range = time.maxCoeff() - time.minCoeff();
        if (!(range > 0)) throw std::runtime_error("Anchored clock has no range");
        // One global affine display scale preserves relative component phases.
        Eigen::Map<Eigen::VectorXd> result(output_time, n_cells);
        result = (time.array() - time.minCoeff()) / range;
        for (int chart = 0; chart < n_charts; ++chart) {
            const auto& patch = patches[chart];
            double anchor_mean = 0;
            double time_mean = 0;
            double local_energy = 0;
            for (int cell : patch.owned) { anchor_mean += anchor[cell]; time_mean += result[cell]; }
            for (int cell : patch.cells) local_energy += time[cell] * time[cell] / incidence_count[cell];
            local_energy -= coordinates.segment(patch.offset, patch.rank).squaredNorm();
            double* info = output_info + chart_info_columns * chart;
            info[0] = patch.owned.size();
            info[1] = patch.cells.size();
            info[2] = patch.rank;
            info[3] = patch.variance_fraction;
            info[4] = patch.gap_fraction;
            info[5] = anchor_mean / patch.owned.size();
            info[6] = time_mean / patch.owned.size();
            info[7] = std::max(0.0, local_energy);
        }
        return 0;
    } catch (const std::exception& failure) {
        std::snprintf(error, error_size, "%s", failure.what());
        return 1;
    } catch (...) {
        std::snprintf(error, error_size, "%s", "State-clock native solver failed");
        return 1;
    }
}
