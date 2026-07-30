#!/usr/bin/env bash
set -euo pipefail
export HOME=/home/bio
CONDA=/home/bio/miniforge3

echo "=== STEP 1: cleanup broken colon-path install ==="
rm -rf "/home/bio/C:UsersLenovo/miniforge3" || true

echo "=== STEP 2: ensure installer present ==="
if [ ! -f /home/bio/miniforge.sh ]; then
  wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /home/bio/miniforge.sh
fi

echo "=== STEP 3: (re)install base conda to clean path ==="
if [ ! -x "$CONDA/bin/conda" ]; then
  bash /home/bio/miniforge.sh -b -p "$CONDA" -f
fi
"$CONDA/bin/conda" --version
echo "REINSTALL_DONE"

echo "=== STEP 4: create cr2 env (cellrank + slepc4py + petsc4py) ==="
"$CONDA/bin/mamba" create -y -n cr2 -c conda-forge \
  python=3.11 cellrank petsc4py slepc4py h5py scanpy anndata
echo "ENV_DONE"

echo "=== STEP 5: verify imports (slepc must be importable) ==="
"$CONDA/envs/cr2/bin/python" - <<'PY'
import cellrank, scanpy, anndata, numpy
import petsc4py, slepc4py
print("cellrank", cellrank.__version__)
print("scanpy", scanpy.__version__)
print("petsc4py", petsc4py.__version__)
print("slepc4py", slepc4py.__version__)
print("IMPORT_OK")
PY
echo "ALL_DONE"
