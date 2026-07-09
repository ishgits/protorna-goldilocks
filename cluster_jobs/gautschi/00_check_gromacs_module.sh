#!/bin/bash
set -euo pipefail
module --force purge
module load gcc/11.4.1 openmpi/4.1.6
module load gromacs/2024.3
if command -v gmx >/dev/null 2>&1; then GMX=gmx; elif command -v gmx_mpi >/dev/null 2>&1; then GMX=gmx_mpi; else echo "No GROMACS executable found"; exit 1; fi
echo "Using: $(command -v ${GMX})"
${GMX} --version | head -40
