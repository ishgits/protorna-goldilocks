#!/bin/bash
# ============================================================
# slurm_npt_eq.sh - NPT equilibration (1 ns) for an RNA duplex
# proto-RNA Goldilocks - Purdue RCAC Gautschi smallgpu
# ============================================================
#SBATCH --job-name=rna_npt_eq
#SBATCH -A pearce21
#SBATCH --partition=smallgpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/npt_eq_%A_%a.log
#SBATCH --error=logs/npt_eq_%A_%a.err

module --force purge
module load gcc/11.4.1 openmpi/4.1.6
module load gromacs/2024.3

if command -v gmx >/dev/null 2>&1; then
  GMX="gmx"
elif command -v gmx_mpi >/dev/null 2>&1; then
  GMX="gmx_mpi"
else
  echo "ERROR: no GROMACS executable found after loading gromacs/2024.3" >&2
  echo "Try: module spider gromacs/2024.3" >&2
  exit 1
fi
${GMX} --version | head -20
mkdir -p logs reports work


MDP=gromacs/mdp
DUPLEXES=(gc_anchor mixed_anchor au_anchor)
D=${DUPLEXES[$SLURM_ARRAY_TASK_ID]}
TOPOL="work/topol_${D}.top"

echo "=== [$D] NPT equilibration (1 ns, 300 K / 1 bar) ==="
${GMX} grompp -f ${MDP}/npt_eq.mdp -c work/nvt_${D}.gro -r work/nvt_${D}.gro \
           -t work/nvt_${D}.cpt -p ${TOPOL} -o work/npt_eq_${D}.tpr -maxwarn 2
${GMX} mdrun -v -deffnm work/npt_eq_${D} -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK} \
          -gpu_id 0 -nb gpu -pme gpu
echo "[$D] NPT eq complete. Check density:"
echo "  ${GMX} energy -f work/npt_eq_${D}.edr -o work/density_${D}.xvg   (select Density)"
