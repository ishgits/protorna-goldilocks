#!/bin/bash
# ============================================================
# slurm_em_nvt.sh - EM + NVT equilibration for an RNA duplex
# proto-RNA Goldilocks - adapted from Pearce-group Negishi scripts
# Target: Purdue RCAC Negishi (ROCm/AMD GPU, container GROMACS; account=eapsdept)
# ============================================================
#SBATCH --job-name=rna_em_nvt
#SBATCH -A eapsdept
#SBATCH --partition=gpu
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-2                 # gc_anchor, mixed_anchor, au_anchor
#SBATCH --output=logs/em_nvt_%A_%a.log
#SBATCH --error=logs/em_nvt_%A_%a.err

module --force purge
ml rocmcontainers gromacs

if command -v gmx >/dev/null 2>&1; then
  GMX="gmx"
elif command -v gmx_mpi >/dev/null 2>&1; then
  GMX="gmx_mpi"
else
  echo "ERROR: no GROMACS executable found after loading rocmcontainers gromacs" >&2
  exit 1
fi
${GMX} --version | head -20
mkdir -p logs reports work


MDP=gromacs/mdp
DUPLEXES=(gc_anchor mixed_anchor au_anchor)
D=${DUPLEXES[$SLURM_ARRAY_TASK_ID]}
TOPOL="work/topol_${D}.top"
GRO="work/${D}.gro"                 # solvated + ionised system from build_system.sh
mkdir -p logs

# --- ONE-TIME per system: temperature-coupling index groups -------------------
# The mdp uses tc-grps = "RNA Water_and_ions". If ${GMX} make_ndx does not expose
# those default groups for your OL3/acpype residue names, create them once:
#   ${GMX} make_ndx -f ${GRO} -o work/${D}.ndx
#     > keep the nucleic-acid residues as group "RNA" and the rest as "Water_and_ions"
# then add:  -n work/${D}.ndx   to the grompp calls below.
# ------------------------------------------------------------------------------

echo "=== [$D] Energy minimisation ==="
${GMX} grompp -f ${MDP}/em.mdp -c ${GRO} -p ${TOPOL} -o work/em_${D}.tpr -maxwarn 2
${GMX} mdrun -v -deffnm work/em_${D} -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK}   # CPU: EM is fast
if [ $? -ne 0 ]; then echo "[$D] EM FAILED"; exit 1; fi
echo "[$D] EM done. $(grep -m1 'Fmax' work/em_${D}.log)"

echo "=== [$D] NVT equilibration (500 ps, 300 K) ==="
${GMX} grompp -f ${MDP}/nvt_eq.mdp -c work/em_${D}.gro -r work/em_${D}.gro \
           -p ${TOPOL} -o work/nvt_${D}.tpr -maxwarn 2
${GMX} mdrun -v -deffnm work/nvt_${D} -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK} \
          -gpu_id 0 -nb gpu -pme gpu
echo "[$D] NVT complete -> work/nvt_${D}.gro"
