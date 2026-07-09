#!/bin/bash
# ============================================================
# slurm_npt_prod.sh - NPT production (50 ns) for an RNA duplex
# proto-RNA Goldilocks - Purdue RCAC Negishi (account=eapsdept)
# Estimated wall time: ~20-40 h for 50 ns (small duplex + water)
# ============================================================
#SBATCH --job-name=rna_prod
#SBATCH -A eapsdept
#SBATCH --partition=gpu
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/prod_%A_%a.log
#SBATCH --error=logs/prod_%A_%a.err

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

echo "=== [$D] NPT production (50 ns) ==="
# grompp only if the .tpr doesn't exist yet (so re-submits just restart)
if [ ! -f work/prod_${D}.tpr ]; then
  ${GMX} grompp -f ${MDP}/npt_prod.mdp -c work/npt_eq_${D}.gro -t work/npt_eq_${D}.cpt \
             -p ${TOPOL} -o work/prod_${D}.tpr -maxwarn 2
fi
${GMX} mdrun -v -deffnm work/prod_${D} -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK} \
          -gpu_id 0 -nb gpu -pme gpu \
          -cpi work/prod_${D}.cpt        # auto-restart if job is interrupted / re-queued
echo "[$D] Production complete -> work/prod_${D}.xtc"

# --- analysis: backbone RMSD + inter-strand H-bond count (stability observables) ---
echo "=== [$D] analysis ==="
# backbone RMSD vs. the equilibrated start (group: RNA)
echo "RNA RNA" | ${GMX} rms -s work/prod_${D}.tpr -f work/prod_${D}.xtc \
                         -o work/rmsd_${D}.xvg -tu ns
# inter-strand hydrogen bonds over time (base-pairing proxy). Select strand1 then strand2.
# Define index groups strand1/strand2 once with:  ${GMX} make_ndx -f work/prod_${D}.tpr
echo "strand1 strand2" | ${GMX} hbond -s work/prod_${D}.tpr -f work/prod_${D}.xtc \
                         -num work/hbnum_${D}.xvg 2>/dev/null || \
  echo "  (create strand1/strand2 index groups for hbond; see BUILD_RECIPE.md)"
echo "[$D] wrote work/rmsd_${D}.xvg  work/hbnum_${D}.xvg"
