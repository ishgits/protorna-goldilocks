#!/bin/bash
#SBATCH --job-name=protorna-md
#SBATCH --account=YOUR_RCAC_ALLOCATION      # <-- set Purdue RCAC allocation
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-2                          # one duplex per array task
#SBATCH --output=logs/md_%A_%a.out

# Purdue RCAC SLURM template. The pipeline dispatches jobs like this via sbatch and pulls
# outputs back into the session, keeping large trajectories on the cluster.

module load anaconda cuda                    # adjust to RCAC module names
source activate protorna

SYSTEMS=(gc_anchor mixed_anchor au_anchor)   # zp_hachimoji excluded: not yet built/parameterized
SYS=${SYSTEMS[$SLURM_ARRAY_TASK_ID]}

python src/run_md.py \
  --system work/${SYS} \
  --ns 50 \
  --platform CUDA \
  --out work/md

echo "finished ${SYS}"
