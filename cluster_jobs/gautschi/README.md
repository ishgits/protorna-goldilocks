# Gautschi GROMACS submission scripts

Use these from the repository root after `gromacs/build_system.sh` has created `work/*.gro` and `work/topol_*.top`.

Key SLURM settings:
- account: `pearce21`
- partition: `smallgpu`
- CPUs per task: `16`
- GPUs per node: `1`
- GROMACS module stack:
  - `gcc/11.4.1`
  - `openmpi/4.1.6`
  - `gromacs/2024.3`
- executable auto-detection: script uses `gmx` if found, otherwise `gmx_mpi`

First check the module:

```bash
bash cluster_jobs/gautschi/00_check_gromacs_module.sh
```

Manual submit sequence:

```bash
mkdir -p logs reports work
sbatch cluster_jobs/gautschi/slurm_em_nvt.sh
# after all array tasks complete successfully:
sbatch cluster_jobs/gautschi/slurm_npt_eq.sh
# after all array tasks complete successfully:
sbatch cluster_jobs/gautschi/slurm_npt_prod.sh
```

Monitor:

```bash
squeue -u $USER
sacct -j <JOBID> --format=JobID,JobName%24,State,Elapsed,MaxRSS,ExitCode
```

Analyze after production:

```bash
python gromacs/analyze_gmx.py --work work
```
