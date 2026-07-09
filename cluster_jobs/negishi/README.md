# Negishi GROMACS submission scripts

Use these from the repository root after `gromacs/build_system.sh` has created `work/*.gro` and `work/topol_*.top`.

Key SLURM settings:
- account: `eapsdept`
- partition: `gpu`
- QoS: `normal`
- CPUs per task: `16`
- GPUs per node: `1`
- GROMACS module: `ml rocmcontainers gromacs`
- production walltime: `24:00:00`

First check the module:

```bash
bash cluster_jobs/negishi/00_check_gromacs_module.sh
```

Manual submit sequence:

```bash
mkdir -p logs reports work
sbatch cluster_jobs/negishi/slurm_em_nvt.sh
# after all array tasks complete successfully:
sbatch cluster_jobs/negishi/slurm_npt_eq.sh
# after all array tasks complete successfully:
sbatch cluster_jobs/negishi/slurm_npt_prod.sh
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
