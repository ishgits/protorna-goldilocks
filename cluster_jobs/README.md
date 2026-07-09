# Cluster-specific GROMACS job scripts

This folder contains two SLURM variants for the same GROMACS workflow:

- `negishi/` — account `eapsdept`, partition `gpu`, ROCm container GROMACS.
- `gautschi/` — account `pearce21`, partition `smallgpu`, module stack `gcc/11.4.1 openmpi/4.1.6 gromacs/2024.3`.

Run all commands from the repository root so paths like `gromacs/mdp`, `work`, `logs`, and `reports` resolve correctly.

Current scope: three canonical validation anchors only:

```text
gc_anchor
mixed_anchor
au_anchor
```

The `zp_hachimoji` placeholder was intentionally removed from these submission arrays because the current `build_system.sh` does not yet build or parameterize Z/P residues.
