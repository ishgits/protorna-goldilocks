# real_anchor — one genuinely real anchor ordering (no synthetic input, no cluster)

This is the honest, physics-based replacement for the synthetic gate in
`src/validate_analysis.py`. It runs **real OL3 MD on real A-form RNA coordinates** and feeds
the resulting trajectories through the repo's tested `src/analyze.py`, entirely on CPU with
lightweight compute (OpenMM + implicit GBn2 solvent). No pNAB, no AmberTools, no GPU cluster.

## What it does

1. **Scaffold.** `scaffold_1rna.pdb` is the central 6 base pairs of PDB **1RNA**
   (Dock-Bregeon et al., *J. Mol. Biol.* 1989 — `[U(UA)6A]2`, an A-form RNA duplex). Chain A =
   `AUAUAU`, chain B its antiparallel Watson–Crick complement.
2. **Build the exact config anchors by base mutation** (`build_anchor.py`). Mutating A→G and
   U→C on **both** strands preserves complementarity, so from one identical A-form backbone we
   get the three `configs/duplexes.yaml` anchors:
   `au = AUAUAU`, `gc = GCGCGC`, `mixed = GCAUGC`. PDBFixer rebuilds the base heavy atoms and
   hydrogens; OpenMM's `amber14/RNA.OL3.xml` supplies the force field. This is a **controlled**
   comparison — identical starting helix and identical protocol, only base identity differs.
3. **MD.** Short OL3 / GBn2 implicit-solvent Langevin MD at 340 K (per anchor, one CPU call).
4. **Readouts** (`analyze_real.py`):
   - Watson–Crick pair geometry (purine-N1 ⋯ pyrimidine-N3) → `src/analyze.py`
     base-pair persistence + `stability_score` (the repo's own logic).
   - backbone RMSD (mdtraj).
   - **MM-GBSA-style interstrand interaction energy** `E(duplex) − E(strandA) − E(strandB)` on
     the OL3/GBn2 potential — the "MD-binding-energy" used as the stability proxy.

## Result (see `work/anchor_ordering.json`, `../results/real_anchor_ordering.png`)

| anchor | seq | interstrand interaction (kcal/mol) | Turner NN dG37 (kcal/mol) |
|---|---|---|---|
| gc    | GCGCGC | **−69.8 ± 4.0** | −10.46 |
| mixed | GCAUGC | **−57.8 ± 3.3** | −7.64 |
| au    | AUAUAU | **−40.8 ± 2.7** | −0.54 |

The **energy readout recovers the known ordering GC > mixed > AU** from real dynamics, and it
**rank-agrees with independent Turner nearest-neighbor thermodynamics** (`src/nn_predict.py`).

**Honest caveat (important):** at this *short* sampling the **geometric** readout
(`stability_score` from bp-persistence + RMSD) does **not** resolve the ordering — all pairs
stay closed over 20 ps, so `ordering_reproduced` is `False` for the geometric composite. The
ordering signal here comes from the **energy** observable, which separates 3-H-bond G:C from
2-H-bond A:U with little sampling. Reproducing the ordering *geometrically* needs long
(multi-ns) MD and is the thing the cluster runs are for. This is exactly the gap between the
"reproduces known relative stabilities end-to-end" wording and what short/cheap MD delivers.

## Reproduce

```bash
pip install openmm pdbfixer mdtraj          # CPU is fine
cd real_anchor
for a in gc mixed au; do
  python build_anchor.py --scaffold scaffold_1rna.pdb --anchor $a --out work --ns 0.02 --temp 340
done
python analyze_real.py                       # writes work/anchor_ordering.json
```

Longer runs (`--ns 5 --temp 300`, needs a GPU or patience) let the AU-rich anchor fray and
should make the *geometric* `stability_score` reproduce the ordering too — the recommended
next step before submission.
