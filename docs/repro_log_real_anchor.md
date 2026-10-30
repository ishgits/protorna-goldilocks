# Reproducibility dry-run — `real_anchor/`

**Date:** 2026-07-09 · **Verdict: PASS** (energy ordering
reproduces; energies match committed values within stated error). One expected-and-informative
difference in the *geometric* readout is documented below — it strengthens, rather than undermines,
the central honesty claim.

---

## 1. What was run

A clean-room regeneration of the headline artifact, following the reproduce recipe in
`docs/METHODS.md` §1B. The committed `real_anchor/work/` was left untouched; all output was
written to a fresh `real_anchor/repro_work/`.

```bash
# fresh conda env: python 3.11, conda-forge
#   openmm pdbfixer mdtraj numpy scipy matplotlib pyyaml
cd real_anchor
for a in gc mixed au; do
  python build_anchor.py --scaffold scaffold_1rna.pdb --anchor $a --out repro_work --ns 0.02 --temp 340
done
# analyze with WORK pointed at repro_work/ (module import; WORK="repro_work")
python -c "import analyze_real as AR; AR.WORK='repro_work'; AR.main()"
```

Each anchor: central 6 bp of PDB **1RNA** (A-form), base-mutated on one identical backbone
(gc = GCGCGC, mixed = GCAUGC, au = AUAUAU), OL3 + GBn2 implicit solvent, **~20 ps (4999 × 4 fs) at
340 K on CPU**, 66 analysis frames. Interstrand interaction energy is a single-trajectory MM-GBSA-style
readout: `E(complex) − E(strandA) − E(strandB)` with `amber14-all.xml` + `implicit/gbn2.xml`,
`NoCutoff`.

## 2. Environment (auditable record)

- platform: macOS-15.1.1-arm64-arm-64bit (arm64)
- python: 3.11.15
- openmm: 8.5.2.dev-36a30cb
- mdtraj: 1.11.1
- numpy: 2.4.6
- scipy: 1.17.1
- pdbfixer: 1.11 (conda-forge)
- compute platform used: **CPU** (Reference/CPU/OpenCL available; no CUDA)

> The committed numbers were produced on a comparable OpenMM 8.5 CPU build. arm64 vs the original
> machine is exactly the expected case: small energy drift is expected and fine; a *changed
> ordering* would be the real finding. The ordering did **not** change.

## 3. Result — interstrand interaction energy (the readout that discriminates)

| anchor | regenerated Eint (kcal/mol) | committed Eint (kcal/mol) | Δ | \|Δ\|/SD |
|---|---|---|---|---|
| gc    | **−69.63 ± 4.23** | −69.81 ± 4.01 | +0.19 | 0.05 |
| mixed | **−58.34 ± 3.46** | −57.82 ± 3.35 | −0.52 | 0.15 |
| au    | **−40.14 ± 3.54** | −40.81 ± 2.72 | +0.67 | 0.25 |

- **Ordering reproduced: GC > mixed > AU** — identical to committed, identical to `expected_ordering`.
- Every energy matches the committed value to **well under 1 kcal/mol**, and every delta is **≤ 0.25×**
  the committed standard deviation. This is agreement, not coincidence.
- The gaps between anchors (~11 and ~18 kcal/mol) remain large relative to the per-frame fluctuations (GC–mixed ≈ 3.2, mixed–AU > 5 pooled SD), so the rank separation survives
  the run-to-run noise. **PASS.**

## 4. The geometric `stability_score` — expected non-reproduction, and why it matters

The pipeline computes a second, purely geometric readout (`stability_score = bp_persistence ×
rmsd_stability`). It is **not** the discriminating observable, and this dry-run makes that concrete:

| anchor | stability_score (this run) | stability_score (committed) |
|---|---|---|
| gc    | 0.4927 | 0.4447 |
| mixed | 0.4701 | 0.4596 |
| au    | 0.4418 | 0.4419 |
| **spread (max−min)** | **0.0508** | **0.0177** |
| **`ordering_reproduced` flag** | **True** | **False** |

The `ordering_reproduced` flag — which is driven by the geometric score, not the energy — **flipped
between the two runs**. In this run the geometric score happened to fall GC > mixed > AU; in the
committed run it did not (mixed edged out gc). The between-anchor spread is only **0.018–0.051** on a
~0.44–0.49 baseline — the anchors are geometrically indistinguishable at 20 ps, all base pairs stay
closed, and the resulting rank is dominated by thermal noise in a short trajectory.

**This is the central finding, reproduced live:** at short CPU sampling the **energy** readout
robustly separates the anchors (deltas < 0.25 SD, gaps ≈ 3.2 / > 5 pooled SD) while the **geometric** readout does
not (rank flips between statistically-identical runs). The honest statement — *"recovers GC > mixed >
AU via the MD interaction-energy readout; the geometric score needs multi-ns GPU sampling to resolve
the ordering"* — is exactly right, and the committed JSON's `ordering_reproduced: false` should be
read as a property of the geometric score at 20 ps, **not** a failure of the result.

> **Do not** "fix" this by shipping the run where the geometric flag came out True — that would be
> cherry-picking noise. The correct artifact is the one that reports the energy ordering as the result
> and labels the geometric score as unresolved at short sampling.

## 5. Files produced

- `real_anchor/repro_work/{gc,mixed,au}.dcd` — regenerated trajectories
- `real_anchor/repro_work/{gc,mixed,au}_top.pdb`, `*_complete.pdb` — topologies / built structures
- `real_anchor/repro_work/anchor_ordering.json` — regenerated summary (energies above)

These clean-room reproduction files are written under `real_anchor/repro_work/`, are regenerated on demand, and are **not committed** by default (they are `.gitignore`d). A reviewer regenerates them by running the six lines in §1.

## 6. Bottom line

The headline artifact **regenerates end-to-end from a clean checkout**: build → MD → analyze recovers
GC > mixed > AU on the interstrand-energy readout, with energies matching the committed values within
error. A reviewer opening the repo can run six lines and watch the result fall out — it is
demonstrated, not asserted.
