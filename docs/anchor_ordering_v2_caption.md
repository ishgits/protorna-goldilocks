# Figure caption — `anchor_ordering_v2`

**The two independent readouts agree on the relative stability rank of three canonical RNA 6-mer
anchors.** **(A)** MD interstrand interaction energy (single-trajectory MM-GBSA-style,
OL3/GBn2 implicit solvent, ~20 ps at 340 K on CPU; mean ± SD over 66 frames) for the three anchors
built by base mutation on one identical A-form backbone (central 6 bp of PDB 1RNA):
GC (GCGCGC) −69.8 ± 4.0 > mixed (GCAUGC) −57.8 ± 3.3 > AU (AUAUAU) −40.8 ± 2.7 kcal/mol; adjacent gaps are large relative to the per-frame fluctuations (GC–mixed ≈ 3.2 pooled SD, mixed–AU > 5 pooled SD). **(B)** Turner nearest-neighbor ΔG°₃₇ for the same sequences (Xia 1998, via NNDB;
computed by `src/nn_predict.py`): −10.46 > −7.64 > −0.54 kcal/mol. Colors identify the anchor and are
shared across panels; on both axes more negative = more stable.

**This figure shows *relative* rank-ordering only — not absolute Tm and not absolute ΔG.** The two
readouts use different scales and are not directly comparable in magnitude; the point is that both
place the anchors in the same order (GC > mixed > AU). The AU anchor is robustly weakest here but, as
an isolated 6-mer (ΔG°₃₇ = −0.54 kcal/mol), is only marginally a duplex — a validation anchor, not a
calibrated low-stability gate. The purely geometric `stability_score` is deliberately **not** plotted:
at ~20 ps it does not resolve the ordering (the run-to-run flag flips; see
`docs/repro_log_real_anchor.md`), so plotting it here would misrepresent it as discriminating.

*Reproduce:* `python results/plot_anchor_ordering_v2.py` (from repo root); inputs are
`real_anchor/work/anchor_ordering.json` and `data/turner_benchmarks.csv`. Colorblind-safe Okabe-Ito
palette.
