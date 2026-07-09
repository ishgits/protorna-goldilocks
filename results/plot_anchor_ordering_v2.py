#!/usr/bin/env python3
"""Two-panel anchor-ordering figure.

Left  (A): MD interstrand interaction energy per anchor (mean +/- SD, 66 frames),
           from real_anchor/work/anchor_ordering.json.
Right (B): Turner nearest-neighbor dG37 per anchor, from data/turner_benchmarks.csv.

Both readouts recover the SAME rank order GC > mixed > AU. The figure shows
*relative* ordering only -- NOT absolute Tm or absolute dG. The geometric
stability_score is deliberately NOT plotted: at ~20 ps it does not resolve the
ordering (see docs/repro_log_real_anchor.md).

Reproduce:
    conda env: python 3.11, numpy matplotlib  (see docs/repro_log_real_anchor.md for versions)
    python results/plot_anchor_ordering_v2.py       # run from repo root
Outputs: results/anchor_ordering_v2.png (300 dpi) + results/anchor_ordering_v2.svg
"""
import json, csv, os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---- data ----
J = json.load(open(os.path.join(REPO, "real_anchor/work/anchor_ordering.json")))
pa = J["per_anchor"]
order = ["gc", "mixed", "au"]
md_e  = {a: pa[a]["interstrand_interaction_kcal"] for a in order}
md_sd = {a: pa[a]["interstrand_interaction_sd"]  for a in order}
turner = {}
for row in csv.DictReader(open(os.path.join(REPO, "data/turner_benchmarks.csv"))):
    key = {"gc_anchor": "gc", "mixed_anchor": "mixed", "au_anchor": "au"}.get(row["name"])
    if key:
        turner[key] = float(row["nn_dG37_kcal_mol"])

# ---- style ----
mpl.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 120, "savefig.bbox": "tight",
})
GREY = "#6E6E6E"
COL = {"gc": "#0072B2", "mixed": "#009E73", "au": "#D55E00"}  # Okabe-Ito, CVD-safe
xlab = ["GC\n(GCGCGC)", "mixed\n(GCAUGC)", "AU\n(AUAUAU)"]
x = np.arange(3)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.2, 3.6))

# Panel A -- MD interstrand energy
for i, a in enumerate(order):
    axL.plot([i, i], [0, md_e[a]], color=COL[a], lw=1.2, alpha=0.5, zorder=1)
    axL.errorbar(i, md_e[a], yerr=md_sd[a], fmt='o', ms=9, color=COL[a],
                 capsize=4, elinewidth=1.4, mec='white', mew=0.8, zorder=3)
    axL.annotate(f"{md_e[a]:.0f}", (i, md_e[a]), textcoords="offset points",
                 xytext=(10, 0), va='center', fontsize=7, color=COL[a], fontweight='bold')
axL.axhline(0, color=GREY, lw=0.8, ls='--')
axL.set_xticks(x); axL.set_xticklabels(xlab)
axL.set_ylabel("MD interstrand interaction energy\n(kcal/mol, mean \u00b1 SD)")
axL.set_title("A.  MD readout (this work)", loc='left')
axL.set_ylim(-82, 6); axL.margins(x=0.15)
axL.text(0.02, 0.03, "more negative = more stable", transform=axL.transAxes,
         fontsize=6.5, style='italic', color=GREY, ha='left', va='bottom')

# Panel B -- Turner NN dG37
for i, a in enumerate(order):
    axR.plot([i, i], [0, turner[a]], color=COL[a], lw=1.2, alpha=0.5, zorder=1)
    axR.plot(i, turner[a], 'o', ms=9, color=COL[a], mec='white', mew=0.8, zorder=3)
    axR.annotate(f"{turner[a]:.2f}", (i, turner[a]), textcoords="offset points",
                 xytext=(10, 0), va='center', fontsize=7, color=COL[a], fontweight='bold')
axR.axhline(0, color=GREY, lw=0.8, ls='--')
axR.set_xticks(x); axR.set_xticklabels(xlab)
axR.set_ylabel("Turner NN $\\Delta G^{\\circ}_{37}$ (kcal/mol)")
axR.set_title("B.  Turner NN thermodynamics (independent)", loc='left')
axR.set_ylim(-12.5, 1.8); axR.margins(x=0.15)
axR.text(0.02, 0.03, "more negative = more stable", transform=axR.transAxes,
         fontsize=6.5, style='italic', color=GREY, ha='left', va='bottom')

fig.suptitle("Both readouts recover the same rank order:  GC > mixed > AU",
             fontsize=9, y=1.00)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(REPO, "results/anchor_ordering_v2.png"), dpi=300)
fig.savefig(os.path.join(REPO, "results/anchor_ordering_v2.svg"))
print("wrote results/anchor_ordering_v2.{png,svg}")
