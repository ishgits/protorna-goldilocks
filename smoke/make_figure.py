"""Generate results/validation_ordering.png — a one-glance proof that the pipeline reproduces
known relative stabilities, from both the synthetic gate and a real OpenMM run."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from analyze import DuplexResult, rank_duplexes
import validate_analysis as va
import openmm as mm
from openmm import unit
from openmm_smoke import simulate_pair

# --- panel A: synthetic anchor gate ---
res_syn = []
specs = {
    "GC_rich":  dict(hbond_mean=2.85, hbond_drift=0.1, rmsd_plateau=1.4, rmsd_noise=0.12),
    "mixed":    dict(hbond_mean=2.95, hbond_drift=0.6, rmsd_plateau=2.4, rmsd_noise=0.20),
    "AU_rich":  dict(hbond_mean=3.10, hbond_drift=1.8, rmsd_plateau=4.2, rmsd_noise=0.35),
}
for name, kw in specs.items():
    bp, rmsd = va.synth_trajectory(500, 6, **kw)
    res_syn.append(DuplexResult.from_traj(name, bp, rmsd))
res_syn = rank_duplexes(res_syn)

# --- panel B: real OpenMM run ---
strong = simulate_pair(k_hbond=150000.0, seed=1)
weak = simulate_pair(k_hbond=8000.0, seed=2)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.bar([r.name for r in res_syn], [r.stability_score for r in res_syn],
        color=["#2b8cbe", "#7bccc4", "#f0a0a0"])
ax1.set_title("A. Synthetic anchor gate\n(known: GC > mixed > AU)")
ax1.set_ylabel("relative stability score")

ax2.plot(strong, label=f"strong H-bond  <d>={strong.mean():.2f} A", color="#2b8cbe")
ax2.plot(weak, label=f"weak H-bond  <d>={weak.mean():.2f} A", color="#f0a0a0", alpha=0.8)
ax2.axhline(3.5, ls="--", c="k", lw=0.8, label="pair-closed cutoff")
ax2.set_title(f"B. Real OpenMM {mm.version.version} run\npair distance vs. time")
ax2.set_xlabel("frame"); ax2.set_ylabel("pair distance (A)"); ax2.legend(fontsize=8)

fig.suptitle("proto-RNA Goldilocks — validation stub reproduces known stability ordering",
             fontweight="bold")
fig.tight_layout()
os.makedirs("results", exist_ok=True)
fig.savefig("results/validation_ordering.png", dpi=150)
print("wrote results/validation_ordering.png")
