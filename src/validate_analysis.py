"""
validate_analysis.py — self-contained proof that the ranking logic reproduces KNOWN ordering.

No GPU, no MD engine, no exotic deps (NumPy only). Runs in seconds anywhere.

We synthesise three controlled duplex "trajectories" whose relative stabilities are known by
construction — a stable GC-rich anchor, an intermediate mixed duplex, and an unstable AU-rich
anchor that frays — and assert the analysis pipeline ranks them in the correct order. This is
the epistemic gate the whole screen depends on: if the readout can't rank the certain cases,
it can't be trusted on novel bases.
"""
from __future__ import annotations
import numpy as np
from analyze import DuplexResult, rank_duplexes, ordering_matches, report
import json

RNG = np.random.default_rng(20260901)  # deterministic


def synth_trajectory(n_frames: int, n_pairs: int, *, hbond_mean: float, hbond_drift: float,
                     rmsd_plateau: float, rmsd_noise: float):
    """Fabricate (bp_distances, rmsd_series) with a controlled stability level.

    hbond_mean   : equilibrium H-bond distance (A). ~2.9 = tightly paired.
    hbond_drift  : linear opening over the run (A). >0 = fraying (unstable).
    rmsd_plateau : late-window backbone RMSD (A). higher = more distorted.
    """
    t = np.linspace(0, 1, n_frames)
    # base-pair H-bond distances: mean + fraying drift + per-pair thermal noise
    drift = np.outer(t, np.full(n_pairs, hbond_drift))
    noise = RNG.normal(0.0, 0.15, size=(n_frames, n_pairs))
    bp = hbond_mean + drift + noise
    # backbone RMSD: rises then plateaus, with noise
    rmsd = rmsd_plateau * (1 - np.exp(-5 * t)) + RNG.normal(0.0, rmsd_noise, size=n_frames)
    return np.abs(bp), np.abs(rmsd)


def main() -> int:
    specs = {
        # name          hbond_mean  drift  rmsd_plateau  rmsd_noise   (known stability)
        "GC_rich_anchor": dict(hbond_mean=2.85, hbond_drift=0.1, rmsd_plateau=1.4, rmsd_noise=0.12),
        "mixed_duplex":   dict(hbond_mean=2.95, hbond_drift=0.6, rmsd_plateau=2.4, rmsd_noise=0.20),
        "AU_rich_anchor": dict(hbond_mean=3.10, hbond_drift=1.8, rmsd_plateau=4.2, rmsd_noise=0.35),
    }
    expected = ["GC_rich_anchor", "mixed_duplex", "AU_rich_anchor"]  # most -> least stable

    results = []
    for name, kw in specs.items():
        bp, rmsd = synth_trajectory(n_frames=500, n_pairs=6, **kw)
        results.append(DuplexResult.from_traj(name, bp, rmsd))

    ranked = rank_duplexes(results)
    rep = report(results, expected=expected)
    print(json.dumps(rep, indent=2))

    ok = ordering_matches(ranked, expected)
    print("\n" + ("PASS" if ok else "FAIL") +
          f": pipeline {'reproduced' if ok else 'did NOT reproduce'} the known stability ordering "
          f"({' > '.join(expected)}).")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
