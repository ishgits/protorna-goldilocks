"""
analyze.py — duplex-stability readouts and relative ranking.

Design principle: the *ranking logic* is pure NumPy and depends on nothing exotic, so it can
be unit-tested anywhere (see src/validate_analysis.py). Reading real trajectories (MDAnalysis)
is an optional adapter layered on top. The scientific observable is RELATIVE rank-ordering,
never absolute dG.

Core observables per duplex trajectory:
  * base-pair persistence : fraction of frames each Watson-Crick pair stays "closed"
                            (heavy-atom H-bond donor/acceptor distance below a cutoff)
  * RMSD stability        : 1 / (1 + late-window backbone-RMSD spread) -- low, steady RMSD = stable
  * stability_score       : a bounded composite used only to RANK duplexes against each other
"""
from __future__ import annotations
import json
from dataclasses import dataclass, asdict
from typing import Sequence
import numpy as np

WC_HBOND_CUTOFF_ANG = 3.5  # heavy-atom donor..acceptor cutoff for a "closed" pair


# ----------------------------------------------------------------------------- core (pure numpy)
def base_pair_persistence(bp_distances: np.ndarray, cutoff: float = WC_HBOND_CUTOFF_ANG) -> np.ndarray:
    """bp_distances: (n_frames, n_pairs) H-bond distances (Angstrom).
    Returns per-pair fraction of frames the pair is closed."""
    bp_distances = np.asarray(bp_distances, dtype=float)
    closed = bp_distances <= cutoff
    return closed.mean(axis=0)


def rmsd_stability(rmsd_series: np.ndarray, late_fraction: float = 0.5) -> float:
    """rmsd_series: (n_frames,) backbone RMSD (Angstrom) vs. the starting A-form helix.
    Uses the late window (post-equilibration). Higher = more conformationally stable."""
    r = np.asarray(rmsd_series, dtype=float)
    late = r[int(len(r) * (1 - late_fraction)):]
    spread = late.mean() + late.std()          # penalise both drift and fluctuation
    return 1.0 / (1.0 + spread)


def stability_score(bp_distances: np.ndarray, rmsd_series: np.ndarray,
                    cutoff: float = WC_HBOND_CUTOFF_ANG) -> float:
    """Bounded composite in (0, 1]. For RANKING duplexes against one another only."""
    persistence = float(np.mean(base_pair_persistence(bp_distances, cutoff)))
    stab = rmsd_stability(rmsd_series)
    return float(persistence * stab)


@dataclass
class DuplexResult:
    name: str
    mean_bp_persistence: float
    rmsd_stability: float
    stability_score: float

    @classmethod
    def from_traj(cls, name: str, bp_distances: np.ndarray, rmsd_series: np.ndarray) -> "DuplexResult":
        return cls(
            name=name,
            mean_bp_persistence=float(np.mean(base_pair_persistence(bp_distances))),
            rmsd_stability=float(rmsd_stability(rmsd_series)),
            stability_score=stability_score(bp_distances, rmsd_series),
        )


def rank_duplexes(results: Sequence[DuplexResult]) -> list[DuplexResult]:
    """Most stable first."""
    return sorted(results, key=lambda r: r.stability_score, reverse=True)


def ordering_matches(ranked: Sequence[DuplexResult], expected_names_most_to_least: Sequence[str]) -> bool:
    """Validation gate: did we reproduce the KNOWN relative ordering?"""
    got = [r.name for r in ranked]
    return got == list(expected_names_most_to_least)


def report(results: Sequence[DuplexResult], expected: Sequence[str] | None = None) -> dict:
    ranked = rank_duplexes(results)
    out = {"ranking_most_to_least_stable": [asdict(r) for r in ranked]}
    if expected is not None:
        out["expected_ordering"] = list(expected)
        out["ordering_reproduced"] = ordering_matches(ranked, expected)
    return out


# ------------------------------------------------------------------- optional MDAnalysis adapter
def load_trajectory_observables(topology: str, trajectory: str, pair_atom_indices):
    """Compute (bp_distances, rmsd_series) from real MD output. Requires MDAnalysis.

    pair_atom_indices: list of (donor_idx, acceptor_idx) atom-index tuples defining the WC H-bonds.
    Kept import-guarded so the pure ranking logic above stays testable without MDAnalysis.
    """
    import MDAnalysis as mda
    from MDAnalysis.analysis.rms import RMSD

    u = mda.Universe(topology, trajectory)
    backbone = u.select_atoms("nucleicbackbone")
    rmsd = RMSD(backbone, backbone, ref_frame=0).run().results.rmsd[:, 2]

    pairs = np.array(pair_atom_indices)
    dists = np.empty((len(u.trajectory), len(pairs)))
    for i, _ in enumerate(u.trajectory):
        pos = u.atoms.positions
        d = pos[pairs[:, 0]] - pos[pairs[:, 1]]
        dists[i] = np.linalg.norm(d, axis=1)
    return dists, rmsd


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Rank duplexes by relative MD stability.")
    p.add_argument("--work", required=True, help="dir of per-duplex .npz files (bp_distances, rmsd)")
    p.add_argument("--report", required=True, help="output JSON path")
    args = p.parse_args()

    import glob, os
    results = []
    for f in sorted(glob.glob(os.path.join(args.work, "*.npz"))):
        d = np.load(f)
        name = os.path.splitext(os.path.basename(f))[0]
        results.append(DuplexResult.from_traj(name, d["bp_distances"], d["rmsd"]))
    rep = report(results)
    with open(args.report, "w") as fh:
        json.dump(rep, fh, indent=2)
    print(json.dumps(rep, indent=2))
