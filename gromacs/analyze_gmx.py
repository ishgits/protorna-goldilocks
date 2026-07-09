"""
analyze_gmx.py - turn GROMACS .xvg outputs into the SAME stability ranking used everywhere else.

Reads, per duplex:
    work/rmsd_<id>.xvg   (gmx rms  -> backbone RMSD vs. equilibrated start, nm)
    work/hbnum_<id>.xvg  (gmx hbond -num -> inter-strand H-bond count over time)

and reuses src/analyze.py's tested logic (rmsd_stability + a persistence proxy) to rank the
duplexes and check the ordering against configs/duplexes.yaml. The scientific readout is
RELATIVE rank-ordering, identical in spirit to the OpenMM path -- only the data source differs.
"""
from __future__ import annotations
import os, sys, glob, json, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from analyze import rmsd_stability, DuplexResult, rank_duplexes, ordering_matches  # noqa: E402


def read_xvg(path: str) -> np.ndarray:
    """Return the y-column (col 1) of a GROMACS .xvg, skipping @/# comment lines."""
    ys = []
    with open(path) as fh:
        for line in fh:
            if line[0] in "#@" or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                ys.append(float(parts[1]))
    return np.asarray(ys, dtype=float)


def build_results(work: str, ids: list[str]) -> list[DuplexResult]:
    # first pass: gather H-bond series to set a common normalisation (relative persistence)
    hb_series, rmsd_series = {}, {}
    for d in ids:
        rmsd_series[d] = read_xvg(os.path.join(work, f"rmsd_{d}.xvg")) * 10.0   # nm -> Angstrom
        hp = os.path.join(work, f"hbnum_{d}.xvg")
        hb_series[d] = read_xvg(hp) if os.path.exists(hp) else None

    all_hb = np.concatenate([v for v in hb_series.values() if v is not None]) \
        if any(v is not None for v in hb_series.values()) else None
    hb_ref = np.percentile(all_hb, 95) if all_hb is not None and all_hb.size else 1.0

    results = []
    for d in ids:
        rstab = rmsd_stability(rmsd_series[d])
        if hb_series[d] is not None and hb_ref > 0:
            persistence = float(np.clip(hb_series[d].mean() / hb_ref, 0.0, 1.0))
        else:
            persistence = 1.0  # RMSD-only fallback
        score = persistence * rstab
        results.append(DuplexResult(name=d, mean_bp_persistence=persistence,
                                    rmsd_stability=float(rstab), stability_score=float(score)))
    return results


def main() -> int:
    p = argparse.ArgumentParser(description="Rank RNA duplexes from GROMACS .xvg outputs.")
    p.add_argument("--work", default="work")
    p.add_argument("--ids", nargs="*", default=None,
                   help="duplex ids; default = auto-detect rmsd_*.xvg in --work")
    p.add_argument("--expected", nargs="*", default=["gc_anchor", "mixed_anchor", "au_anchor"],
                   help="known ordering, most -> least stable (validation gate)")
    p.add_argument("--report", default="reports/gmx_ordering.json")
    args = p.parse_args()

    ids = args.ids or sorted(
        os.path.basename(f)[5:-4] for f in glob.glob(os.path.join(args.work, "rmsd_*.xvg")))
    if not ids:
        print(f"no rmsd_*.xvg found in {args.work}/ -- run production + analysis first.")
        return 2

    results = build_results(args.work, ids)
    ranked = rank_duplexes(results)
    expected = [e for e in args.expected if e in ids]  # only score anchors we actually have
    ok = ordering_matches([r for r in ranked if r.name in expected], expected) if expected else None

    rep = {"ranking_most_to_least_stable": [r.__dict__ for r in ranked],
           "expected_ordering": expected, "ordering_reproduced": ok}
    os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
    with open(args.report, "w") as fh:
        json.dump(rep, fh, indent=2)
    print(json.dumps(rep, indent=2))
    if ok is not None:
        print(("PASS" if ok else "FAIL") + ": GROMACS pipeline "
              + ("reproduced" if ok else "did NOT reproduce") + " the known ordering "
              + " > ".join(expected))
    return 0 if (ok or ok is None) else 1


if __name__ == "__main__":
    raise SystemExit(main())
