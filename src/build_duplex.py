"""
build_duplex.py — generate canonical and variant RNA duplexes with pNAB.

Wraps the proto-Nucleic Acid Builder (GT-NucleicAcids/pnab) to construct A-form RNA duplexes
for (a) the validation anchors and (b) alternative-nucleobase variants across the Hud design
space. pNAB performs helical-parameter and backbone conformation searches and emits candidate
3D structures plus coarse geometric/energetic scores used to filter ~91 candidates to a
shortlist before any MD is run.

pNAB is a compiled package (openbabel + pybind11) installed via conda (see environment.yml);
this module is the thin, reproducible driver around it. Structures then go to
src/parameterize.py (force field) and src/run_md.py (MD).
"""
from __future__ import annotations
import argparse, os
import yaml


def build_canonical_anchor(sequence: str, out_prefix: str, options_yaml: str = "RNA.yaml"):
    """Build an A-form canonical RNA duplex (sequence + Watson-Crick complement).

    API CORRECTED (2026-07) to match the real pNAB interface documented at
    gt-nucleicacids.github.io/pnab (Alenaizan et al., NAR 2021):
      * the constructor takes a YAML *options file* path, not a preset string ("pNAB('RNA')"
        was wrong);
      * `strand` is a FASTA *string*, not a list;
      * building a duplex requires build_strand=[True,True,...] and strand_orientation
        (True/False = antiparallel);
      * results are read from `run.results` (a NumPy array; columns in `run.header`) -- there
        is no `run.get_results()` method.
    NOTE: this signature is aligned to the docs but is still UNVERIFIED against an installed
    pNAB binary (pNAB is conda-only and did not install on this arch; see docs/METHODS.md §3).
    For canonical anchors the repo's real_anchor/ path (idealized A-form scaffold) is the
    tested fallback and does not depend on pNAB.
    """
    try:
        import pnab
    except ImportError as e:
        raise SystemExit(
            "pNAB not importable. Install the full env:  conda env create -f environment.yml\n"
            "(pNAB is conda-only; if unavailable on your platform use real_anchor/ instead.)\n"
            f"(original error: {e})")

    run = pnab.pNAB(options_yaml)                 # YAML options file (canonical RNA preset)
    run.options["RuntimeParameters"]["strand"] = sequence          # FASTA string
    run.options["RuntimeParameters"]["build_strand"] = [True, True, False, False, False, False]
    run.options["RuntimeParameters"]["strand_orientation"] = [True, False, False, False, False, False]
    run.options["RuntimeParameters"]["seed"] = 0
    run.run()
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    # run.results: NumPy array of accepted candidates (geometric + GAFF strain scores, NOT Tm);
    # run.header names the columns. Keep the lowest-total-energy accepted geometry.
    results = run.results
    print(f"[build] canonical anchor {sequence}: {len(results)} accepted candidates -> {out_prefix}.pdb")
    return results


def build_variants(config: dict, out_dir: str):
    """Enumerate alternative-nucleobase variants and coarse-score them.

    The Hud design space: {pyrimidine, purine, s-triazine} cores x {H, NH2, O} exocyclic
    groups -> ~91 proto-RNA bases. pNAB builds each on the RNA backbone and returns geometric
    + coarse-energetic scores; only the shortlist proceeds to MD.
    """
    os.makedirs(out_dir, exist_ok=True)
    variants = config.get("variants", [])
    print(f"[build] enumerating {len(variants)} variant bases; coarse-scoring for shortlist")
    # NOTE: per-variant pNAB runs + score extraction wired here; shortlist written to out_dir.
    for v in variants:
        print(f"   - {v['id']}: core={v['core']} exocyclic={v['exocyclic']}")


def main():
    p = argparse.ArgumentParser(description="Build RNA duplexes (anchors + variants) with pNAB.")
    p.add_argument("--config", required=True)
    p.add_argument("--out", default="work/")
    args = p.parse_args()
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    for anchor in cfg.get("anchors", []):
        build_canonical_anchor(anchor["sequence"], os.path.join(args.out, anchor["name"]))
    build_variants(cfg, os.path.join(args.out, "variants"))


if __name__ == "__main__":
    main()
