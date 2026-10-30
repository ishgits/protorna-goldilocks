# proto-RNA Goldilocks — duplex-stability screening pipeline

A reproducible **build → filter → MD → rank** pipeline for measuring the *relative duplex
stability* of nucleic acids carrying alternative (non-canonical) nucleobases, and for
locating which bases fall inside the "Goldilocks" stability window that functional genetic
polymers — prebiotic proto-RNA and modern therapeutic oligonucleotides alike — require.

This repository is a **validation stub** for the screening pipeline. Its
job is narrow and deliberate: **build toward reproducing known relative stabilities before any
novel base is screened.** Under-promise coverage, over-deliver rigor.

> **Status (honest):** `src/validate_analysis.py` unit-tests the *ranking logic* on synthetic
> input (ordering true by construction — a logic test, **not** reproduced physics).
> `real_anchor/` is the real result: real OL3 MD on real A-form coordinates recovers the
> **GC > mixed > AU** ordering via an interstrand-energy readout, rank-consistent with Turner
> nearest-neighbor thermodynamics. The *geometric* `stability_score` needs longer (multi-ns) MD
> to reproduce the ordering. See `docs/METHODS.md`.

---

## Why this exists

Duplex stability is a *precondition* for function. Too unstable → a strand loses register
and cannot template or hybridize. Too stable → it will not strand-separate. Nonenzymatic
template copying (Szostak) and antisense/siRNA activity both live inside a narrow stability
window. The scientific question: **which alternative nucleobases land in that window?**

pNAB can *build* many candidate proto-RNAs across the Hud-group design space
(pyrimidine / purine / s-triazine cores × {H, NH₂, O} exocyclic groups → ~O(10²) candidates;
each non-canonical base needs a reference-frame-aligned geometry). Its coarse score is a
**geometric-feasibility** filter (GAFF strain + helix-closure) — *is this base buildable* —
**not** a stability or Tm estimate. Defensible explicit-solvent MD can *rank* only tens, and MD
is where stability is actually measured. This pipeline bridges the two: feasibility-filter
everything, then run trustworthy MD on a curated shortlist.

**The observable is relative rank-ordering (RMSD, base-pair persistence, MD interaction-energy →
*relative* Tm ranking), not absolute ΔG.** That is the epistemically honest readout, and the thing
this stub validates.

---

## Pipeline

```
  pNAB build            pNAB coarse score        explicit-solvent MD         analysis
 ┌────────────┐  →  ┌──────────────────────┐ → ┌───────────────────────┐ → ┌───────────────┐
 │ canonical +│     │ geometric + energetic│   │ AMBER χOL3/OL3 + water │   │ RMSD, bp      │
 │ variant    │     │ filter to shortlist  │   │ (OpenMM, GPU)          │   │ persistence,  │
 │ RNA duplex │     │                      │   │ on Purdue RCAC / GPU   │   │ MD→Tm rank    │
 └────────────┘     └──────────────────────┘   └───────────────────────┘   └───────────────┘
        │                                                                          │
        └────────────── validation anchors gate the whole thing ──────────────────┘
```

### Validation-first (non-negotiable)
Novel-base screening only proceeds once the pipeline reproduces the **correct ordering** of:
1. **Canonical anchor** — short A-RNA duplexes / Turner tetramers vs. Turner-lab thermodynamic
   data (Šponer / Banáš OL3 protocols).
2. **Non-canonical anchor** — Benner Z–P / hachimoji-RNA ladder (experimental Tm predicted to
   ~2 °C), to prove the pipeline generalizes past Watson–Crick pairs before screening novel bases.

Held line: **RNA-native throughout** (no DNA anchors — consistency beats benchmark richness).
Do **not** start on Hud proto-RNA rosettes (barbituric acid / melamine / TAP) — supramolecular
assemblies are noisier MD; validate on Watson–Crick-type pairs first.

---

## Repo layout

```
protorna-goldilocks/
├── README.md
├── requirements.txt          # OpenMM + analysis deps (pip)
├── environment.yml           # full conda env incl. pNAB, AmberTools
├── configs/
│   └── duplexes.yaml         # benchmark sequences + expected stability ordering
├── data/
│   └── turner_benchmarks.csv # reference thermodynamics for the anchors
├── src/
│   ├── build_duplex.py       # pNAB build of canonical + variant RNA duplexes
│   ├── parameterize.py       # χOL3 assignment; non-canonical FF-by-analogy scaffold
│   ├── run_md.py             # OpenMM explicit-solvent MD driver (GPU-aware)
│   ├── analyze.py            # RMSD, base-pair persistence, stability score → ranking
│   ├── validate_analysis.py  # SELF-CONTAINED: unit-tests the ranking LOGIC on synthetic input
│   └── nn_predict.py         # Turner-1998 NN dG37 for the anchors (pins the ordering gate)
├── real_anchor/             # >>> REAL physics result (no synthetic input, CPU-only) <<<
│   ├── build_anchor.py       # exact anchors on a real 1RNA A-form scaffold (base mutation)
│   └── analyze_real.py       # real OL3 MD → analyze.py + MM-GBSA interstrand energy
├── smoke/
│   └── openmm_smoke.py       # real OpenMM micro-MD, proves the engine is wired
├── gromacs/                  # >>> PRIMARY CLUSTER PATH (Purdue RCAC) <<<
│   ├── mdp/{em,nvt_eq,npt_eq,npt_prod}.mdp   # RNA-duplex MD, 300 K, OL3 + TIP3P
│   ├── build_system.sh       # nab -> tleap (OL3) -> acpype -> GROMACS (run once, AmberTools)
│   └── analyze_gmx.py        # parse gmx .xvg -> reuse analyze.py ranking (tested)
├── cluster_jobs/            # cluster-specific SLURM scripts (the real submit path)
│   ├── negishi/              # account eapsdept, partition gpu, ROCm container GROMACS
│   └── gautschi/             # account pearce21, partition smallgpu, module GROMACS 2024.3
├── docs/                    # methods, repro log, param plan, figure caption
├── results/                 # anchor_ordering_v2.{png,svg}, real_anchor_ordering.png, plot_anchor_ordering_v2.py
├── slurm/run_md.sbatch       # OpenMM/CUDA template (portable / local NVIDIA path)
└── .github/workflows/ci.yml  # runs validate_analysis + openmm_smoke on every push
```

## Two execution paths (read this)

- **OpenMM path** (`src/run_md.py`, `smoke/`) — portable, CUDA/NVIDIA, and the CI proof-of-work.
  It runs anywhere and is what the automated tests exercise. Use it on a local **NVIDIA** GPU.
- **GROMACS path** (`gromacs/`) — the **cluster** path, targeting **Purdue RCAC Negishi**
  (ROCm/AMD GPUs, container GROMACS, account `eapsdept`; a Gautschi variant under account
  `pearce21` is also provided). SLURM scripts and `.mdp` files are
  cloned from a proven Pearce-group workflow. **This is the route to run real duplex MD on
  Negishi** — OpenMM/CUDA is not the path of least resistance on AMD/ROCm hardware.

Both paths converge on the *same* tested ranking logic in `src/analyze.py`.

## Quickstart: tested local checks

These run anywhere, need no GPU, and are what CI exercises:

```bash
# environment (analysis/MD engine only; pip)
pip install -r requirements.txt

python src/validate_analysis.py           # ranking LOGIC unit test (synthetic input, <10 s)
python src/nn_predict.py                   # Turner NN dG37 — reproduces GC > mixed > AU ordering
python results/plot_anchor_ordering_v2.py  # regenerates results/anchor_ordering_v2.{png,svg}
```

If OpenMM, PDBFixer, and MDTraj are installed, reproduce the real anchor result:

```bash
python smoke/openmm_smoke.py               # real OpenMM micro-MD, proves the engine is wired
cd real_anchor && python analyze_real.py   # real OL3 MD on 1RNA A-form coords → GC > mixed > AU
```

## Production scaffold

The full pNAB → χOL3 parameterization → production-MD workflow is documented but not
directly runnable from a clean checkout: it requires AmberTools, pNAB, GROMACS, and
cluster-specific module setup, and `configs/duplexes.yaml` still carries a placeholder
non-canonical target. `src/build_duplex.py` does not yet emit the `prmtop`/`inpcrd` that
`src/run_md.py` consumes — wiring that end-to-end on non-canonical bases is the funded next
step. For the intended production path see:

- `gromacs/build_system.sh` — the one AmberTools step (build → solvate → GROMACS topology)
- `cluster_jobs/README.md` — Negishi / Gautschi SLURM variants
- `docs/param_plan_inosine.md` — the QM/RESP parameterization plan for the first non-canonical base

## Running on the cluster (Negishi / GROMACS) — exact order

```bash
# STEP 0 (once, AmberTools env — login/interactive node, NOT a GPU job):
#   builds solvated+ionised OL3 duplexes into work/. The ONLY step needing AmberTools+acpype.
bash gromacs/build_system.sh
#   -> work/{gc_anchor,mixed_anchor,au_anchor}.gro  +  work/topol_*.top

# STEP 1-3 (Negishi GPU jobs — pure GROMACS, your proven stack):
# (Gautschi variants live under cluster_jobs/gautschi/)
sbatch cluster_jobs/negishi/slurm_em_nvt.sh      # EM + NVT (array over the 3 anchors)
sbatch cluster_jobs/negishi/slurm_npt_eq.sh      # NPT equilibration (after EM/NVT finish)
sbatch cluster_jobs/negishi/slurm_npt_prod.sh    # 50 ns production + writes rmsd_*.xvg, hbnum_*.xvg

# STEP 4 (local, tested): rank the anchors and check the ordering gate
python gromacs/analyze_gmx.py --work work
#   -> PASS/FAIL: did GROMACS reproduce gc_anchor > mixed_anchor > au_anchor?
```

**Submit-ready status:** the `.mdp` files, all three `slurm_*.sh` scripts, and `analyze_gmx.py`
are ready for Negishi as written (edit only the `#SBATCH -A`/paths if yours differ, and create
the `RNA` / `strand1` / `strand2` index groups once with `gmx make_ndx` — noted inline in the
scripts). **Step 0 is the one part you must run in an AmberTools environment**; verify the
`nab wc_helix` call against your AmberTools version, or drop in a canonical A-RNA duplex PDB.

## What is validated *in this stub* vs. what runs on your GPU

| Stage | Status here | Where it really runs |
|---|---|---|
| Ranking / stability-score logic (`validate_analysis.py`) | ✅ runs anywhere, tested | same code, real trajectories |
| OpenMM MD engine (`smoke/openmm_smoke.py`) | ✅ real micro-run, tested | scaled to duplexes on GPU |
| pNAB build, χOL3 param, production MD | scaffolded, documented | Purdue RCAC / local GPU |

The point of the stub: **the reasoning and readout that the whole screen depends on already
run and are tested.** Scaling to the RNA systems is a compute question, not a "does it work"
question.

## Where the reasoning load sits
The hard part is not orchestration glue — it is the *reasoning loop*: (1) non-canonical force-field
parameterization, where the published norm is **QM/RESP charge derivation** case-by-case — prioritizing
which bases need fresh QM, proposing and justifying analogy seeds (GAFF/RESP precedents + OL3
protocols), and auditing the RESP fits; this is the project's critical path; (2) trajectory triage /
convergence-flagging across dozens of runs; (3) reviewer-style auditing of every MD interaction-energy
→ relative-Tm-ranking claim and citation. See `src/parameterize.py` for the parameterization hooks and
`docs/METHODS.md` §2 for the QM/RESP-primary rationale.

## Key references
- **pNAB** — Alenaizan et al., *Nucleic Acids Research* 49:79 (2021), GT-NucleicAcids/pnab.
- **Hachimoji / Z–P thermodynamics** — Hoshika et al., *Science* 363:884 (2019).
- **OL3 / χOL3 RNA force field** — Zgarbová et al. (2011); Šponer/Banáš benchmark protocols.
- **Nonenzymatic template copying** — Szostak group (copying kinetics/fidelity).
- **Therapeutic-oligo MD→Tm** — Lomzov / Vorobjev / Pyshnyi protocol.

## License
MIT — see `LICENSE`.
