# Methods & findings — proto-RNA Goldilocks

Technical reference for the validation results in this repository: how the real anchor
result was produced and how to reproduce it, the parameterization approach for
non-canonical bases, and the scope/limitations of each pipeline stage. Every quantitative
claim below traces to a committed command or a cited source.

---

## 1. The real anchor result (`real_anchor/`)

The one physics result in this repo recovers the known **GC > mixed > AU** relative-stability
ordering of three canonical RNA 6-mers, produced with lightweight CPU compute only
(OpenMM + PDBFixer + mdtraj; no AmberTools/pNAB, no GPU, no cluster).

- **Scaffold.** Real A-form backbone = central 6 bp of **PDB 1RNA** (A-RNA duplex). The three
  `configs/duplexes.yaml` anchors are built on one identical backbone by base mutation
  (A→G, U→C on both strands preserves Watson–Crick complementarity) — a *controlled* comparison.
- **MD + readout.** Real OL3/GBn2 implicit-solvent MD (~20 ps, 340 K) fed through the repo's
  tested `src/analyze.py`, plus an **MM-GBSA-style interstrand interaction energy** (the
  "MD interaction-energy → relative-Tm" readout).
- **Result.** Interstrand energy **gc −69.8 ± 4.0 > mixed −57.8 ± 3.3 > au −40.8 ± 2.7 kcal/mol**
  — correct ordering, adjacent gaps large relative to per-frame fluctuations (GC–mixed ≈ 3.2
  pooled SD, mixed–AU > 5 pooled SD), rank-matching Turner NN ΔG°₃₇.
- **Honest limitation.** At ~20 ps the purely *geometric* `stability_score` does **not** resolve
  the ordering (`ordering_reproduced: false`) — all base pairs stay closed; the discriminating
  signal is energetic, not geometric. This is informative, not a failure: it shows the geometric
  composite needs longer (multi-ns) MD, while the interstrand-energy readout is the
  sampling-robust discriminator at short sampling.

### 1B. Reproduce recipe

Clean-room regeneration of the headline artifact. Write to a fresh directory so the committed
`real_anchor/work/` is left untouched (see `docs/repro_log_real_anchor.md` for a logged run):

```bash
pip install openmm pdbfixer mdtraj
cd real_anchor
for a in gc mixed au; do
  python build_anchor.py --scaffold scaffold_1rna.pdb --anchor $a --out work --ns 0.02 --temp 340
done
python analyze_real.py            # -> work/anchor_ordering.json ; figure at results/real_anchor_ordering.png
```

Expected: energy ordering `gc > mixed > au` reproduces; energies match committed values within
the stated per-frame SD. A small energy drift between runs is expected and fine; a *changed
ordering* would be the meaningful failure.

### 1C. The ordering gate (Turner NN)

`data/turner_benchmarks.csv` carries NN-predicted ΔG°₃₇ for all anchors, computed reproducibly by
`src/nn_predict.py` (Turner 1998 / NNDB): **GCGCGC −10.46 > GCAUGC −7.64 > AUAUAU −0.54 kcal/mol**,
with per-row source and verify flags. This makes the ordering gate falsifiable.

- Reproduce: `python src/nn_predict.py`
- **Caveat.** AUAUAU at ΔG°₃₇ ≈ −0.54 kcal/mol is a *marginal* duplex (very low Tm). It is a clean
  "least-stable extreme" validation anchor, not a calibrated low-stability gate; a slightly stiffer
  AU-containing 6-mer would melt in a more measurable range.

---

## 2. Force-field parameterization for non-canonical bases (QM/RESP-primary)

For non-canonical / XNA bases the published norm is **QM/RESP charge derivation, case-by-case**;
analogy is used to *seed* and *audit*, never to finalize the charges. Charge-by-analogy alone is
not a defensible parameterization method for a novel base.

- Self-consistent non-canonical DNA params derive charges by QM/RESP —
  [Life 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9143812/).
- The Ducque XNA builder generates force-field params from QM —
  [NAR 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11014352/).

The practical protocol: QM/RESP is primary; analogy reasoning is used to (i) prioritize which bases
need fresh QM, (ii) propose and justify bonded-term seeds from canonical neighbours, and (iii) audit
the RESP fit atom-by-atom against a published parameter set. See `docs/param_plan_inosine.md` for a
concrete worked plan (inosine), and `src/parameterize.py` for the parameterization hooks.

---

## 3. Pipeline-stage scope (what each stage measures)

**pNAB build + coarse score = geometric feasibility, not stability.** pNAB can build antiparallel
RNA duplexes and arbitrary bases, but each non-canonical base needs a **user-supplied 3D geometry
aligned to the standard reference frame** (only A/G/C/U/T + two hexad bases are built in). Its coarse
score is a **GAFF bond/angle/torsion/vdW strain + geometric helix-closure** filter — *is this base
buildable into a periodic helix* — **not** a Tm or stability estimate. The funnel is therefore:
pNAB filters on **geometric feasibility**; **MD does the stability ranking**.
[pNAB, NAR 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC7797056/) ·
[pNAB examples](https://gt-nucleicacids.github.io/pnab/html/examples.html).

> **Build note.** pNAB and AmberTools are conda-only and do not install on arm64; for the canonical
> anchors, `real_anchor/` is the tested, pNAB-free path. `src/build_duplex.py` was corrected to the
> documented pNAB interface (YAML options file; `strand` FASTA string; `build_strand` /
> `strand_orientation` for duplexes; `run.results` NumPy array) but is still unverified against an
> installed binary.

**OL3 short-construct caveats.** OL3 is the standard RNA force field but has documented failure modes
on very short constructs — spurious intercalated states and over-stabilization in some tetramers
(GACC, CAAU, CCCC), addressed by later reparameterization. Prefer the **6-mer Watson–Crick anchors**
as the primary gate; treat tetramers as secondary and state the intercalation caveat.
[Zgarbová 2011 revision](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5312698/) ·
[OL3 base-phosphate fix, JCTC 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10687871/).

**MD interaction-energy → Tm is *relative* only.** The MD-binding-energy → Tm correlation is real but
largely DNA / RNA-DNA-hybrid derived and calibrated on families of similar oligos — good for
*relative* ranking, not absolute Tm to a few °C. The honest readout throughout is relative
rank-ordering, never absolute ΔG or absolute Tm.
[Lomzov JCIM 2011](https://pubs.acs.org/doi/abs/10.1021/ci200141j) ·
[RNA/DNA MD ΔG, 2021](https://link.springer.com/article/10.1134/S002689332105006X).

---

## 4. Prior-art positioning

The closest prior art is **backbone/sugar-modified** analog MD (LNA/PNA/morpholino/phosphorothioate)
and **XNA** pipelines (Ducque), plus modified-ribonucleotide MD — none is a nucleobase-space stability
screen over the Hud proto-RNA design space. The precise framing: a **nucleobase-space stability screen
for proto-RNA**, distinct from the well-trodden backbone/sugar-analog MD, with Ducque as the
methodological neighbour (and as evidence QM params are standard).
[Ducque, NAR 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11014352/) ·
[comparative MD of modified analogs, NARGAB 2024](https://academic.oup.com/nargab/article/6/4/lqae155/7912069).

---

*Standard held: under-promise coverage, over-deliver rigor. Every number in `real_anchor/` and
`data/` is reproducible from a committed command; each claim is labelled known / inferred / uncertain
in the docs that use it.*
