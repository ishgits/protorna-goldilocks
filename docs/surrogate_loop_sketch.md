# Surrogate active-learning loop — design sketch (stretch goal)

**Motivation.** The `real_anchor/` result shows that **~20 ps of CPU MD already separates duplexes on
the interstrand-energy readout** (GC −69.8 > mixed −57.8 > AU −40.8 kcal/mol, gaps ≈ 3.2 / > 5 pooled SD). That is
the key enabling observation: the discriminating signal is cheap. Full production MD (multi-ns,
explicit-solvent, GPU) on all ~O(10²) geometrically-feasible candidates is the cost bottleneck. An
active-learning surrogate turns the cheap signal into a ranking that lets full MD run only where it
changes the answer.

## The loop

```
                    ┌─────────────────────────────────────────────┐
                    │  ~O(10²) geometrically-feasible candidate     │
                    │  bases   (pNAB build + feasibility filter)    │
                    └───────────────────────┬─────────────────────┘
                                            │  (1) featurize
                                            ▼
   (5) retrain on          ┌───────────────────────────────┐
   new labels              │  SURROGATE  f: features → Ê_int │◄──────────┐
   ┌──────────────────────►│  (GP / gradient-boosted trees)  │           │
   │                       └───────────────┬───────────────┘           │
   │                                       │  (3) predict Ê_int + σ      │
   │                                       ▼   for ALL candidates        │
   │              ┌─────────────────────────────────────────┐          │
   │              │  acquisition = rank by Ê_int, break ties  │          │
   │              │  by σ (uncertainty)  → pick next batch     │          │
   │              └───────────────────────┬───────────────────┘          │
   │                                      │  (4) run the EXPENSIVE label  │
   │                                      ▼   on the top decile only      │
   │        ┌──────────────────────────────────────────────┐            │
   └────────┤  short-MD interaction energy (the ground truth)├────────────┘
            │  ~20 ps OL3/GBn2, the real_anchor recipe       │
            └────────────────────────────────────────────────┘
   Seed (t=0): the 3 anchors + any candidates with published Tm  →  initial labels
   Stop when:  top-decile membership is stable across two iterations
               OR σ over the unlabeled top quartile falls below the anchor gap (~12 kcal/mol)
```

## The five stages

1. **Featurize.** Cheap, physically-motivated descriptors computed without MD: number and geometry of
   Watson–Crick H-bond donors/acceptors on the pairing edge; ring count / aromatic surface area;
   dipole and partial-charge summary from the RESP fit (the parameterization stage already produces these per base);
   stacking-relevant polarizability proxy; pNAB geometric-closure score. ~10–20 features per base.
2. **Train a cheap surrogate.** Regress features → **short-MD interstrand interaction energy** (the
   `real_anchor` observable, *not* absolute Tm). A **Gaussian process** is the natural first choice
   because it returns a calibrated per-candidate uncertainty σ that the acquisition step needs; a
   gradient-boosted-tree ensemble with quantile outputs is the fallback for larger feature sets.
3. **Surrogate-rank all candidates.** Predict Ê_int and σ for the full ~O(10²) set — milliseconds,
   versus GPU-hours per full-MD label.
4. **Label the top decile only.** Run the short-MD interaction-energy calculation (the exact
   `real_anchor/analyze_real.py` readout) on the highest-ranked ~10 candidates — plus a few
   high-σ candidates for exploration (uncertainty sampling), so the loop does not just exploit its
   current belief.
5. **Retrain and repeat.** Fold the new labels back into the training set, refit, re-rank. The full
   *production* multi-ns GPU MD runs only on the candidates that survive to the final stable top
   decile — the top decile of a validated pre-rank.

**Training signal:** short (~20 ps) MD interstrand interaction energy — cheap, already shown to
discriminate. **Stopping criterion:** top-decile membership stable across two consecutive iterations,
or predictive σ over the unlabeled top quartile drops below the ~12 kcal/mol anchor-to-anchor gap
(below which the surrogate can no longer reorder the shortlist).

## Why this stays auditable

The loop is automated end-to-end — featurization from the RESP fit and geometry, dispatching
the short-MD label jobs, refitting the surrogate, and selecting the next batch are the same
orchestration the pipeline already uses, closed into a loop. And it stays auditable: any
ranking that promotes a candidate into the MD queue must be traceable to features and labels,
not to an overfit model. The surrogate never *replaces* MD; it decides *where MD runs*.

## Honest risk

The surrogate is trained on ~20 ps interaction energies for a shortlist of mostly canonical-like
chemistries. It may **not generalize to novel heterocycles** whose stacking or H-bonding is
qualitatively outside the training distribution — precisely the candidates the project most wants to
find. Mitigations: keep an uncertainty-sampling fraction in every batch (label some high-σ candidates
regardless of predicted rank), and treat the surrogate as a *triage* that reorders the MD queue, never
as a final stability verdict. If the top decile churns instead of converging, that is the signal the
features are insufficient — a finding, not a failure.

## Summary

Because ~20 ps of CPU MD already separates the anchors on the interaction-energy readout, a cheap
uncertainty-aware surrogate trained on short-MD energies can pre-rank the full ~O(10²) candidate set,
so production multi-ns GPU MD is spent only on the top decile — a reviewer-auditable
active-learning loop that scales the screen without diluting its rigor.
