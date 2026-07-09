# QM/RESP parameterization plan — inosine (hypoxanthine riboside) in the OL3/χOL3 framework

**Scope.** A concrete, reproducible protocol for deriving a defensible AMBER-compatible parameter set
for **inosine** as an RNA residue, compatible with the `leaprc.RNA.OL3` (χOL3) force field the rest of
the pipeline uses. It is the concrete "QM/RESP on one non-canonical base" worked example for the
QM/RESP-primary approach in `docs/METHODS.md` §2: **QM/RESP is primary; analogy is a
seed and an auditor, never a replacement for the charge derivation.**

**Why inosine.** Cleaner precedent than 2-aminopurine, and common in real RNA (wobble position,
A-to-I editing). Its base — **hypoxanthine** — is a 9-substituted purine that is structurally
guanine minus the C2 exocyclic amino group (C2 bears H instead of –NH₂), with the same C6=O6
carbonyl and N1–H imidazole/pyrimidine skeleton. That gives two independent canonical neighbours to
seed bonded terms from (guanine for the six-membered/carbonyl side, adenine for the aromatic C2–H),
**and** a published RESP parameter set to validate the final charges against (Aduri 2007).

> **Honest framing throughout.** Every atom-typing and bonded decision below is labelled
> **[SEED-by-analogy]** (inherited, then verified) or **[QM-required]** (must be derived/scanned, no
> defensible analogue). Partial charges are **[QM-required] in full** — analogy may sanity-check them
> but must not set them.

---

## 0. What is *not* re-derived (inherited verbatim from OL3)

The ribose and phosphodiester backbone of inosine are chemically identical to canonical ribonucleotides.

| Component | Decision | Reason |
|---|---|---|
| Ribose ring bonds/angles/torsions (C1′–C2′–C3′–C4′–O4′, 2′-OH) | **[SEED-by-analogy]** — inherit OL3 verbatim | Sugar is unmodified; re-deriving would only add noise and break backbone consistency. |
| Phosphodiester backbone (P, OP1/OP2, O5′, O3′, α/β/γ/δ/ε/ζ torsions) | **[SEED-by-analogy]** — inherit OL3 verbatim | Identical to canonical RNA; χOL3 backbone refit applies unchanged. |
| Backbone + sugar **partial charges** | **[SEED-by-analogy]** — constrain to canonical ff99/OL3 values during the RESP fit (see §2) | Keeps the base swappable into an OL3 helix with no charge discontinuity at the glycosidic bond. |

Only the **base (hypoxanthine) charges** and any **base-specific bonded terms** are new work.

---

## 1. Level of theory for RESP charge derivation

**Use the canonical AMBER/ff99→OL3 recipe, unchanged, so the new base is charge-consistent with the
force field it is inserted into.** Mixing a "better" level of theory would break additivity.

1. **Geometry optimization:** HF/6-31G\* (gas phase). *(Rationale below; a B3LYP/6-31G\* optimized
   geometry is an acceptable variant provided the ESP single point stays HF/6-31G\*.)*
2. **ESP single point:** **HF/6-31G\*** with a Merz–Singh–Kollman (MK) grid.
   - This is the deliberate choice of the AMBER charge model (Bayly 1993; Cornell 1995): HF/6-31G\*
     systematically **over-polarizes** gas-phase dipoles by ~10–20%, which empirically mimics the
     mean condensed-phase polarization a *fixed-charge, non-polarizable* force field cannot represent
     explicitly. Using a larger basis or a correlated method would *under*-polarize relative to the
     rest of OL3 and is the classic mistake when extending AMBER to a new residue.
3. **Charge fit:** **two-stage RESP** (Bayly 1993) — stage 1 hyperbolic restraint weight *a* = 0.0005,
   stage 2 *a* = 0.001, restraining non-hydrogen charges toward zero; stage 2 refits only the
   equivalenced groups.

**Precedent that this is the norm for non-canonical bases (not analogy):**
- **Aduri et al., *J. Chem. Theory Comput.* 3:1464–1475 (2007), DOI 10.1021/ct600329w** — AMBER
  parameters for all **107** naturally occurring modified nucleosides in RNA, **including inosine**
  (given the residue-code letter "I"). Charges derived exactly this way: **HF/6-31G(d) RESP,
  deliberately chosen for consistency with the original Cornell/Cieplak charge model rather than a
  higher level of theory**, with the sugar moiety charge-constrained to the canonical value. Their
  own canonical-nucleoside charges reproduced Cieplak et al. to within a **0.036 e** standard
  deviation — which is the yardstick for the §3 acceptance bounds. Bonded (bond/angle/dihedral) terms
  in Aduri were taken from **GAFF**, so treat their *charges* as the validation target and re-derive
  bonded terms under the modern OL3 refit (see next bullet).
- **modXNA — Wang/Cheatham-style modular parametrization, *J. Chem. Theory Comput.* (2024),
  DOI 10.1021/acs.jctc.4c01164** — explicitly notes Aduri's 2007 set predates OL3 (2011) and lacks the
  updated torsions/radii; provides a modern, OL3-consistent route. Use modXNA's inosine charges as a
  **second** validation cross-check where available; if the fitted set agrees with *both* Aduri and
  modXNA it is well-anchored.
- **Ducque XNA builder, *Nucleic Acids Res.* (2024)** and the self-consistent non-canonical-DNA
  parameter work (*Life* 12:746, 2022) — both derive charges by QM/RESP case-by-case; analogy is used
  only to seed. See `docs/METHODS.md` §2 for the QM/RESP-primary rationale and citations.

---

## 2. Bonded terms — inherit by analogy vs. re-derive, with the physical reasoning

Atom naming follows the standard purine convention (N1, C2, N3, C4, C5, C6, O6, N7, C8, N9; H at
N1, C2, C8; glycosidic N9–C1′).

| Term | Decision | Physical reasoning |
|---|---|---|
| Imidazole ring (N7, C8, N9, C4, C5) bonds/angles | **[SEED-by-analogy]** from guanine/adenine (identical 5-membered ring) | The imidazole ring is chemically identical across all purines; OL3 guanine terms transfer exactly. |
| C6=O6 carbonyl + N1–H (bond, angle, in-plane) | **[SEED-by-analogy]** from **guanine** | Hypoxanthine's 6-oxo/N1–H edge is the guanine 6-oxo/N1–H edge; same hybridization, same lactam tautomer at physiological pH. |
| **C2–H2 (aromatic CH)** and the N1–C2, C2–N3 ring bonds | **[SEED-by-analogy]** from **adenine C2–H2** | Guanine has –NH₂ at C2, so it gives no C2–H term; adenine's C2–H2 is the correct analogue for the aromatic CH. This is the one place a single canonical neighbour is insufficient — hence the guanine⊕adenine chimera. |
| Ring bond *lengths/force constants* around C2 (the deletion site) | **[SEED-by-analogy], then VERIFY** with `parmchk2`; escalate any `ATTN`/penalty-flagged term to a short QM scan | Removing the amino group changes the electronic density in the pyrimidine ring; equilibrium bond lengths shift slightly. `parmchk2` will flag terms it had to guess — those, and only those, get QM. |
| **Glycosidic torsion χ (O4′–C1′–N9–C4)** | **[QM-required] — run a relaxed χ scan** and compare to the OL3 canonical-purine χ profile before adopting it | χOL3 refit χ *for canonical purines*. Inosine's base electronics differ enough (no C2 amino donor) that the χ profile should be checked, not assumed. If the scan overlays the OL3 purine profile within ~0.5 kcal/mol, adopt OL3 χ **[then it becomes SEED]**; if it deviates, refit. This is the single most important bonded check for a purine analogue and the most common silent error. |
| Any bonded term `parmchk2` cannot assign from GAFF2/parm | **[QM-required]** frcmod entry from the QM Hessian | No defensible analogue exists; derive it rather than let `parmchk2` interpolate. |

**Charges (all base atoms): [QM-required] in full.** You cannot delete guanine's N2 charge and
re-normalize — removing the amino group redistributes ESP density across the *entire* ring, not just
at C2. The whole base must be re-fit by RESP. Analogy's only role here is auditing (§3).

**Two hard constraints on the fit (both auditable):**
- **Net integer charge.** The capped fragment fits to an exact integer (e.g. −1 for the
  5′-monophosphate nucleotide, or a neutral capped nucleoside with backbone charges constrained to the
  canonical set). Verify: `sum(charges)` rounds to the integer within <1e-4 e.
- **Symmetry equivalence.** Symmetry-equivalent atoms get equal charges: methyl/methylene hydrogens on
  any phosphate cap, and the two 5′ H5′/H5′′ (equivalenced in stage-2 RESP). The hypoxanthine base
  itself has no internal symmetry, so no base atoms are equivalenced — flag it if the fit accidentally
  ties two ring atoms.

---

## 3. Validation check — compare fitted charges against published inosine parameters

**Primary target: Aduri et al. 2007** inosine residue (distributed with AmberTools as the modified-RNA
library, historically `all_modrna08.lib` / residue code often `I`/`RI`). Same level of theory, same
capping philosophy → agreement is expected and *disagreement is diagnostic*.

**Procedure (do this in `parameterize.py`'s audit hook, not by eye):**
1. Extract the per-atom base charges from your RESP fit and from the Aduri set.
2. Report, per atom: `q_fit`, `q_Aduri`, `Δ = q_fit − q_Aduri`.
3. Summary statistics: max |Δ|, RMSD over base atoms, and Σq (both should be the same integer).

**Acceptance / what the deltas mean:**
- **RMSD ≲ 0.036 e, max |Δ| ≲ 0.08 e:** consistent with the published set — pass. The 0.036 e figure
  is the exact Aduri-vs-Cieplak canonical-nucleoside agreement, so it is the natural threshold.
  Residual differences are the known RESP conformer-dependence (single- vs multi-conformer fit),
  the Aduri-uses-GAFF-bonded-terms difference, and capping details.
- **max |Δ| ≳ 0.1 e localized at C2/N3/O6:** re-examine the optimized geometry (tautomer? a spurious
  intramolecular H-bond in the gas-phase conformer inflates local charges) before trusting the fit.
- **A whole-ring systematic offset:** almost always a backbone-charge-constraint mismatch (you fit the
  base without pinning the sugar to canonical values) — refit with the constraint of §2.

Report the comparison table in the parameterization note; a fit that cannot be reconciled with Aduri
within these bounds does **not** proceed to production MD.

---

## 4. Exact command sequence

Two interchangeable QM front-ends; the RESP fit and downstream AMBER steps are identical.

### 4a. QM ESP — Gaussian route
```bash
# 1. geometry optimization (HF/6-31G*), capped inosine nucleoside/nucleotide fragment
#    inosine_opt.gjf header:
#    %chk=inosine.chk
#    #P HF/6-31G* Opt SCF=Tight
g16 < inosine_opt.gjf > inosine_opt.log

# 2. ESP single point on the optimized geometry, Merz-Kollman grid written for RESP
#    inosine_esp.gjf header (read geometry from chk):
#    #P HF/6-31G* SCF=Tight Pop=MK IOp(6/33=2,6/42=6,6/50=1)
g16 < inosine_esp.gjf > inosine_esp.log
```

### 4b. QM ESP — Psi4 route (open-source; matches the Modal/Q7 burst)
```python
# inosine_resp.py  — geometry opt + MK ESP + two-stage RESP via the `resp` package
import psi4, resp
mol = psi4.geometry(open("inosine_capped.xyz").read())
psi4.set_options({"basis": "6-31G*"})
psi4.optimize("scf", molecule=mol)                       # HF/6-31G* geometry
opts = {"VDW_SCALE_FACTORS": [1.4, 1.6, 1.8, 2.0],       # MK-like shells
        "VDW_POINT_DENSITY": 1.0, "RESP_A": 0.0005,       # stage-1 restraint
        "RESP_B": 0.1, "BASIS_ESP": "6-31G*", "METHOD_ESP": "scf"}
resp.resp([mol], opts)                                    # stage 1
opts["RESP_A"] = 0.001                                    # stage-2 restraint
# (equivalence H5'/H5'' and any cap methyls here, then re-run stage 2)
```

### 4c. Fit + assemble the AMBER residue (AmberTools, both routes converge here)
```bash
# two-stage RESP with backbone/sugar charges constrained to canonical (respgen writes the constraints)
antechamber -i inosine_esp.log -fi gout -o inosine.ac -fo ac -c resp -at amber -rn INO -nc -1
# missing bonded terms (flags anything without a defensible analogue -> QM those)
parmchk2 -i inosine.ac -f ac -o inosine.frcmod -a Y      # -a Y prints the analogy source per term
# build a test duplex and confirm it assembles under OL3
tleap -f - <<'EOF'
source leaprc.RNA.OL3
INO = loadmol2 inosine.mol2
loadamberparams inosine.frcmod
# ...load a short duplex with INO substituted, check charge & missing-param report...
check INO
saveamberparm testduplex testduplex.prmtop testduplex.inpcrd
quit
EOF
```

**Watch-outs baked into the commands:**
- `parmchk2 -a Y` prints the analogy source for every bonded term it assigns — that log **is** the
  seed-vs-QM audit trail; any term sourced from a distant atom type is a QM-escalation candidate.
- `-nc -1` (net charge) must match the capping; a non-integer RESP sum means the fragment/cap is wrong.
- Run the χ scan (§2) as a separate relaxed dihedral scan; do not trust the inherited χOL3 profile until
  it overlays.

---

## 5. One-paragraph summary for Q4/Q6

> For inosine, the ribose–phosphate backbone and the imidazole/6-oxo edges are inherited verbatim from
> χOL3 and guanine; the aromatic C2–H is seeded from adenine; and the complete base charge set plus the
> glycosidic χ torsion are derived by HF/6-31G\* RESP and a relaxed χ scan — the published standard
> (Aduri 2007; Ducque 2024), against which the fitted charges are audited atom-by-atom. Claude's role is
> to *prioritize* which terms need QM (χ and the full base charges) versus which inherit safely, to
> *propose and justify* the guanine⊕adenine analogy seeds, and to *audit* the RESP fit against the
> published set — not to replace the QM.
