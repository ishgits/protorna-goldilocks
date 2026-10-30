"""
parameterize.py — force-field assignment, incl. the non-canonical critical path.

Canonical RNA residues use AMBER chiOL3/OL3 directly. The hard, timeline-critical work is
assigning DEFENSIBLE bonded + partial-charge (RESP) parameters to novel heterocycles for which
no validated force field exists. This module is the seam where careful expert review is load-bearing:
it drafts by-analogy parameter sets, flags where analogy breaks and QM charge derivation is
required, and hands off to a RESP workflow.

Two paths per non-canonical base:
  (1) BY-ANALOGY  — map atoms to the closest canonical/GAFF environments; inherit bonded terms;
                    estimate charges from analogous groups. Fast, needs expert justification.
  (2) QM/RESP     — antechamber/psi4 -> RESP fit for charges; frcmod for missing bonded terms.
                    Slower, but embarrassingly parallel across bases and terms.

Review hook (see docs): a proposed parameter set must be audited against GAFF/RESP precedent
+ Sponer/Banas OL3 protocols before use; the audit either approves the analogy set or flags
it to escalate to QM/RESP.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Literal


@dataclass
class BaseParamPlan:
    base_id: str
    core: Literal["pyrimidine", "purine", "s-triazine"]
    exocyclic: list[str]                       # e.g. ["NH2", "O", "H"]
    strategy: Literal["by_analogy", "qm_resp"] = "by_analogy"
    analogy_map: dict[str, str] = field(default_factory=dict)   # novel_atom -> canonical/GAFF type
    needs_qm: bool = False
    justification: str = ""


def propose_by_analogy(base_id: str, core: str, exocyclic: list[str]) -> BaseParamPlan:
    """First-pass analogy mapping. The mapping and justification produced here are a scaffold
    for expert review: each assignment is checked against the GAFF/OL3 atom-type tables +
    Sponer/Banas protocol notes, and escalated to QM/RESP (needs_qm=True) where the electronic
    environment has no good analogue."""
    plan = BaseParamPlan(base_id=base_id, core=core, exocyclic=exocyclic, strategy="by_analogy")
    # first-pass heuristic; a reviewed analogy map supersedes this in production:
    if "O" in exocyclic and core == "s-triazine":
        plan.needs_qm = True
        plan.strategy = "qm_resp"
        plan.justification = ("s-triazine + carbonyl O has no close canonical analogue; "
                              "escalate to QM/RESP for charges and ring bonded terms.")
    else:
        plan.justification = "exocyclic groups map onto canonical purine/pyrimidine environments."
    return plan


def to_openmm_residue_xml(plan: BaseParamPlan, out_path: str):
    """Emit an OpenMM residue-template XML for this base (consumed by run_md.py).
    For qm_resp plans, charges come from the RESP fit; for by_analogy, from inherited types."""
    # writer wired to antechamber/parmed output; stub records intent for the scaffold.
    with open(out_path, "w") as fh:
        fh.write(f"<!-- residue template for {plan.base_id} "
                 f"(strategy={plan.strategy}, needs_qm={plan.needs_qm}) -->\n")
    return out_path


if __name__ == "__main__":
    demo = propose_by_analogy("triazine_OOO", "s-triazine", ["O", "O", "O"])
    print(demo)
