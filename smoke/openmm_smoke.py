"""
openmm_smoke.py — REAL OpenMM dynamics, end-to-end into the analysis code.

Purpose: prove (a) the MD engine is correctly installed and integrates dynamics, and (b) the
same analyze.py ranking used on real duplexes consumes MD output and reproduces a known
ordering. This is a deliberately reduced model, NOT the RNA screen: two "base-pair" particles
in a Langevin thermostat, one held by a STRONG harmonic H-bond mimic, one by a WEAK one. The
strong pair should stay tighter (higher persistence, more stable) — and the pipeline must rank
strong > weak. The production pipeline swaps this toy System for pNAB-built RNA duplexes under
AMBER chiOL3 (see src/run_md.py); the analysis layer is identical.

Runs on CPU in seconds. GPU not required for the smoke test.
"""
from __future__ import annotations
import sys, os
import numpy as np

import openmm as mm
from openmm import unit

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from analyze import DuplexResult, rank_duplexes, ordering_matches  # noqa: E402

KB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA


def simulate_pair(k_hbond, r0_nm=0.29, n_steps=20000, report_every=20, temperature=300.0, seed=1):
    """Two particles joined by a harmonic 'H-bond' of stiffness k_hbond (kJ/mol/nm^2).
    Returns per-frame pair distance in ANGSTROM (to match analyze.py's cutoff units)."""
    system = mm.System()
    system.addParticle(120.0)   # amu, lumped "base" mass
    system.addParticle(120.0)

    bond = mm.HarmonicBondForce()
    bond.addBond(0, 1, r0_nm * unit.nanometer, k_hbond * unit.kilojoule_per_mole / unit.nanometer**2)
    system.addForce(bond)

    integrator = mm.LangevinMiddleIntegrator(temperature * unit.kelvin,
                                             5.0 / unit.picosecond,
                                             2.0 * unit.femtosecond)
    integrator.setRandomNumberSeed(seed)

    platform = mm.Platform.getPlatformByName("Reference")  # deterministic, no GPU needed
    context = mm.Context(system, integrator, platform)
    context.setPositions([[0, 0, 0], [r0_nm, 0, 0]] * unit.nanometer)
    context.setVelocitiesToTemperature(temperature * unit.kelvin, seed)

    dists = []
    for _ in range(n_steps // report_every):
        integrator.step(report_every)
        pos = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        d_nm = np.linalg.norm(pos[0] - pos[1])
        dists.append(d_nm * 10.0)  # nm -> Angstrom
    return np.array(dists)


def main() -> int:
    # strong vs weak H-bond mimic -> known ordering: strong is more stable
    strong = simulate_pair(k_hbond=150000.0, seed=1)
    weak = simulate_pair(k_hbond=8000.0, seed=2)

    # feed into the SAME analysis used on real duplexes. RMSD proxy = deviation from r0 (2.9 A).
    def as_result(name, d):
        bp = d.reshape(-1, 1)                 # (frames, 1 pair)
        rmsd = np.abs(d - 2.9)                # distortion from equilibrium pairing
        return DuplexResult.from_traj(name, bp, rmsd)

    results = [as_result("strong_pair", strong), as_result("weak_pair", weak)]
    ranked = rank_duplexes(results)
    expected = ["strong_pair", "weak_pair"]

    print(f"OpenMM {mm.version.version} — Reference platform")
    for r in ranked:
        print(f"  {r.name:12s}  persistence={r.mean_bp_persistence:.3f}  "
              f"rmsd_stab={r.rmsd_stability:.3f}  score={r.stability_score:.4f}")
    print(f"  <d> strong={strong.mean():.2f} A   <d> weak={weak.mean():.2f} A")

    ok = ordering_matches(ranked, expected)
    print(("PASS" if ok else "FAIL") +
          f": real OpenMM run -> analyze.py ranked {' > '.join(r.name for r in ranked)} "
          f"(expected {' > '.join(expected)}).")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
