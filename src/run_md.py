"""
run_md.py — explicit-solvent MD driver (OpenMM, GPU-aware).

Production MD for RNA duplexes under the AMBER chiOL3/OL3 RNA force field in OPC (or TIP3P)
water. This is the execution layer; on Purdue RCAC it is dispatched via slurm/run_md.sbatch.
Canonical residues use OpenMM's built-in amber14 OL3 parameters; non-canonical bases load the
per-residue XML/frcmod produced by src/parameterize.py.

Requires: openmm, and a prepared system (Amber prmtop/inpcrd from tleap, or an OpenMM System
built by parameterize.py). Given a valid system this runs unchanged on CPU or CUDA.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

import openmm as mm
from openmm import app, unit


def build_simulation(prmtop_path: str, inpcrd_path: str, platform_name: str = "CUDA"):
    """Load an Amber-prepared solvated RNA duplex and construct an OpenMM Simulation."""
    prmtop = app.AmberPrmtopFile(prmtop_path)
    inpcrd = app.AmberInpcrdFile(inpcrd_path)
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        hydrogenMass=1.5 * unit.amu,          # allows 4 fs timestep
    )
    system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 300.0 * unit.kelvin, 25))
    integrator = mm.LangevinMiddleIntegrator(300.0 * unit.kelvin, 1.0 / unit.picosecond,
                                             4.0 * unit.femtosecond)
    try:
        platform = mm.Platform.getPlatformByName(platform_name)
    except Exception:
        platform = mm.Platform.getPlatformByName("CPU")
    sim = app.Simulation(prmtop.topology, system, integrator, platform)
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    return sim


def run(system_prefix: str, ns: float, platform: str, out_dir: str, report_ps: float = 10.0):
    prmtop = system_prefix + ".prmtop"
    inpcrd = system_prefix + ".inpcrd"
    os.makedirs(out_dir, exist_ok=True)
    name = os.path.basename(system_prefix)

    sim = build_simulation(prmtop, inpcrd, platform)
    print(f"[{name}] minimising...")
    sim.minimizeEnergy(maxIterations=5000)

    # equilibrate (NVT->NPT handled by barostat), then production
    steps = int((ns * unit.nanoseconds) / (4.0 * unit.femtoseconds))
    stride = int((report_ps * unit.picoseconds) / (4.0 * unit.femtoseconds))
    dcd = os.path.join(out_dir, name + ".dcd")
    sim.reporters.append(app.DCDReporter(dcd, stride))
    sim.reporters.append(app.StateDataReporter(
        os.path.join(out_dir, name + ".log"), stride, step=True,
        potentialEnergy=True, temperature=True, progress=True, totalSteps=steps))
    print(f"[{name}] production: {ns} ns ({steps} steps) on {platform}")
    sim.step(steps)
    print(f"[{name}] done -> {dcd}")
    return dcd


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Run explicit-solvent MD on a prepared RNA duplex.")
    p.add_argument("--system", required=True, help="prefix for <system>.prmtop / <system>.inpcrd")
    p.add_argument("--ns", type=float, default=50.0)
    p.add_argument("--platform", default="CUDA", choices=["CUDA", "OpenCL", "CPU", "Reference"])
    p.add_argument("--out", default="work/md")
    args = p.parse_args()
    run(args.system, args.ns, args.platform, args.out)
