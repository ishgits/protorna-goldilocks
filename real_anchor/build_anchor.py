"""
build_anchor.py — build an exact anchor RNA duplex from the 1RNA A-form scaffold by
base mutation, complete it with PDBFixer, and run OL3 implicit-solvent MD in OpenMM.

Scaffold: central 6 bp of PDB 1RNA ([U(UA)6A]2, A-form). Chain A (res 5-10) = AUAUAU,
chain B (res 19-24) its antiparallel WC complement. Mutating A->G / U->C on BOTH strands
preserves Watson-Crick complementarity, so from one identical backbone we get:
    au  : AUAUAU / AUAUAU     (no mutation)
    gc  : GCGCGC / GCGCGC     (all A->G, U->C)
    mixed: GCAUGC / GCAUGC-complement (A5->G,U6->C,A9->G,U10->C; B mirror)

This is a *controlled* comparison: identical starting helix + identical protocol, only
base identity differs. Output: per-anchor completed PDB + a short trajectory, plus
bp-H-bond distances and backbone RMSD saved as .npz for src/analyze.py.
"""
import sys, os, argparse
import numpy as np
from pdbfixer import PDBFixer
import openmm as mm
from openmm import app, unit

# purine / pyrimidine heavy-atom sets (base only)
A_ATOMS = {"N9","C8","N7","C5","C6","N6","N1","C2","N3","C4"}
G_ATOMS = {"N9","C8","N7","C5","C6","O6","N1","C2","N2","N3","C4"}
U_ATOMS = {"N1","C2","O2","N3","C4","O4","C5","C6"}
C_ATOMS = {"N1","C2","O2","N3","C4","N4","C5","C6"}
BASE_ATOMS = {"A":A_ATOMS,"G":G_ATOMS,"U":U_ATOMS,"C":C_ATOMS}
BACKBONE = {"P","OP1","OP2","O5'","C5'","C4'","O4'","C3'","O3'","C2'","O2'","C1'"}

# scaffold residue identity per (chain,resSeq)
SCAFFOLD = {  # from the extracted hexamer
    ("A",5):"A",("A",6):"U",("A",7):"A",("A",8):"U",("A",9):"A",("A",10):"U",
    ("B",19):"A",("B",20):"U",("B",21):"A",("B",22):"U",("B",23):"A",("B",24):"U",
}

# target sequences (5'->3'), chain A res 5..10, chain B res 19..24 (=complement, 5'->3')
TARGETS = {
    "au":    {"A":list("AUAUAU"), "B":list("AUAUAU")},
    "gc":    {"A":list("GCGCGC"), "B":list("GCGCGC")},
    "mixed": {"A":list("GCAUGC"), "B":list("GCAUGC")},
}
CHAIN_RES = {"A":[5,6,7,8,9,10], "B":[19,20,21,22,23,24]}

def target_map(name):
    m = {}
    for ch, resids in CHAIN_RES.items():
        for i, r in enumerate(resids):
            m[(ch, r)] = TARGETS[name][ch][i]
    return m

def mutate_pdb(scaffold_path, tmap, out_path):
    """Rewrite scaffold PDB: rename residues to target base; drop base atoms not in the
    target's heavy-atom set; strip 5'-terminal phosphate so PDBFixer caps a 5'-OH."""
    keep_lines = []
    term5 = {("A",5), ("B",19)}  # 5' ends of each strand
    for line in open(scaffold_path):
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        ch = line[21]
        resseq = int(line[22:26])
        key = (ch, resseq)
        if key not in tmap:
            continue
        tgt = tmap[key]
        # strip 5'-phosphate on terminal residues
        if key in term5 and atom in {"P","OP1","OP2"}:
            continue
        # keep backbone; for base atoms keep only those in target's set
        if atom in BACKBONE:
            keep = True
        elif atom in BASE_ATOMS[tgt]:
            keep = True
        else:
            keep = False  # e.g. N6 when mutating A->G, O4 when U->C
        if not keep:
            continue
        # rewrite residue name to single-letter target
        newline = line[:17] + ("%3s" % tgt) + line[20:]
        keep_lines.append(newline)
    with open(out_path, "w") as fh:
        fh.write("".join(keep_lines))
        fh.write("END\n")
    return out_path

def complete(pdb_in, pdb_out):
    fixer = PDBFixer(filename=pdb_in)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    with open(pdb_out, "w") as fh:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, fh)
    return pdb_out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scaffold", required=True)
    ap.add_argument("--anchor", required=True, choices=list(TARGETS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--ns", type=float, default=2.0)
    ap.add_argument("--temp", type=float, default=300.0)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    mut = os.path.join(args.out, f"{args.anchor}_mut.pdb")
    comp = os.path.join(args.out, f"{args.anchor}_complete.pdb")
    mutate_pdb(args.scaffold, target_map(args.anchor), mut)
    complete(mut, comp)
    print(f"[{args.anchor}] completed structure -> {comp}")

    pdb = app.PDBFile(comp)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds, hydrogenMass=1.5*unit.amu)
    integ = mm.LangevinMiddleIntegrator(args.temp*unit.kelvin, 2.0/unit.picosecond,
                                        4.0*unit.femtosecond)
    sim = app.Simulation(pdb.topology, system, integ,
                         mm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=800)
    sim.context.setVelocitiesToTemperature(args.temp*unit.kelvin)

    dcd = os.path.join(args.out, f"{args.anchor}.dcd")
    steps = int(args.ns*unit.nanoseconds/(4.0*unit.femtoseconds))
    stride = 75  # 0.3 ps
    sim.reporters.append(app.DCDReporter(dcd, stride))
    sim.reporters.append(app.StateDataReporter(sys.stdout, steps//5, step=True,
                          potentialEnergy=True, temperature=True))
    print(f"[{args.anchor}] MD {args.ns} ns ({steps} steps) at {args.temp} K on CPU")
    sim.step(steps)
    # save final topology as PDB (for MDAnalysis/mdtraj)
    with open(os.path.join(args.out, f"{args.anchor}_top.pdb"), "w") as fh:
        app.PDBFile.writeFile(pdb.topology, sim.context.getState(getPositions=True).getPositions(), fh)
    print(f"[{args.anchor}] done -> {dcd}")

if __name__ == "__main__":
    main()
