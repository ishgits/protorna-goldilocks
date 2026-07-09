"""
analyze_real.py — turn the real anchor trajectories into the repo's stability readouts.

For each anchor (au/mixed/gc) built on the identical 1RNA A-form backbone:
  * Watson-Crick pair geometry: purine-N1 ... pyrimidine-N3 distance per bp per frame
    -> fed to src/analyze.py base_pair_persistence + stability_score (the repo's own logic)
  * backbone RMSD vs frame 0 (mdtraj) -> rmsd_stability
  * MM-GBSA-style interstrand interaction energy:  E(duplex) - E(strandA) - E(strandB)
    on the OL3/GBn2 potential, single-trajectory, averaged over frames. This is the
    "MD-binding-energy" used as the Tm proxy; unlike fraying it separates
    3-H-bond G:C from 2-H-bond A:U with little sampling.
Outputs anchor_ordering.json and per-anchor npz.
"""
import os, sys, json
import numpy as np
import mdtraj as md
import openmm as mm
from openmm import app, unit

REPO = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(REPO, "src"))
import analyze as A  # the repo's tested ranking logic

ANCHORS = ["gc", "mixed", "au"]
WORK = "work"
PAIRS = [(0, 11), (1, 10), (2, 9), (3, 8), (4, 7), (5, 6)]  # chainA res k <-> chainB res 11-k
PURINE = {"A", "G"}

def wc_atom(top, res_index):
    """Return (atom_index, is_purine) for the WC central atom: purine N1 / pyrimidine N3."""
    res = list(top.residues)[res_index]
    is_pur = res.name in PURINE
    want = "N1" if is_pur else "N3"
    for a in res.atoms:
        if a.name == want:
            return a.index, is_pur
    raise KeyError(f"{want} not in residue {res.name}{res_index}")

def bp_distances(traj):
    top = traj.topology
    idx = []
    for i, j in PAIRS:
        ai, ip = wc_atom(top, i)
        aj, jp = wc_atom(top, j)
        # ensure purine-N1 ... pyrimidine-N3
        if ip and not jp:
            idx.append((ai, aj))
        elif jp and not ip:
            idx.append((aj, ai))
        else:  # both purine (shouldn't happen for WC) -> use N1..N1
            idx.append((ai, aj))
    idx = np.array(idx)
    d = md.compute_distances(traj, idx) * 10.0  # nm -> Angstrom
    return d  # (n_frames, n_pairs)

def interstrand_energy(pdb_path, traj):
    """MM-GBSA single-trajectory interaction energy (kcal/mol), mean over frames."""
    pdb = app.PDBFile(pdb_path)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    chains = list(pdb.topology.chains())
    resA = [r for r in chains[0].residues()]
    resB = [r for r in chains[1].residues()]

    def ctx(delete_res):
        mod = app.Modeller(pdb.topology, pdb.positions)
        mod.delete(delete_res)
        s = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff)
        integ = mm.VerletIntegrator(1*unit.femtosecond)
        c = mm.Context(s, integ, mm.Platform.getPlatformByName("CPU"))
        # map: indices of kept atoms in original topology order
        kept = [a.index for a in pdb.topology.atoms()
                if a.residue not in delete_res]
        return c, kept

    c_all, kept_all = ctx([])
    c_A, kept_A = ctx(resB)   # strand A only
    c_B, kept_B = ctx(resA)   # strand B only

    xyz = traj.xyz  # nm, (n_frames, n_atoms, 3)
    KJ2KCAL = 1.0 / 4.184
    dE = []
    for f in range(traj.n_frames):
        pos = xyz[f]
        def energy(c, kept):
            c.setPositions(pos[kept] * unit.nanometer)
            return c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
                unit.kilojoule_per_mole)
        e = energy(c_all, kept_all) - energy(c_A, kept_A) - energy(c_B, kept_B)
        dE.append(e * KJ2KCAL)
    return float(np.mean(dE)), float(np.std(dE))

def main():
    results = []
    summary = {}
    for a in ANCHORS:
        top = os.path.join(WORK, f"{a}_top.pdb")
        dcd = os.path.join(WORK, f"{a}.dcd")
        traj = md.load(dcd, top=top)
        bbnames = {"P", "O5'", "C5'", "C4'", "C3'", "O3'"}
        bb = np.array([a.index for a in traj.topology.atoms if a.name in bbnames])
        rmsd = md.rmsd(traj, traj, 0, atom_indices=bb) * 10.0  # Angstrom
        bpd = bp_distances(traj)
        dr = A.DuplexResult.from_traj(a, bpd, rmsd)
        results.append(dr)
        eint, estd = interstrand_energy(top, traj)
        summary[a] = dict(mean_bp_N1N3_ang=float(bpd.mean()),
                          bp_persistence=dr.mean_bp_persistence,
                          rmsd_stability=dr.rmsd_stability,
                          stability_score=dr.stability_score,
                          interstrand_interaction_kcal=eint,
                          interstrand_interaction_sd=estd,
                          n_frames=traj.n_frames)
        np.savez(os.path.join(WORK, f"{a}.npz"), bp_distances=bpd, rmsd=rmsd)

    ranked = A.rank_duplexes(results)
    expected = ["gc", "mixed", "au"]
    rep = A.report(results, expected=expected)
    rep["per_anchor"] = summary
    rep["energy_ranking_most_to_least_stable"] = sorted(
        summary, key=lambda k: summary[k]["interstrand_interaction_kcal"])
    with open(os.path.join(WORK, "anchor_ordering.json"), "w") as fh:
        json.dump(rep, fh, indent=2)
    print(json.dumps(rep, indent=2))

if __name__ == "__main__":
    main()
