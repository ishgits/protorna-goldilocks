#!/bin/bash
# ============================================================
# build_system.sh - build solvated, ionised A-RNA duplex systems for GROMACS
#   nab (A-RNA duplex) -> tleap (OL3 + TIP3P + ions) -> acpype (-> GROMACS .top/.gro)
#
# RUN ONCE, on a machine with AmberTools (nab, tleap) + acpype available.
# On Negishi this is a login-node / interactive step, NOT a GPU job. Output lands in ./work/
# and is then consumed by the cluster_jobs/<cluster>/slurm_*.sh production scripts.
#
# This is the ONE step that must run in an AmberTools environment. Everything downstream
# (EM/NVT/NPT/production/analysis) is pure GROMACS on your proven Negishi setup.
# ============================================================
set -euo pipefail

# --- environment (edit to your setup) ---------------------------------------------------------
# e.g.  module load anaconda && conda activate ambertools   # must provide nab, tleap, acpype
#       (acpype: pip install acpype  /  conda install -c conda-forge acpype)
command -v nab    >/dev/null || { echo "nab not found (AmberTools)"; exit 1; }
command -v tleap  >/dev/null || { echo "tleap not found (AmberTools)"; exit 1; }
command -v acpype >/dev/null || { echo "acpype not found (pip/conda install acpype)"; exit 1; }

mkdir -p work && cd work

# --- anchor duplexes: id -> strand-1 sequence (5'->3'). All self-complementary. -------------
declare -A SEQS=( [gc_anchor]=gcgcgc [mixed_anchor]=gcaugc [au_anchor]=auauau )
# Non-canonical generalisation anchor (Benner Z-P) is added once Z/P residues are prepared.

revcomp () { echo "$1" | rev | tr 'acgu' 'ugca'; }   # RNA reverse complement

for D in "${!SEQS[@]}"; do
  S1=${SEQS[$D]}; S2=$(revcomp "$S1")
  echo "=== building $D : 5'-${S1}-3' / 5'-${S2}-3' ==="

  # 1) ideal A-RNA Watson-Crick duplex via nab's wc_helix
  #    NOTE: verify the wc_helix signature against your AmberTools version (AmberTools/nab manual,
  #    'Building Nucleic Acids'). Args below: (seq1,"",type1, seq2,"",type2, xoff, incl, twist, rise, opts)
  cat > ${D}.nab <<EOF
molecule m;
m = wc_helix( "${S1}", "", "arna", "${S2}", "", "arna", 2.25, -4.96, 32.7, 2.81, "" );
putpdb( "${D}_raw.pdb", m );
EOF
  nab ${D}.nab -o ${D}_nab
  ./${D}_nab
  # Fallback if wc_helix is unavailable: drop in a canonical A-RNA duplex PDB as ${D}_raw.pdb.

  # 2) OL3 force field + TIP3P solvation + neutralise + ~0.1 M NaCl
  cat > ${D}.leap <<EOF
source leaprc.RNA.OL3
source leaprc.water.tip3p
m = loadpdb ${D}_raw.pdb
solvatebox m TIP3PBOX 12.0
addions m Na+ 0                 # neutralise the backbone charge
addions m Na+ 14 Cl- 14        # ~0.1 M NaCl (adjust count to your box volume)
saveamberparm m ${D}.prmtop ${D}.inpcrd
savepdb m ${D}_solv.pdb
quit
EOF
  tleap -f ${D}.leap

  # 3) AMBER -> GROMACS
  acpype -p ${D}.prmtop -x ${D}.inpcrd -b ${D} -o gmx
  cp ${D}.amb2gmx/${D}_GMX.gro  ${D}.gro
  cp ${D}.amb2gmx/${D}_GMX.top  topol_${D}.top
  # posre.itp for -DPOSRES: acpype writes posre_${D}.itp; ensure topol includes it under
  # #ifdef POSRES. (acpype adds this automatically for the solute.)
  echo "   -> work/${D}.gro , work/topol_${D}.top"
done

echo "All anchor systems built. Next: submit the appropriate script under cluster_jobs/<cluster>/ (e.g. sbatch cluster_jobs/negishi/slurm_em_nvt.sh)"
