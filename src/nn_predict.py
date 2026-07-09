"""
nn_predict.py — Turner 1998 nearest-neighbor prediction of RNA duplex dG37 for the
validation anchors, so the ordering gate is pinned to real thermodynamics.

Parameters: Watson-Crick NN dG37 (kcal/mol), Xia et al. Biochemistry 1998, 37:14719
(as tabulated in the Turner NNDB, rna.urmc.rochester.edu/NNDB). Initiation, symmetry,
and terminal-AU penalties per the same model. These are the STANDARD Turner values;
exact figures should be confirmed against the NNDB before use as a hard gate.
"""
# NN step keyed by 5'->3' top strand dinucleotide (value = its reverse-complement equiv.)
NN = {"AA":-0.93,"UU":-0.93,"AU":-1.10,"UA":-1.33,
      "CU":-2.08,"AG":-2.08,"CA":-2.11,"UG":-2.11,
      "GU":-2.24,"AC":-2.24,"GA":-2.35,"UC":-2.35,
      "CG":-2.36,"GG":-3.26,"CC":-3.26,"GC":-3.42}
INIT=4.09; SYM=0.43; TERM_AU=0.45

def revcomp(s):
    c={"A":"U","U":"A","G":"C","C":"G"}; return "".join(c[b] for b in reversed(s))

def dG37(seq):
    dg=INIT
    for i in range(len(seq)-1):
        dg+=NN[seq[i:i+2]]
    if seq==revcomp(seq): dg+=SYM
    for end in (seq[0],seq[-1]):
        if end in "AU": dg+=TERM_AU
    return round(dg,2)

if __name__=="__main__":
    for s in ["GCGCGC","GCAUGC","AUAUAU","GACC","CAAU"]:
        print(f"{s:8s} dG37 = {dG37(s):6.2f} kcal/mol   (revcomp {revcomp(s)})")
