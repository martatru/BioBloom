from pymol import cmd
import sys


'''

mkdir -p /home/marta/Pulpit/pliki_do_dokowania_niezjebane/

REF=/home/marta/Pulpit/ACE_niezjebane/ACE/1O86.pdb
REC=/home/marta/Pulpit/ACE_niezjebane/ACE/ACE_repacked.pdb
PEPDIR=/home/marta/Pulpit/ACE_niezjebane/PEPTYDY/
OUTDIR=/home/marta/Pulpit/pliki_do_dokowania_niezjebane/

for pep in "$PEPDIR"/*.pdb; do
  base=$(basename "$pep" .pdb)
  /usr/bin/pymol -cq place_pep_into_ace.py -- "$REF" "$REC" "$pep" "$OUTDIR/ACE_${base}.pdb"


'''

if len(sys.argv) < 5:
    print("Usage: pymol -cq place_pep_into_ace.py -- REF_with_LPR.pdb REC_repacked.pdb PEP.pdb OUT.pdb")
    sys.exit(1)

ref_with_lpr, receptor_pdb, pep_pdb, out_pdb = sys.argv[1:5]

cmd.load(ref_with_lpr, "ref")       # 1O86 z LPR (tylko jako referencja)
cmd.load(receptor_pdb, "rec")       # czysty ACE
cmd.load(pep_pdb, "pep")            # jeden peptyd

cmd.alter("rec and polymer.protein", 'chain="A"')
cmd.alter("pep and polymer.protein", 'chain="P"')
cmd.sort()

if cmd.count_atoms("ref and resn LPR") > 0:
    cx, cy, cz = cmd.centerofmass("ref and resn LPR")
else:
    cx, cy, cz = cmd.centerofmass("rec")

px, py, pz = cmd.centerofmass("pep")
cmd.translate([cx - px, cy - py, cz - pz], "pep")

cmd.save(out_pdb, "rec or pep")