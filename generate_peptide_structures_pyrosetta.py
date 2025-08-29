#!/usr/bin/env python3
import argparse, os, re
from pyrosetta import init, pose_from_sequence

"""

python generate_peptide_structures_pyrosetta.py /home/marta/Desktop/supplementary_data/list_of_peptides_selected_for_docking.fasta -o pdbs --silent

"""

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")  # standard L-amino acids (1-letter)

def read_fasta(path):
    entries = []
    header = None
    seq_lines = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_lines)))
                header = line[1:].strip() or "seq"
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        entries.append((header, "".join(seq_lines)))
    # If file had no header, treat whole file as a single sequence
    if not entries:
        with open(path, "r") as f:
            raw = re.sub(r"\s+", "", f.read())
        entries = [("seq", raw)]
    return entries

def sanitize_seq_id(s):
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s[:60] if s else "seq"

def validate_seq(seq, name):
    seq = seq.upper().replace("U", "C")  # simple safeguard; edit if you need Sec
    bad = [c for c in seq if c not in VALID_AA]
    if bad:
        uniq = sorted(set(bad))
        print(f"[WARN] {name}: non-standard letters found {uniq}. "
              "They will be kept as-is and may fail in Rosetta.")
    return seq

def main():
    ap = argparse.ArgumentParser(description="Convert FASTA peptides to PDB (PyRosetta).")
    ap.add_argument("fasta", help="Input FASTA file (multi-FASTA supported).")
    ap.add_argument("-o", "--outdir", default="out_pdbs", help="Output directory.")
    ap.add_argument("--minlen", type=int, default=2, help="Minimum length to keep.")
    ap.add_argument("--maxlen", type=int, default=50, help="Maximum length to keep.")
    ap.add_argument("--silent", action="store_true", help="Mute Rosetta output.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Initialize Rosetta
    opts = "-mute all" if args.silent else ""
    init(opts)

    entries = read_fasta(args.fasta)
    if not entries:
        print("[ERR] No sequences found.")
        return

    n_ok = 0
    for raw_name, raw_seq in entries:
        name = sanitize_seq_id(raw_name)
        seq = validate_seq(raw_seq, name)
        if not seq:
            print(f"[SKIP] {name}: empty sequence.")
            continue
        if not (args.minlen <= len(seq) <= args.maxlen):
            print(f"[SKIP] {name}: length {len(seq)} outside [{args.minlen},{args.maxlen}].")
            continue

        try:
            # Build an extended peptide pose with the standard full-atom residue type set
            pose = pose_from_sequence(seq, "fa_standard")
            out_path = os.path.join(args.outdir, f"{name}.pdb")
            pose.dump_pdb(out_path)
            print(f"[OK]  {name}: {len(seq)} aa -> {out_path}")
            n_ok += 1
        except Exception as e:
            print(f"[FAIL] {name}: {e}")

    print(f"\nDone. Wrote {n_ok} PDB(s) to {args.outdir}.")

if __name__ == "__main__":
    main()
