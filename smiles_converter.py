"""
this script converts fasta-->smiles or smiles -->fasta and saves the output into a parsed folder

"""
import os, sys, csv, subprocess
from typing import Optional, Tuple

def choose_direction():
    print("Choose conversion:")
    print("  1) FASTA → SMILES (p2smi/fasta2smi)")
    print("  2) SMILES → FASTA (RDKit)")
    c = input("Enter 1 or 2 (or f2s/s2f): ").strip().lower()
    if c in ("1","f2s","fasta>smiles","fasta->smiles"): return "f2s"
    if c in ("2","s2f","smiles>fasta","smiles->fasta"): return "s2f"
    sys.exit("Invalid choice.")

def ask_existing_path(prompt):
    p = input(prompt).strip()
    if not os.path.exists(p):
        sys.exit(f"Error: path does not exist: {p}")
    return p

def ask_outdir():
    outdir = input("Enter the output folder: ").strip()
    if not outdir:
        sys.exit("Error: output folder cannot be empty.")
    os.makedirs(outdir, exist_ok=True)
    return outdir

# ---------- FASTA→SMILES using p2smi CLI ----------
def fasta_to_smiles(fasta_path: str, out_smi: str) -> None:
    cmd = ["fasta2smi", "-i", fasta_path, "-o", out_smi]
    print("Running:", " ".join(cmd))
    res = subprocess.run(cmd)
    if res.returncode != 0:
        sys.exit(f"fasta2smi failed with code {res.returncode}")
    print(f"OK: saved → {out_smi}")

# ---------- SMILES→FASTA using RDKit ----------
def load_rdkit():
    try:
        from rdkit import Chem
        return Chem
    except Exception as e:
        sys.exit("Error: RDKit not available in this environment.")

# Amino-acid residue SMARTS/SMILES patterns for canonical AAs
# (taken from the RDKit discussion prototype; order is important)  :contentReference[oaicite:1]{index=1}
AA_SMILES = {
    'ALA': 'C[C@H](N)C=O',
    'CYS': 'N[C@H](C=O)CS',
    'ASP': 'N[C@H](C=O)CC(=O)O',
    'GLU': 'N[C@H](C=O)CCC(=O)O',
    'PHE': 'N[C@H](C=O)Cc1ccccc1',
    'GLY': 'NCC=O',
    'HIS': 'N[C@H](C=O)Cc1c[nH]cn1',
    'ILE': 'CC[C@H](C)[C@H](N)C=O',
    'LYS': 'NCCCC[C@H](N)C=O',
    'LEU': 'CC(C)C[C@H](N)C=O',
    'MET': 'CSCC[C@H](N)C=O',
    'ASN': 'NC(=O)C[C@H](N)C=O',
    'PRO': 'O=C[C@@H]1CCCN1',
    'GLN': 'NC(=O)CC[C@H](N)C=O',
    'ARG': 'N=C(N)NCCC[C@H](N)C=O',
    'SER': 'N[C@H](C=O)CO',
    'THR': 'C[C@@H](O)[C@H](N)C=O',
    'VAL': 'CC(C)[C@H](N)C=O',
    'TRP': 'N[C@H](C=O)Cc1c[nH]c2ccccc12',
    'TYR': 'N[C@H](C=O)Cc1ccc(O)cc1',
}
AA_ORDER = ['GLY','ALA','VAL','CYS','ASP','GLU','PHE','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','TRP','TYR']

def mol_to_sequence_via_rdkit(Chem, mol) -> Optional[str]:
    if mol is None: return None
    # Mark peptide C-alpha atoms so RDKit can build a sequence later
    CA_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[C:0](=[O:1])[C:2][N:3]"))
    for atoms in CA_matches:
        a = mol.GetAtomWithIdx(atoms[2])
        info = Chem.AtomPDBResidueInfo()
        info.SetName(" CA ")  # spaces are important
        a.SetMonomerInfo(info)

    # Assign residue names to marked CA atoms by matching AA templates
    for aa in AA_ORDER:
        tmpl = Chem.MolFromSmiles(AA_SMILES[aa])
        for match in mol.GetSubstructMatches(tmpl):
            for idx in match:
                a = mol.GetAtomWithIdx(idx)
                mi = a.GetMonomerInfo()
                if mi is not None and mi.GetName() == " CA ":
                    info = Chem.AtomPDBResidueInfo()
                    info.SetName(" CA "); info.SetResidueName(aa)
                    a.SetMonomerInfo(info)

    # Renumber backbone atoms to get correct order, then extract sequence
    bb = "O" + "C(=O)CN"*len(mol.GetSubstructMatches(Chem.MolFromSmiles(AA_SMILES["GLY"])))
    bb_mol = Chem.MolFromSmiles(bb)
    bb_match = mol.GetSubstructMatches(bb_mol)
    if not bb_match:
        return None
    id_list = list(bb_match[0])[::-1]  # reverse
    all_ids = [a.GetIdx() for a in mol.GetAtoms()]
    for idx in all_ids:
        if idx not in id_list:
            id_list.append(idx)
    m2 = Chem.RenumberAtoms(mol, newOrder=id_list)
    try:
        seq = Chem.MolToSequence(m2)
        return seq if seq else None
    except Exception:
        return None

def smi_line_parser(line: str) -> Tuple[str, str]:
    # Typical .smi: "SMILES<space>name", tolerate tabs
    parts = line.strip().split()
    if not parts: return ("","")
    smi = parts[0]
    name = parts[1] if len(parts) > 1 else ""
    return (smi, name)

def smiles_to_fasta(smi_path: str, out_fasta: str) -> None:
    Chem = load_rdkit()
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')

    ok, fail = 0, 0
    with open(smi_path, "r", encoding="utf-8", errors="ignore") as fin, \
         open(out_fasta, "w", encoding="utf-8") as fout:
        for i, line in enumerate(fin, 1):
            if not line.strip(): continue
            smi, name = smi_line_parser(line)
            if not smi: continue
            mol = Chem.MolFromSmiles(smi)
            seq = mol_to_sequence_via_rdkit(Chem, mol)
            if seq:
                header = name if name else f"mol_{i}"
                fout.write(f">{header}\n{seq}\n")
                ok += 1
            else:
                fail += 1
                print(f"[warn] line {i}: could not derive sequence", file=sys.stderr)

    print(f"Done. Wrote {ok} sequences to {out_fasta}. Failed: {fail}")

def main():
    mode = choose_direction()
    if mode == "f2s":
        fasta = ask_existing_path("Enter the path to the FASTA file with peptides: ")
        outdir = ask_outdir()
        base = os.path.splitext(os.path.basename(fasta))[0]
        out_smi = os.path.join(outdir, base + ".smi")
        if os.path.exists(out_smi):
            ans = input(f"Output exists ({out_smi}). Overwrite? [y/N]: ").strip().lower()
            if ans not in ("y","yes"): sys.exit("Aborted.")
        fasta_to_smiles(fasta, out_smi)
    else:
        smi = ask_existing_path("Enter the path to the SMILES (.smi) file: ")
        outdir = ask_outdir()
        base = os.path.splitext(os.path.basename(smi))[0]
        out_fa = os.path.join(outdir, base + ".fasta")
        if os.path.exists(out_fa):
            ans = input(f"Output exists ({out_fa}). Overwrite? [y/N]: ").strip().lower()
            if ans not in ("y","yes"): sys.exit("Aborted.")
        smiles_to_fasta(smi, out_fa)

if __name__ == "__main__":
    main()
