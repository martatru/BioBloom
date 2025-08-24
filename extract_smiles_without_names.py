import os

# ask user for input file path
inp = input("Enter path to input file: ").strip()

# prepare output path
out_dir = "/home/marta/Desktop/openbabel_output"
os.makedirs(out_dir, exist_ok=True)
base = os.path.splitext(os.path.basename(inp))[0]
out = os.path.join(out_dir, f"{base}2.smi")

n_ok = n_bad = 0
with open(inp, "r", encoding="utf-8") as fin, open(out, "w", encoding="utf-8") as fout:
    for line in fin:
        line = line.strip()
        if not line:
            continue
        # take only the part after colon (SMILES)
        if ":" in line:
            smiles = line.split(":", 1)[1].strip()
        else:
            smiles = line  # fallback if line has only SMILES
        
        if smiles:
            fout.write(smiles + "\n")
            n_ok += 1
        else:
            n_bad += 1

print(f"Extracted {n_ok} SMILES and saved to: {out}")
if n_bad:
    print(f"Skipped {n_bad} invalid/empty lines.")
