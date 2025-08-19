import sys
import os
import re

species = sys.argv[1]  # e.g., 'a_platensis', 't_suecica'
input_folder = "/home/marta/Desktop/biopep_results/batch_processing_results"
output_folder = "/home/marta/Desktop/biopep_results/peptide_lists"

os.makedirs(output_folder, exist_ok=True)

file_path = os.path.join(input_folder, f"{species}.csv")
print(f"Processing file: {file_path}")

try:
    with open(file_path, "r") as f:
        content = f.read()
except Exception as e:
    print(f"Failed to read {file_path}: {e}")
    sys.exit(1)

# extract peptides from "Results of enzyme action" section
# this regex captures the list of peptides, separated by spaces or newlines
match = re.search(r"Results of enzyme action\s*(.*?)\s*Location of released peptides", content, re.DOTALL)
if match:
    peptide_text = match.group(1)
    # split by spaces, strip empty entries
    peptides = [p.strip() for p in re.split(r"\s+", peptide_text) if p.strip()]
else:
    print("No peptide section found in the file.")
    peptides = []

# save peptide list
output_file = os.path.join(output_folder, f"{species}_peptides.txt")
with open(output_file, "w") as f:
    for pep in peptides:
        f.write(f"{pep}\n")

print(f"Saved peptide list to {output_file}")
