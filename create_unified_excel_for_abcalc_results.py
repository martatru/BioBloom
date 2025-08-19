import os
import pandas as pd
import glob
import argparse
from pathlib import Path

# === ARGPARSE ===
parser = argparse.ArgumentParser(description="Parse Biopep results for given species")
parser.add_argument(
    "species",
    nargs="+",
    help="Enter one or more species names to process"
)
args = parser.parse_args()

# === CONFIGURATION ===
root_folder = Path("/home/marta/Desktop/biopep_results")
output_folder = root_folder / "abc_evaluation"
output_folder.mkdir(parents=True, exist_ok=True)

# === LOOP THROUGH SPECIES ===
for species in args.species:
    species_path = root_folder / species
    data = []

    if species_path.is_dir():
        for protein_folder in os.listdir(species_path):
            protein_path = species_path / protein_folder
            if protein_path.is_dir():
                # Extract protein ID from folder name (between first and second '|')
                if '|' in protein_folder:
                    protein_id = protein_folder.split('|')[1]
                else:
                    protein_id = protein_folder

                # Find all Excel files in the folder
                excel_files = glob.glob(str(protein_path / '*.xlsx'))
                for excel_file in excel_files:
                    try:
                        df = pd.read_excel(excel_file, header=None)
                        sequence = df.iloc[0, 1]  # B1
                        A = df.iloc[2, 3]         # D3
                        B = df.iloc[2, 4]         # E3
                        data.append([protein_id, sequence, A, B])
                    except Exception as e:
                        print(f"[!] Error reading file {excel_file}: {e}")

        # Create DataFrame for species
        if data:
            df_species = pd.DataFrame(data, columns=['Protein_ID', 'Sequence', 'A', 'B'])
            excel_path = output_folder / f"abcalculations_results_{species}.xlsx"
            df_species.to_excel(excel_path, index=False)
            print(f"[OK] Data for {species} saved to file: {excel_path}")
        else:
            print(f"[!] No data to save for species: {species}")
    else:
        print(f"[!] No folder found for species: {species}")
