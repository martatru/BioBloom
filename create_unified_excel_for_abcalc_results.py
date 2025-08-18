import os
import pandas as pd
import glob

# Input root folder with species subfolders
root_folder = '/home/marta/Desktop/biopep_results'

# Output folder
output_folder = '/home/marta/Desktop/biopep_results/abc_evaluation'
os.makedirs(output_folder, exist_ok=True)

# Species to analyze
species_list = ["arthrospira_platensis", "microchloropsis_salina", "tetraselmis_suecica"]

for species in species_list:
    data = []
    species_path = os.path.join(root_folder, species)

    if os.path.isdir(species_path):
        for protein_folder in os.listdir(species_path):
            protein_path = os.path.join(species_path, protein_folder)
            if os.path.isdir(protein_path):
                # Extract protein ID from folder name (between first and second '|')
                if '|' in protein_folder:
                    protein_id = protein_folder.split('|')[1]
                else:
                    protein_id = protein_folder

                # Find all Excel files in the folder
                excel_files = glob.glob(os.path.join(protein_path, '*.xlsx'))
                for excel_file in excel_files:
                    try:
                        df = pd.read_excel(excel_file, header=None)
                        sequence = df.iloc[0, 1]  # B1
                        A = df.iloc[2, 3]         # D3
                        B = df.iloc[2, 4]         # E3
                        data.append([protein_id, sequence, A, B])
                    except Exception as e:
                        print(f"Error reading {excel_file}: {e}")

        # Create dataframe for species
        df_species = pd.DataFrame(data, columns=['Protein_ID', 'Sequence', 'A', 'B'])

        # Save to one Excel file per species
        excel_path = os.path.join(output_folder, f"abcalculations_results_{species}.xlsx")
        df_species.to_excel(excel_path, index=False)
        print(f"[OK] Dane dla {species} zapisane w pliku: {excel_path}")
    else:
        print(f"[!] Brak folderu dla gatunku: {species}")
