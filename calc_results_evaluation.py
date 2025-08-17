import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob

# Input root folder with species subfolders
root_folder = '/home/marta/Desktop/biopep_results'

# Output folder
output_folder = '/home/marta/Desktop/biopep_results/abc_evaluation'
os.makedirs(output_folder, exist_ok=True)

# Choose species to analyze
species_to_analyze = 'tetraselmis_suecica'  # put the species folder name here

# Collect all data
data = []

species_path = os.path.join(root_folder, species_to_analyze)
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
                    A = df.iloc[2, 3]  # D3
                    B = df.iloc[2, 4]  # E3
                    data.append([protein_id, A, B])
                except Exception as e:
                    print(f"Error reading {excel_file}: {e}")

# Create dataframe
df_species = pd.DataFrame(data, columns=['Protein_ID', 'A', 'B'])

# Basic statistics
stats = df_species.agg({'A':['mean','median','std','min','max','count'],
                        'B':['mean','median','std','min','max','count']})
stats.loc['CV'] = stats.loc['std'] / stats.loc['mean']

# Percent of proteins with low B (< 0.1)
low_B_threshold = 0.1
percent_low_B = (df_species['B'] < low_B_threshold).mean() * 100

# Top N proteins with low B and high A
N = 5
top_proteins = df_species.sort_values(['B', 'A'], ascending=[True, False]).head(N)

# Save statistics and top proteins to Excel
excel_path = os.path.join(output_folder, f'{species_to_analyze}_ACEI_analysis.xlsx')
with pd.ExcelWriter(excel_path) as writer:
    stats.to_excel(writer, sheet_name='Basic_Statistics')
    pd.DataFrame({'Percent_Low_B':[percent_low_B]}).to_excel(writer, sheet_name='Percent_Low_B')
    top_proteins.to_excel(writer, sheet_name='Top_Proteins', index=False)

# Plots
plt.figure(figsize=(10,5))
sns.histplot(df_species['A'], kde=True, bins=20)
plt.title(f"Histogram of A for {species_to_analyze}")
plt.xlabel("A (frequency of ACE-I fragments)")
plt.ylabel("Count")
plt.savefig(os.path.join(output_folder, f'{species_to_analyze}_Histogram_A.png'))
plt.close()

plt.figure(figsize=(10,5))
sns.histplot(df_species['B'], kde=True, bins=20)
plt.title(f"Histogram of B for {species_to_analyze}")
plt.xlabel("B (potential ACE-I activity, µM⁻¹)")
plt.ylabel("Count")
plt.savefig(os.path.join(output_folder, f'{species_to_analyze}_Histogram_B.png'))
plt.close()

plt.figure(figsize=(10,5))
sns.boxplot(y=df_species['A'])
plt.title(f"Boxplot of A for {species_to_analyze}")
plt.ylabel("A (frequency of ACE-I fragments)")
plt.savefig(os.path.join(output_folder, f'{species_to_analyze}_Boxplot_A.png'))
plt.close()

plt.figure(figsize=(10,5))
sns.boxplot(y=df_species['B'])
plt.title(f"Boxplot of B for {species_to_analyze}")
plt.ylabel("B (potential ACE-I activity, µM⁻¹)")
plt.savefig(os.path.join(output_folder, f'{species_to_analyze}_Boxplot_B.png'))
plt.close()

# Correlation A vs B
correlation = df_species[['A','B']].corr()
correlation.to_excel(os.path.join(output_folder, f'{species_to_analyze}_Correlation_A_vs_B.xlsx'))
