import os
import sys
import pandas as pd

# === Command-line argument ===
if len(sys.argv) != 2:
    print("Usage: python file.py <species_name>")
    sys.exit(1)

species = sys.argv[1].strip()

# === Paths ===
input_folder = '/home/marta/Desktop/biopep_results/abc_evaluation'
input_file = os.path.join(input_folder, f"abcalculations_results_{species}.xlsx")
output_file = os.path.join(input_folder, f"analysis_results_{species}.txt")

if not os.path.exists(input_file):
    print(f"[!] No file found for species: {species}")
    sys.exit(1)

# === Load data ===
df_species = pd.read_excel(input_file)

# === Open output file ===
with open(output_file, 'w') as out_f:

    out_f.write(f"\n=== Analysis for {species} ===\n")

    # Basic statistics
    stats = df_species.agg({'A':['mean','median','std','min','max','count'],
                            'B':['mean','median','std','min','max','count']})
    stats.loc['CV'] = stats.loc['std'] / stats.loc['mean']
    out_f.write("\nBasic Statistics:\n")
    out_f.write(stats.to_string())
    out_f.write("\n")

    # Percent of proteins with low B (<0.1)
    low_B_threshold = 0.1
    percent_low_B = (df_species['B'] < low_B_threshold).mean() * 100
    out_f.write(f"\nPercent of proteins with B < {low_B_threshold}: {percent_low_B:.2f}%\n")

    # Compute A/B ratio and select top 10 proteins
    df_species['Score'] = df_species['A'] / df_species['B']
    top_proteins = df_species.sort_values('Score', ascending=False).head(10)

    out_f.write("\nTop 10 proteins (highest A/B ratio):\n")
    out_f.write(top_proteins.to_string(index=False))
    out_f.write("\n" + "="*50 + "\n")

print(f"[OK] Analysis results saved to: {output_file}")
