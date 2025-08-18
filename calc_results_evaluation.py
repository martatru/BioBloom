import os
import pandas as pd

# Paths
input_folder = '/home/marta/Desktop/biopep_results/abc_evaluation'
species_list = ["arthrospira_platensis", "microchloropsis_salina", "tetraselmis_suecica"]

for species in species_list:
    input_file = os.path.join(input_folder, f"abcalculations_results_{species}.xlsx")
    df_species = pd.read_excel(input_file)

    print(f"\n=== Analysis for {species} ===")

    # Basic statistics
    stats = df_species.agg({'A':['mean','median','std','min','max','count'],
                            'B':['mean','median','std','min','max','count']})
    stats.loc['CV'] = stats.loc['std'] / stats.loc['mean']
    print("\nBasic Statistics:")
    print(stats)

    # Percent of proteins with low B (<0.1)
    low_B_threshold = 0.1
    percent_low_B = (df_species['B'] < low_B_threshold).mean() * 100
    print(f"\nPercent of proteins with B < {low_B_threshold}: {percent_low_B:.2f}%")

    # Best protein: highest A + lowest B
    best_protein = df_species.sort_values(['B','A'], ascending=[True, False]).head(1)
    print("\nBest protein (highest A + lowest B):")
    print(best_protein)
