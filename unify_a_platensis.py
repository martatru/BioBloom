#!/usr/bin/env python3

import os

def main():
    # Ask for 4 input files
    input_files = []
    for i in range(4):
        path = input(f"Enter path for peptide list file {i+1}: ").strip()
        if not os.path.isfile(path):
            print(f"File not found: {path}")
            return
        input_files.append(path)

    # Output directory
    output_dir = "/home/marta/Desktop/biobloom_proteomic_data/aplatensis_unified_peptide_lists"
    os.makedirs(output_dir, exist_ok=True)  # tworzy katalog jeśli nie istnieje

    # Output file
    output_path = os.path.join(output_dir, "unified_peptides.txt")

    # Wczytaj wszystkie peptydy i usuń duplikaty
    peptides = set()
    for file_path in input_files:
        with open(file_path, "r") as f:
            for line in f:
                peptide = line.strip()
                if peptide:  # pomiń puste linie
                    peptides.add(peptide)

    # Sortowanie dla czytelności
    peptides = sorted(peptides)

    # Zapisz do pliku
    with open(output_path, "w") as out:
        out.write("\n".join(peptides))

    print(f"Unified list written to {output_path}")
    print(f"Total unique peptides: {len(peptides)}")

if __name__ == "__main__":
    main()
