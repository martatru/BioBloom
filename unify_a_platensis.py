#!/usr/bin/env python3

def main():
    import os

    # Ask for 4 input files
    input_files = []
    for i in range(4):
        path = input(f"Enter path for peptide list file {i+1}: ").strip()
        if not os.path.isfile(path):
            print(f"File not found: {path}")
            return
        input_files.append(path)

    # Output file
    output_path = "/home/marta/Desktop/biobloom_proteomic_data/aplatensis_unified_peptide_lists"

    peptides = set()
    for file_path in input_files:
        with open(file_path, "r") as f:
            for line in f:
                peptide = line.strip()
                if peptide:  # skip empty lines
                    peptides.add(peptide)

    # Sort for reproducibility
    peptides = sorted(peptides)

    # Write unified list
    with open(output_path, "w") as out:
        out.write("\n".join(peptides))

    print(f"Unified list written to {output_path}")
    print(f"Total unique peptides: {len(peptides)}")


if __name__ == "__main__":
    main()
