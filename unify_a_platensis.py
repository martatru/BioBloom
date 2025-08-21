"""
Merging separated a. platensis peptide files and deleting peptide duplicates script

"""

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
    os.makedirs(output_dir, exist_ok=True)  # create directory if it doesn't exist

    # Output file
    output_path = os.path.join(output_dir, "unified_peptides.txt")

    # Load all peptides and remove duplicates
    peptides = set()
    for file_path in input_files:
        with open(file_path, "r") as f:
            for line in f:
                peptide = line.strip()
                if peptide:  # skip empty lines
                    peptides.add(peptide)

    # Sort for readability
    peptides = sorted(peptides)

    # Write to file
    with open(output_path, "w") as out:
        out.write("\n".join(peptides))

    print(f"Unified list written to {output_path}")
    print(f"Total unique peptides: {len(peptides)}")

if __name__ == "__main__":
    main()
