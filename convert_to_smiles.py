import os
import subprocess

def main():
    # Ask the user for the input FASTA file
    input_fasta = input("Enter the path to the FASTA file with peptides: ").strip()
    if not os.path.exists(input_fasta):
        print(f"Error: file {input_fasta} does not exist.")
        return

    # Output folder
    output_dir = "/home/marta/Desktop/openbabel_output"
    os.makedirs(output_dir, exist_ok=True)

    # Output filename (same as input, but with .smi extension)
    base_name = os.path.splitext(os.path.basename(input_fasta))[0]
    output_file = os.path.join(output_dir, f"{base_name}.smi")

    try:
        # Call the p2smi CLI (fasta2smi)
        subprocess.run(
            ["fasta2smi", "-i", input_fasta, "-o", output_file],
            check=True
        )
        print(f"SMILES successfully saved to: {output_file}")
    except subprocess.CalledProcessError as e:
        print("Error running fasta2smi:", e)

if __name__ == "__main__":
    main()
