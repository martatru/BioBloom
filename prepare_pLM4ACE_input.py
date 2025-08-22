import os

# Ask user for input file path
input_file = input("Enter the full path to your peptide .txt file: ").strip()

# Check if the file exists
if not os.path.isfile(input_file):
    print("File not found. Please check the path and try again.")
    exit()

# Extract the filename without extension to use in output file
species_name = os.path.splitext(os.path.basename(input_file))[0]

# Path to Desktop folder
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
output_folder = os.path.join(desktop_path, "pLM4ACE_input")

# Create the folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Output file path with species name included
output_file = os.path.join(output_folder, f"{species_name}_formatted.txt")

# Read sequences from the input file
with open(input_file, "r") as f:
    lines = f.readlines()

# Format sequences with '>' on a separate line
formatted_lines = []
for line in lines:
    seq = line.strip()
    if seq:  # skip empty lines
        formatted_lines.append(">\n")
        formatted_lines.append(seq + "\n")

# Write to the output file
with open(output_file, "w") as f:
    f.writelines(formatted_lines)

print(f"Formatted peptides saved to {output_file}")
