import pandas as pd
from pathlib import Path

# Get folder path from the user
folder_path = input("Enter the path to the folder with .xlsx files: ").strip()
folder = Path(folder_path)

if not folder.is_dir():
    print("The provided path is not a folder.")
    exit()

# Create a list to store high-activity sequences
high_activity_sequences = []

# Iterate over all .xlsx files in the folder
for file_path in folder.glob("*.xlsx"):
    try:
        df = pd.read_excel(file_path)
        # Check if required columns exist
        if 'sequence' in df.columns and 'activity' in df.columns:
            high_seqs = df[df['activity'].str.contains("high activity", case=False, na=False)]['sequence'].tolist()
            high_activity_sequences.extend(high_seqs)
        else:
            print(f"File {file_path.name} does not contain the required columns.")
    except Exception as e:
        print(f"Error reading file {file_path.name}: {e}")

# Check if any sequences were found
if high_activity_sequences:
    # Create output folder
    output_folder = Path("/home/marta/Desktop/pLM4ACE_results/txt_files_summarized")
    output_folder.mkdir(parents=True, exist_ok=True)
    
    # Set output file name based on the folder name
    output_file = output_folder / f"{folder.name}.txt"
    
    # Write sequences to the file
    with open(output_file, "w") as f:
        for seq in high_activity_sequences:
            f.write(seq + "\n")
    
    print(f"Saved {len(high_activity_sequences)} sequences to file: {output_file}")
else:
    print("No high-activity sequences were found.")
