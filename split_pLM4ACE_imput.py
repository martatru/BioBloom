import os

def split_fasta(input_path, output_dir, num_parts):
    # Read the file into peptide pairs
    with open(input_path, "r") as f:
        lines = f.readlines()

    # Group lines into peptide entries (header + sequence)
    entries = []
    for i in range(0, len(lines), 2):
        header = lines[i].strip()
        seq = lines[i + 1].strip()
        entries.append((header, seq))

    total_entries = len(entries)
    part_size = total_entries // num_parts
    remainder = total_entries % num_parts

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    start = 0
    for i in range(num_parts):
        # Distribute remainder (so some files may have one extra entry)
        end = start + part_size + (1 if i < remainder else 0)
        part_entries = entries[start:end]
        start = end

        output_path = os.path.join(output_dir, f"split_part_{i+1}.txt")
        with open(output_path, "w") as out_f:
            for header, seq in part_entries:
                out_f.write(f"{header}\n{seq}\n")

        print(f"Created: {output_path} with {len(part_entries)} entries")

if __name__ == "__main__":
    input_path = input("Enter the path to the input .txt file: ").strip()
    num_parts = int(input("Enter the number of parts to split into: ").strip())
    output_dir = "/home/marta/Desktop/pLM4ACE_input/separated_input"
    split_fasta(input_path, output_dir, num_parts)
