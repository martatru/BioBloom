import os
from odf.opendocument import load
from odf.table import Table, TableRow, TableCell
from odf.text import P

def read_ods_cells(path):
    """Returns all cells from an .ods file as a list of lists of strings"""
    doc = load(path)
    texts = []
    for table in doc.getElementsByType(Table):
        for row in table.getElementsByType(TableRow):
            row_texts = []
            for cell in row.getElementsByType(TableCell):
                cell_text = "".join(
                    t.data for t in cell.getElementsByType(P)
                    for t in t.childNodes if hasattr(t, 'data')
                )
                row_texts.append(cell_text.strip())
            if row_texts:
                texts.append(row_texts)
    return texts

def extract_ace_peptides_from_file(filepath, activity_col="MNO", sequence_col="EF"):
    """Extracts only peptides with ACE inhibitory activity"""
    cells = read_ods_cells(filepath)
    ace_peptides = []

    # Convert column letters to indices (A=0, B=1, ..., Z=25, AA=26, AB=27, etc.)
    def col_to_index(col):
        index = 0
        for c in col:
            index = index * 26 + (ord(c.upper()) - ord('A') + 1)
        return index - 1

    activity_idx = col_to_index(activity_col)
    sequence_idx = col_to_index(sequence_col)

    for row in cells:
        # Skip rows that are too short
        if len(row) <= max(activity_idx, sequence_idx):
            continue
        activity = row[activity_idx].strip().lower()
        sequence = row[sequence_idx].strip()
        if "ace inhibitor" in activity and sequence:
            ace_peptides.append(sequence)
    
    return ace_peptides

def main():
    species = input("Enter the species name (species_name): ").strip()
    
    input_folder = "/home/marta/Desktop/biopep_results/batch_processing_results/odsy/"
    output_folder = "/home/marta/Desktop/biopep_results/peptide_lists/peptides_search_for_bioac_frag_list"
    os.makedirs(output_folder, exist_ok=True)
    
    all_peptides = []
    
    # Select only .ods files that contain the species name
    files = [file for file in os.listdir(input_folder)
             if file.endswith(".ods") and species in file]
    
    for file in files:
        filepath = os.path.join(input_folder, file)
        peptides = extract_ace_peptides_from_file(filepath)
        if peptides:
            all_peptides.extend(peptides)
    
    outpath = os.path.join(output_folder, f"{species}_ace_inhibitors.txt")
    with open(outpath, "w") as f:
        f.write("\n".join(all_peptides))
    
    print(f"Saved {len(all_peptides)} ACE inhibitory peptides to: {outpath}")

if __name__ == "__main__":
    main()
