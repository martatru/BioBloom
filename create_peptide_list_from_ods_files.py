import os
from odf.opendocument import load
from odf.table import Table, TableRow, TableCell
from odf.text import P

#t_suecica
#a_platensis
#m_salina

def read_ods_cells(path):
    """Returns all cells from an .ods file as a list of strings"""
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

def extract_peptides_from_file(filepath):
    """Extracts all peptides under 'Results of enzyme action' sections"""
    cells = read_ods_cells(filepath)
    peptides = []
    i = 0
    while i < len(cells):
        row = cells[i]
        for cell in row:
            if cell.strip().lower() == "results of enzyme action":
                # search for the first non-empty cell in the following rows
                j = i + 1
                while j < len(cells):
                    found = False
                    for next_cell in cells[j]:
                        if next_cell.strip():
                            raw = next_cell.replace(";", "").strip()
                            new_peptides = [pep.strip() for pep in raw.split(" - ") if pep.strip() and len(pep.strip()) > 1]
                            peptides.extend(new_peptides)
                            found = True
                            break
                    if found:
                        break
                    j += 1
        i += 1
    return peptides

def main():
    species = input("Enter the species name (species_name): ").strip()
    
    input_folder = "/home/marta/Desktop/biopep_results/batch_processing_results/odsy/"
    output_folder = "/home/marta/Desktop/biopep_results/peptide_lists/all_peptides"
    os.makedirs(output_folder, exist_ok=True)
    
    all_peptides = []
    
    # Select only .ods files that contain the species name
    files = [file for file in os.listdir(input_folder)
             if file.endswith(".ods") and species in file]
    
    for file in files:
        filepath = os.path.join(input_folder, file)
        peptides = extract_peptides_from_file(filepath)
        if peptides:
            all_peptides.extend(peptides)
    
    outpath = os.path.join(output_folder, f"{species}.txt")
    with open(outpath, "w") as f:
        f.write("\n".join(all_peptides))
    
    print(f"Saved {len(all_peptides)} peptides to: {outpath}")

if __name__ == "__main__":
    main()
