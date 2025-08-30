# BioBloom Documentation 📚

Full documentation for the **BioBloom** project repository.  

---

## Project Overview

BioBloom explores **microalgae proteomes** as a source of bioactive compounds with pharmacological potential.  
We focus on **bioactive peptides** (esp. ACE inhibitors for cardiovascular health) by:

1. Mining microalgal proteomes,  
2. Querying the **BIOPEP-UWM** database,  
3. Filtering peptide lists & removing duplicates,  
4. Screening ADMET properties with AdmetLab 3.0,  
5. Converting between FASTA ⇄ SMILES ⇄ PDB,  
6. Generating peptide 3D structures,  
7. Preparing receptor structures,  
8. Producing input/output for pLM4ACE and docking pipelines.  

---

## Repository Structure

### `biopep_uwm/`
Automation of the [BIOPEP-UWM](https://biochemia.uwm.edu.pl/biopep/start_biopep.php) batch processing:
- `selenium_biopep_batch_processing.py` – ACE inhibitory activity screening.  
- `selenium_batch_processing_scraper.py` – enzyme action analysis.  
- `search_for_novel_peptides.py` – compare BIOPEP results with known ACE inhibitors.  
- `unify_a_platensis_biopep_output.py` – merge & deduplicate *A. platensis* BIOPEP peptide lists.  

### `admet/`
Filtering and unifying ADMET screening results:
- `unify_admet_output_a_platensis.py` – merge species-specific ADMET outputs.  

### `molecular_docking/`
Peptide/receptor structure preparation:
- `generate_peptide_structures_pyrosetta.py` – build PDBs from FASTA with PyRosetta.  
- `repack_receptor_pyrosetta.py` – side-chain repacking of ACE receptor.  
- `place_pep_into_ace.py` – peptide placement into ACE binding pocket.  
- `select_top_peptides_for_molecular_docking.py` – rank/filter peptides by ADMET criteria.  

### `pLM4ACE/`
Formatting utilities for the pLM4ACE model:
- `prepare_pLM4ACE_input.py` – create input files.  
- `split_pLM4ACE_input.py` – split input into batches.  
- `unify_pLM4ACE_output.py` – merge model outputs.  

### `smiles_conversion/`
Cheminformatics utilities:
- `smiles_converter.py` – bidirectional FASTA ⇄ SMILES.  
- `create_fasta_input_for_smiles_conversion.py` – prep FASTA for conversion.  
- `extract_smiles_without_names.py` – clean SMILES files (remove labels, keep strings only).  

---

## Requirements

- **Python 3.9+**  
- **PyRosetta** (licensed, manual install)  
- Chrome + Chromedriver (auto-installed with `webdriver-manager`)  

### Python dependencies
```bash
pip install biopython selenium tqdm openpyxl pandas numpy rdkit-pypi
