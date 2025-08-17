import os
import time
import argparse
from pathlib import Path
from Bio import SeqIO
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from tqdm import tqdm  # progress bar

# === ARGPARSE ===
parser = argparse.ArgumentParser(description="Download Biopep UWM results for a given species")
parser.add_argument(
    "species",
    choices=["arthrospira_platensis", "microchloropsis_salina", "tetraselmis_suecica"],
    help="Species name (choose: arthrospira_platensis, microchloropsis_salina, tetraselmis_suecica)"
)
args = parser.parse_args()

# === CONFIGURATION ===
fasta_file = f"/home/marta/Desktop/biobloom_proteomic_data/{args.species}.fasta"
base_dir = Path.home() / "Desktop" / "biopep_results"

# Download folder
download_dir_path = base_dir / "downloads"
download_dir_path.mkdir(parents=True, exist_ok=True)
download_dir = str(download_dir_path)

# Chrome settings for automatic Excel download
chrome_options = Options()
prefs = {
    "download.default_directory": download_dir,
    "download.prompt_for_download": False,
    "download.directory_upgrade": True,
    "safebrowsing.enabled": True
}
chrome_options.add_experimental_option("prefs", prefs)
driver = webdriver.Chrome(options=chrome_options)

# === FUNCTION TO PROCESS A SINGLE PROTEIN ===
def process_protein(species, seq_id, sequence):
    protein_dir = base_dir / species / seq_id
    protein_dir.mkdir(parents=True, exist_ok=True)

    driver.get("https://biochemia.uwm.edu.pl/biopep/finding_aby.php?sel_activity=ACE+inhibitor")

    # Wait for textarea and input sequence
    try:
        textarea = WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.NAME, "txt_seq"))
        )
        textarea.clear()
        textarea.send_keys(sequence)
    except:
        print(f"[ERROR] Textarea not found for {seq_id} ({species})")
        return

    # Click "Report" button
    try:
        report_button = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.NAME, "cmdOK2"))
        )
        report_button.click()
    except:
        print(f"[ERROR] Report button not clickable for {seq_id} ({species})")
        return

    # Wait for "Export to EXCEL" button and click via JavaScript
    try:
        export_btn = WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.XPATH, "//input[@value='export to excel']"))
        )
        driver.execute_script("arguments[0].click();", export_btn)
    except:
        print(f"[ERROR] Could not retrieve results for {seq_id} ({species})")
        return

    # Wait for file download
    downloaded_file = None
    for _ in range(30):  # wait up to 30 seconds
        files = list(download_dir_path.glob("*.xls*"))
        if files:
            downloaded_file = max(files, key=lambda f: f.stat().st_mtime)
            break
        time.sleep(1)

    if downloaded_file:
        target = protein_dir / downloaded_file.name
        try:
            downloaded_file.rename(target)
            print(f"[OK] Saved: {target}")
        except Exception as e:
            print(f"[ERROR] Could not move file for {seq_id}: {e}")
    else:
        print(f"[ERROR] No file downloaded for {seq_id} ({species})")

# === MAIN LOOP WITH PROGRESS BAR ===
records = list(SeqIO.parse(fasta_file, "fasta"))
total_proteins = len(records)

print(f"Processing {total_proteins} proteins for species '{args.species}'...\n")

for record in tqdm(records, desc="Proteins processed"):
    seq_id = record.id
    sequence = str(record.seq)
    species = args.species
    process_protein(species, seq_id, sequence)

driver.quit()
