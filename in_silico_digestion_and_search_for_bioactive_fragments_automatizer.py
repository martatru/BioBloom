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
from tqdm import tqdm

#potrzebuje debugowania hihi 

# === ARGPARSE ===
parser = argparse.ArgumentParser(description="Download Biopep UWM digestion + fragments results")
parser.add_argument(
    "species",
    choices=["arthrospira_platensis", "microchloropsis_salina", "tetraselmis_suecica"],
    help="Species name (choose: arthrospira_platensis, microchloropsis_salina, tetraselmis_suecica)"
)
args = parser.parse_args()

# === CONFIGURATION ===
fasta_file = f"/home/marta/Desktop/biobloom_proteomic_data/{args.species}.fasta"
base_dir = Path.home() / "Desktop" / "biopep_results"

stage2_dir = base_dir / "stage2_in_sillico_digestion_results" / args.species
stage3_dir = base_dir / "stage3_seach_for_active_fragments_results" / args.species

# Download folder for Excel
excel_download_dir = base_dir / "downloads"
excel_download_dir.mkdir(parents=True, exist_ok=True)

# Chrome driver configuration
chrome_options = Options()
prefs = {
    "download.default_directory": str(excel_download_dir),
    "download.prompt_for_download": False,
    "download.directory_upgrade": True,
    "safebrowsing.enabled": True
}
chrome_options.add_experimental_option("prefs", prefs)
driver = webdriver.Chrome(options=chrome_options)

# === FUNCTION TO PROCESS A SINGLE PROTEIN ===
def process_protein(species, seq_id, sequence):
    # Stage2: digestion results directory
    digestion_file = stage2_dir / f"{seq_id}.txt"
    digestion_file.parent.mkdir(parents=True, exist_ok=True)

    # Stage3: fragments results directory
    fragments_dir = stage3_dir / seq_id
    fragments_dir.mkdir(parents=True, exist_ok=True)

    # Go to enzyme digestion form
    driver.get("https://biochemia.uwm.edu.pl/biopep/enzymes_act_for_your_seq.php?txt_seq_e=&txt_seq=&enz1=&but_vprotein2.x=33&but_vprotein2.y=11&enz2=&enz3=&prot=&e2=&e3=&ktory=1&e1=eee")

    # Insert protein sequence
    try:
        textarea = WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.NAME, "txt_seq"))
        )
        textarea.clear()
        textarea.send_keys(sequence)
    except:
        print(f"[ERROR] Textarea not found for {seq_id}")
        return

    # Click report button
    try:
        report_btn = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.NAME, "but_report"))
        )
        report_btn.click()
    except:
        print(f"[ERROR] Could not click report button for {seq_id}")
        return

    # Wait for results textarea and save content
    try:
        results_area = WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.NAME, "report"))
        )
        with open(digestion_file, "w") as f:
            f.write(results_area.get_attribute("value"))
        print(f"[OK] Stage2 results saved: {digestion_file}")
    except:
        print(f"[ERROR] Could not save digestion results for {seq_id}")
        return

    # Click "Search for active fragments"
    try:
        frag_btn = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, "//input[@src='images/search_for_act_frag0_normal.gif']"))
        )
        frag_btn.click()
    except:
        print(f"[ERROR] Could not click 'Search for active fragments' for {seq_id}")
        return

    # Click export to excel
    try:
        export_btn = WebDriverWait(driver, 20).until(
            EC.element_to_be_clickable((By.XPATH, "//input[@src='images/export.gif']"))
        )
        driver.execute_script("arguments[0].click();", export_btn)
    except:
        print(f"[ERROR] Could not click export for {seq_id}")
        return

    # Wait for download
    downloaded_file = None
    for _ in range(30):
        files = list(excel_download_dir.glob("*.xls*"))
        if files:
            downloaded_file = max(files, key=lambda f: f.stat().st_mtime)
            break
        time.sleep(1)

    if downloaded_file:
        target = fragments_dir / downloaded_file.name
        try:
            downloaded_file.rename(target)
            print(f"[OK] Stage3 results saved: {target}")
        except Exception as e:
            print(f"[ERROR] Could not move excel file for {seq_id}: {e}")
    else:
        print(f"[ERROR] No excel file downloaded for {seq_id}")

# === MAIN ===
records = list(SeqIO.parse(fasta_file, "fasta"))
total = len(records)
print(f"Processing {total} proteins for {args.species}...\n")

for record in tqdm(records, desc="Proteins processed"):
    process_protein(args.species, record.id, str(record.seq))

driver.quit()
