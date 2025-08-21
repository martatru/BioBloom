"""
BIOPEP UWM Selenium Automation Script
Automates peptide analysis on the BIOPEP UWM website
"""

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException, WebDriverException
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
import os
from datetime import datetime

DEFAULT_TIMEOUT = 1000

class BIOPEPAnalyzer:
    def __init__(self, headless=False, timeout=DEFAULT_TIMEOUT):
        self.driver = None
        self.wait = None
        self.results = []
        self.headless = headless
        self.timeout = timeout
        self.base_folder = "/home/marta/Desktop/biobloom_proteomic_data/accurate_proteomes"

    def setup_driver(self):
        try:
            print("Initializing browser driver...")
            options = webdriver.ChromeOptions()
            options.add_argument('--no-sandbox')
            options.add_argument('--disable-dev-shm-usage')
            options.add_argument('--disable-gpu')
            options.add_argument('--user-agent=Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36')
            if self.headless:
                options.add_argument('--headless=new')
                print("Running in headless mode")
            else:
                print("Running with graphical interface")
            service = Service(ChromeDriverManager().install())
            self.driver = webdriver.Chrome(service=service, options=options)
            self.driver.command_executor.set_timeout(self.timeout)
            self.driver.set_page_load_timeout(self.timeout + 60)
            self.wait = WebDriverWait(self.driver, self.timeout)
            print("Browser driver initialized successfully")
        except WebDriverException as e:
            print(f"Error initializing browser: {e}")
            raise

    def get_file_path(self):
        while True:
            filename = input("Enter FASTA file name (without extension): ").strip()
            if not filename:
                print("File name cannot be empty!")
                continue
            full_path = os.path.join(self.base_folder, f"{filename}.fasta")
            if os.path.exists(full_path):
                print(f"File found: {full_path}")
                return full_path
            else:
                print(f"File {full_path} does not exist!")
                retry = input("Try again? (y/n): ").lower()
                if retry != 'y':
                    return None

    def navigate_to_biopep(self):
        url = "https://biochemia.uwm.edu.pl/biopep/batch_processing.php?but_processing.x=61&but_processing.y=16"
        try:
            print(f"Navigating to: {url}")
            self.driver.get(url)
            self.wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
            print("BIOPEP page loaded")
        except TimeoutException:
            print("Timeout while loading BIOPEP page")
            raise

    def read_fasta_file(self, file_path):
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            print(f"Read FASTA file: {os.path.basename(file_path)} ({len(content)} characters)")
            return content
        except Exception as e:
            print(f"Error reading FASTA file: {e}")
            raise

    def configure_analysis_settings(self, file_path):
        try:
            print("Configuring analysis settings...")
            fasta_content = self.read_fasta_file(file_path)

            # Uncheck unnecessary checkboxes
            checkboxes_to_uncheck = ["profiles", "profilemarker", "calculationsatas", "profileepi"]
            for checkbox_value in checkboxes_to_uncheck:
                try:
                    checkbox = self.driver.find_element(By.CSS_SELECTOR, f"input[name='wybor[]'][value='{checkbox_value}']")
                    if checkbox.is_selected():
                        checkbox.click()
                        print(f"Unchecked: {checkbox_value}")
                except NoSuchElementException:
                    pass

            # Check required checkboxes
            required_checkboxes = ["calculationsab", "enzymesaction", "searchfragments", "calculationsdvb"]
            for option in required_checkboxes:
                try:
                    checkbox = self.driver.find_element(By.CSS_SELECTOR, f"input[name='wybor[]'][value='{option}']")
                    if not checkbox.is_selected():
                        checkbox.click()
                    print(f"Checked: {option}")
                except NoSuchElementException:
                    pass

            # Activity selection
            try:
                activity_select = Select(self.driver.find_element(By.NAME, "sel_activity"))
                activity_select.select_by_value("ACE inhibitor")
                print("Selected activity: ACE inhibitor")
            except NoSuchElementException:
                pass

            # Input sequence
            try:
                textarea = self.wait.until(EC.presence_of_element_located((By.NAME, "txt_seq")))
                self.driver.execute_script("arguments[0].value = arguments[1];", textarea, fasta_content)
                print(f"Entered sequence ({len(fasta_content)} characters)")
            except TimeoutException:
                print("Textarea for sequence not found")
                raise

            # Enzyme selection
            enzyme_settings = {"enz1": "13", "enz2": "12", "enz3": "11"}
            enzyme_names = {"enz1": "pepsin (pH 1.3)", "enz2": "trypsin", "enz3": "chymotrypsin (A)"}
            for enz_name, enz_value in enzyme_settings.items():
                try:
                    enzyme_select = Select(self.driver.find_element(By.NAME, enz_name))
                    enzyme_select.select_by_value(enz_value)
                    print(f"Selected {enz_name}: {enzyme_names[enz_name]} (value: {enz_value})")
                except NoSuchElementException:
                    pass

            # Database selection
            try:
                database_select = self.driver.find_element(By.NAME, "sel_database")
                Select(database_select).select_by_value("pep")
                print("Selected database: pep")
            except NoSuchElementException:
                pass

            print("Settings configured successfully")
        except Exception as e:
            print(f"Error configuring settings: {e}")
            raise

    def start_analysis(self):
        try:
            print("Starting analysis...")
            report_button = self.wait.until(EC.element_to_be_clickable((By.NAME, "but_report")))
            report_button.click()
            print("Clicked button to start analysis")
            self.wait.until(
                EC.presence_of_element_located(
                    (By.CSS_SELECTOR, 'form[name="fm_action"][action="batch_report_cutting_for_seq"] table.table-out')
                )
            )
            print("Result form loaded")
        except TimeoutException:
            print("Timeout waiting for results")
            raise

    def extract_sequences(self):
        try:
            print("Extracting enzyme action results...")
            tables = self.driver.find_elements(By.CSS_SELECTOR, "form[name='fm_action'] table.table-out")
            if not tables:
                print("No result tables found")
                return

            filtered_tables = []
            for table in tables:
                header_cells = table.find_elements(By.TAG_NAME, "th")
                headers = [c.text.strip().lower() for c in header_cells]
                if any("results of enzyme action" in h for h in headers):
                    filtered_tables.append(table)

            print(f"{len(filtered_tables)} tables selected with 'Results of enzyme action' header")

            for table in filtered_tables:
                td_elements = table.find_elements(By.CSS_SELECTOR, "td.info font[size='-1']")
                for td in td_elements:
                    text = td.text.strip()
                    peptides = [p.strip() for p in text.split('-') if len(p.strip()) > 1]
                    for peptide in peptides:
                        if peptide not in self.results:
                            self.results.append(peptide)
                            print(f"Found peptide: {peptide}")

            if not self.results:
                print("No peptides found in enzyme action tables")

        except Exception as e:
            print(f"Error extracting sequences: {e}")
            raise

    def save_results(self, species_name):
        try:
            if not self.results:
                print("No results to save")
                return
            # Save output in the requested folder
            output_folder = os.path.join(
                self.base_folder, 
                "peptide_lists", 
                "biopep_enzymes_action_results"
            )
            os.makedirs(output_folder, exist_ok=True)

            timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
            file_path = os.path.join(output_folder, f"biopep_peptides_{species_name}_{timestamp}.txt")

            with open(file_path, 'w', encoding='utf-8') as f:
                for peptide in self.results:
                    f.write(peptide + "\n")

            print(f"Peptide results saved to: {file_path}")

        except Exception as e:
            print(f"Error saving results: {e}")

    def close_driver(self):
        if self.driver:
            self.driver.quit()
            print("Browser closed")


def main():
    headless_input = input("Run in headless mode? (y/n): ").strip().lower()
    headless_mode = headless_input == 'y'

    analyzer = BIOPEPAnalyzer(headless=headless_mode)
    analyzer.setup_driver()
    try:
        file_path = analyzer.get_file_path()
        if not file_path:
            print("No file provided. Exiting program.")
            return

        species_name = os.path.splitext(os.path.basename(file_path))[0]

        analyzer.navigate_to_biopep()
        analyzer.configure_analysis_settings(file_path)
        analyzer.start_analysis()
        analyzer.extract_sequences()
        analyzer.save_results(species_name)
    finally:
        analyzer.close_driver()


if __name__ == "__main__":
    main()
