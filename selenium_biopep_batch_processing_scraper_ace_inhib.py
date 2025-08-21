#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BIOPEP UWM Selenium Automation Script
Automatyzuje proces analizy peptydów na stronie BIOPEP UWM
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

DEFAULT_TIMEOUT = 1000  # w sekundach

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
            print("🔧 Inicjalizacja drivera przeglądarki...")
            options = webdriver.ChromeOptions()
            options.add_argument('--no-sandbox')
            options.add_argument('--disable-dev-shm-usage')
            options.add_argument('--disable-gpu')
            options.add_argument('--user-agent=Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36')
            if self.headless:
                options.add_argument('--headless=new')
                print("🔧 Uruchamianie w trybie headless")
            else:
                print("🔧 Uruchamianie z interfejsem graficznym")
            service = Service(ChromeDriverManager().install())
            self.driver = webdriver.Chrome(service=service, options=options)

            # 👉 Dodane ustawienie globalnego timeoutu dla komend WebDrivera
            self.driver.command_executor.set_timeout(self.timeout)

            self.driver.set_page_load_timeout(self.timeout + 60)
            self.wait = WebDriverWait(self.driver, self.timeout)
            print("✓ Driver przeglądarki został pomyślnie zainicjalizowany")
        except WebDriverException as e:
            print(f"❌ Błąd podczas inicjalizacji przeglądarki: {e}")
            raise

    def get_file_path(self):
        while True:
            filename = input("Podaj nazwę pliku FASTA (bez rozszerzenia): ").strip()
            if not filename:
                print("❌ Nazwa pliku nie może być pusta!")
                continue
            full_path = os.path.join(self.base_folder, f"{filename}.fasta")
            if os.path.exists(full_path):
                print(f"✓ Znaleziono plik: {full_path}")
                return full_path
            else:
                print(f"❌ Plik {full_path} nie istnieje!")
                retry = input("Czy chcesz spróbować ponownie? (t/n): ").lower()
                if retry != 't':
                    return None

    def navigate_to_biopep(self):
        url = "https://biochemia.uwm.edu.pl/biopep/batch_processing.php?but_processing.x=61&but_processing.y=16"
        try:
            print(f"🌐 Nawigacja do: {url}")
            self.driver.get(url)
            self.wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
            print("✓ Strona BIOPEP została załadowana")
        except TimeoutException:
            print("❌ Timeout podczas ładowania strony BIOPEP")
            raise

    def read_fasta_file(self, file_path):
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            print(f"✓ Przeczytano plik FASTA: {os.path.basename(file_path)} ({len(content)} znaków)")
            return content
        except Exception as e:
            print(f"❌ Błąd podczas czytania pliku FASTA: {e}")
            raise

    def configure_analysis_settings(self, file_path):
        try:
            print("⚙️ Konfigurowanie ustawień analizy...")
            fasta_content = self.read_fasta_file(file_path)

            # Odznacz niepotrzebne checkboxy
            checkboxes_to_uncheck = ["profiles", "profilemarker", "calculationsatas", "profileepi"]
            for checkbox_value in checkboxes_to_uncheck:
                try:
                    checkbox = self.driver.find_element(By.CSS_SELECTOR, f"input[name='wybor[]'][value='{checkbox_value}']")
                    if checkbox.is_selected():
                        checkbox.click()
                        print(f"✓ Odkliknięto: {checkbox_value}")
                except NoSuchElementException:
                    pass

            # Zaznacz wymagane checkboxy
            required_checkboxes = ["calculationsab", "enzymesaction", "searchfragments", "calculationsdvb"]
            for option in required_checkboxes:
                try:
                    checkbox = self.driver.find_element(By.CSS_SELECTOR, f"input[name='wybor[]'][value='{option}']")
                    if not checkbox.is_selected():
                        checkbox.click()
                    print(f"✓ Zaznaczono opcję: {option}")
                except NoSuchElementException:
                    pass

            # Aktywność
            try:
                activity_select = Select(self.driver.find_element(By.NAME, "sel_activity"))
                activity_select.select_by_value("ACE inhibitor")
                print("✓ Wybrano aktywność: ACE inhibitor")
            except NoSuchElementException:
                pass

            # Wprowadzanie sekwencji
            try:
                textarea = self.wait.until(EC.presence_of_element_located((By.NAME, "txt_seq")))
                self.driver.execute_script("arguments[0].value = arguments[1];", textarea, fasta_content)
                print(f"✓ Wprowadzono sekwencje do textarea ({len(fasta_content)} znaków)")
            except TimeoutException:
                print("❌ Nie znaleziono pola textarea dla sekwencji")
                raise

            # Enzymy
            enzyme_settings = {"enz1": "13", "enz2": "12", "enz3": "11"}
            enzyme_names = {"enz1": "pepsin (pH 1.3)", "enz2": "trypsin", "enz3": "chymotrypsin (A)"}
            for enz_name, enz_value in enzyme_settings.items():
                try:
                    enzyme_select = Select(self.driver.find_element(By.NAME, enz_name))
                    enzyme_select.select_by_value(enz_value)
                    print(f"✓ Wybrano {enz_name}: {enzyme_names[enz_name]} (wartość: {enz_value})")
                except NoSuchElementException:
                    pass

            # Baza danych
            try:
                database_select = self.driver.find_element(By.NAME, "sel_database")
                database_sel = Select(database_select)
                database_sel.select_by_value("pep")
                print("✓ Wybrano bazę danych: pep")
            except NoSuchElementException:
                pass

            print("✅ Konfiguracja ustawień zakończona pomyślnie")
        except Exception as e:
            print(f"❌ Błąd podczas konfigurowania ustawień: {e}")
            raise

    def start_analysis(self):
        try:
            print("🚀 Uruchamianie analizy...")
            report_button = self.wait.until(
                EC.element_to_be_clickable((By.NAME, "but_report"))
            )
            report_button.click()
            print("✓ Kliknięto przycisk uruchamiający analizę")
            print("⏳ Oczekiwanie na wyniki analizy...")

            # Czekamy aż pojawi się formularz wynikowy
            self.wait.until(
                EC.presence_of_element_located(
                    (By.CSS_SELECTOR, 'form[name="fm_action"][action="batch_report_cutting_for_seq"] table.table-out')
                )
            )

            print("✓ Formularz wynikowy został załadowany")

        except TimeoutException:
            print("❌ Timeout podczas oczekiwania na wyniki")
            raise

    def extract_sequences(self):
        try:
            print("🔍 Analizowanie wyników...")
            tables = self.driver.find_elements(By.CSS_SELECTOR, "form[name='fm_action'] table.table-out")
            if not tables:
                print("❌ Nie znaleziono tabeli z wynikami")
                return
            print(f"✓ Znaleziono {len(tables)} tabel do analizy")

            filtered_tables = []
            for table in tables:
                rows = table.find_elements(By.CSS_SELECTOR, "tr")
                if not rows:
                    continue

                # znajdź nagłówek
                header_row = None
                header_idx = 0
                for idx, tr in enumerate(rows):
                    if tr.find_elements(By.TAG_NAME, "th"):
                        header_row = tr
                        header_idx = idx
                        break
                if header_row is None:
                    header_row = rows[0]
                    header_idx = 0

                header_cells = header_row.find_elements(By.TAG_NAME, "th")
                if not header_cells:
                    header_cells = header_row.find_elements(By.TAG_NAME, "td")

                headers = [" ".join(c.text.split()).strip().lower() for c in header_cells]
                if any("sequence" in h for h in headers):
                    filtered_tables.append((table, rows, headers, header_idx))

            print(f"✓ Do dalszej analizy wybrano {len(filtered_tables)} tabel z kolumną 'Sequence'")

            for table, rows, headers, header_idx in filtered_tables:
                seq_idx = None
                act_idx = None
                for i, h in enumerate(headers):
                    if seq_idx is None and "sequence" in h:
                        seq_idx = i
                    if act_idx is None and "activity" in h:
                        act_idx = i

                if seq_idx is None or act_idx is None:
                    print("⚠️ Pomijam tabelę: brak wymaganych kolumn 'Sequence' i/lub 'Activity'")
                    continue

                data_rows = rows[header_idx + 1:]
                for row in data_rows:
                    cells = row.find_elements(By.TAG_NAME, "td")
                    if not cells or len(cells) <= max(seq_idx, act_idx):
                        continue

                    activity = cells[act_idx].text.strip()
                    if "ace inhibitor" in activity.lower():
                        sequence = cells[seq_idx].text.strip()
                        if sequence and sequence not in self.results:
                            self.results.append(sequence)
                            print(f"✓ Znaleziono sekwencję: {sequence}")

            if not self.results:
                print("⚠️ Nie znaleziono sekwencji z aktywnością 'ACE inhibitor'")

        except Exception as e:
            print(f"❌ Błąd podczas wyodrębniania sekwencji: {e}")
            raise

    def save_results(self, species_name):
        try:
            if not self.results:
                print("⚠️ Brak wyników do zapisania")
                return
            timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
            # Dodanie nazwy gatunku do nazwy pliku
            file_path = os.path.join(self.base_folder, f"biopep_results_{species_name}_{timestamp}.txt")
            with open(file_path, 'w', encoding='utf-8') as f:
                for seq in self.results:
                    f.write(seq + "\n")
            print(f"💾 Wyniki zapisano w pliku: {file_path}")
        except Exception as e:
            print(f"❌ Błąd podczas zapisywania wyników: {e}")

    def close_driver(self):
        if self.driver:
            self.driver.quit()
            print("🚪 Zamknięto przeglądarkę")


import sys

def main():
    # ✅ Zapytaj użytkownika czy chce tryb headless
    headless_input = input("Czy uruchomić w trybie headless? (t/n): ").strip().lower()
    headless_mode = headless_input == 't'

    analyzer = BIOPEPAnalyzer(headless=headless_mode)
    analyzer.setup_driver()
    try:
        file_path = analyzer.get_file_path()
        if not file_path:
            print("❌ Nie podano pliku. Kończenie programu.")
            return

        # Wyciągnięcie nazwy gatunku z nazwy pliku (bez folderu i rozszerzenia)
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

