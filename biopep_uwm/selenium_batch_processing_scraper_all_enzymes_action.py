#!/usr/bin/env python3
"""
BIOPEP UWM — Selenium batch runner (Enzymes' action → peptide list)

Purpose
-------
Automates BIOPEP UWM "Batch processing" to collect peptides from the
**Results of enzyme action** tables.

This version is fully CLI‑driven (no interactive inputs), uses portable
paths, robust logging, and a context‑managed browser session.

Quick start
-----------
    python biopep_uwm_enzyme_action_refactored.py \
        --fasta data/proteome.fasta \
        --out-dir results/enzymes_action \
        --headless \
        --timeout 600

Outputs
-------
- A newline‑separated TXT with unique peptide sequences extracted from the
  "Results of enzyme action" sections. Filename includes FASTA stem + timestamp.

Notes
-----
- CSS selectors are defensive; if the site markup changes, adjust constants
  in `extract()` accordingly.
- Requires Chrome + webdriver‑manager.
"""
from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException, WebDriverException
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DEFAULT_TIMEOUT = 600
BIOPEP_BATCH_URL = (
    "https://biochemia.uwm.edu.pl/biopep/batch_processing.php?but_processing.x=61&but_processing.y=16"
)
LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("biopep.enzyme_action")


@dataclass
class RunConfig:
    fasta_path: Path
    out_dir: Path
    headless: bool = False
    timeout: int = DEFAULT_TIMEOUT
    url: str = BIOPEP_BATCH_URL

    def validate(self) -> None:
        if not self.fasta_path.exists() or not self.fasta_path.is_file():
            raise FileNotFoundError(f"FASTA not found: {self.fasta_path}")
        self.out_dir.mkdir(parents=True, exist_ok=True)


class BIOPEPEnzymeAction:
    """Encapsulates browser session and scraping logic for enzyme action results."""

    def __init__(self, cfg: RunConfig):
        self.cfg = cfg
        self.driver: Optional[webdriver.Chrome] = None
        self.wait: Optional[WebDriverWait] = None
        self.results: List[str] = []

    # --- lifecycle ---------------------------------------------------------
    def __enter__(self):
        self.setup_driver()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close_driver()
        return False

    def setup_driver(self) -> None:
        logger.info("Initializing Chrome (headless=%s)", self.cfg.headless)
        try:
            options = webdriver.ChromeOptions()
            options.add_argument("--no-sandbox")
            options.add_argument("--disable-dev-shm-usage")
            options.add_argument("--disable-gpu")
            options.add_argument(
                "--user-agent=Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
            )
            if self.cfg.headless:
                options.add_argument("--headless=new")

            service = Service(ChromeDriverManager().install())
            self.driver = webdriver.Chrome(service=service, options=options)
            self.driver.command_executor.set_timeout(self.cfg.timeout)
            self.driver.set_page_load_timeout(self.cfg.timeout + 60)
            self.wait = WebDriverWait(self.driver, self.cfg.timeout)
            logger.info("Chrome ready")
        except WebDriverException as e:
            logger.error("Failed to initialize Chrome: %s", e)
            raise

    # --- flow --------------------------------------------------------------
    def navigate(self) -> None:
        assert self.driver and self.wait
        logger.info("Opening page: %s", self.cfg.url)
        try:
            self.driver.get(self.cfg.url)
            self.wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
            logger.info("BIOPEP page loaded")
        except TimeoutException:
            logger.error("Timeout while loading BIOPEP page")
            raise

    def _read_fasta(self) -> str:
        text = self.cfg.fasta_path.read_text(encoding="utf-8")
        logger.info("Read FASTA %s (%,d chars)", self.cfg.fasta_path.name, len(text))
        return text

    def configure(self) -> None:
        assert self.driver and self.wait
        logger.info("Configuring analysis settings…")
        fasta_content = self._read_fasta()

        # Uncheck non‑needed options if present
        for value in ["profiles", "profilemarker", "calculationsatas", "profileepi"]:
            try:
                el = self.driver.find_element(By.CSS_SELECTOR, f"input[name='wybor[]'][value='{value}']")
                if el.is_selected():
                    el.click()
                    logger.debug("Unchecked %s", value)
            except NoSuchElementException:
                pass

        # Required options: include enzyme action & fragments (as in original)
        for value in ["calculationsab", "enzymesaction", "searchfragments", "calculationsdvb"]:
            try:
                el = self.driver.find_element(By.CSS_SELECTOR, f"input[name='wybor[]'][value='{value}']")
                if not el.is_selected():
                    el.click()
                logger.debug("Checked %s", value)
            except NoSuchElementException:
                pass

        # Activity setting retained (some pages use it as a filter); non‑fatal if absent
        try:
            Select(self.driver.find_element(By.NAME, "sel_activity")).select_by_value("ACE inhibitor")
            logger.debug("Selected activity: ACE inhibitor")
        except NoSuchElementException:
            logger.warning("Activity select not found; continuing")

        # Paste FASTA
        try:
            textarea = self.wait.until(EC.presence_of_element_located((By.NAME, "txt_seq")))
            self.driver.execute_script("arguments[0].value = arguments[1];", textarea, fasta_content)
            logger.info("Pasted sequence (%d chars)", len(fasta_content))
        except TimeoutException:
            logger.error("Sequence textarea not found")
            raise

        # Enzymes: pepsin (pH 1.3), trypsin, chymotrypsin (A)
        for name, value in {"enz1": "13", "enz2": "12", "enz3": "11"}.items():
            try:
                Select(self.driver.find_element(By.NAME, name)).select_by_value(value)
                logger.debug("Selected %s=%s", name, value)
            except NoSuchElementException:
                pass

        # Database: pep
        try:
            Select(self.driver.find_element(By.NAME, "sel_database")).select_by_value("pep")
            logger.debug("Selected database: pep")
        except NoSuchElementException:
            pass

        logger.info("Settings configured")

    def run(self) -> None:
        assert self.wait
        logger.info("Starting analysis…")
        try:
            btn = self.wait.until(EC.element_to_be_clickable((By.NAME, "but_report")))
            btn.click()
            self.wait.until(
                EC.presence_of_element_located(
                    (By.CSS_SELECTOR, "form[name='fm_action'][action='batch_report_cutting_for_seq'] table.table-out")
                )
            )
            logger.info("Result tables present")
        except TimeoutException:
            logger.error("Timeout waiting for results")
            raise

    def extract(self) -> List[str]:
        assert self.driver
        logger.info("Extracting peptides from enzyme action tables…")
        tables = self.driver.find_elements(By.CSS_SELECTOR, "form[name='fm_action'] table.table-out")
        if not tables:
            logger.warning("No result tables found")
            return []

        selected_tables = []
        for table in tables:
            try:
                header_td = table.find_element(By.CSS_SELECTOR, "td.infobold font[size='-1']")
                if "results of enzyme action" in header_td.text.strip().lower():
                    selected_tables.append(table)
            except NoSuchElementException:
                continue

        logger.info("Selected %d enzyme action table(s)", len(selected_tables))

        peptides: set[str] = set()
        for table in selected_tables:
            for td in table.find_elements(By.CSS_SELECTOR, "td.info font[size='-1']"):
                text = td.text.strip()
                # Heuristic: peptides separated by dashes e.g. ALA-GLY-LYS-...
                for p in (frag.strip() for frag in text.split('-')):
                    if len(p) > 1:  # ignore single letters and empties
                        peptides.add(p)

        self.results = sorted(peptides)
        logger.info("Collected %d unique peptide entries", len(self.results))
        return self.results

    def save(self) -> Path:
        species = self.cfg.fasta_path.stem
        ts = datetime.now().strftime("%Y%m%d%H%M%S")
        out_path = self.cfg.out_dir / f"biopep_peptides_{species}_{ts}.txt"
        out_path.write_text("\n".join(self.results) + ("\n" if self.results else ""), encoding="utf-8")
        logger.info("Saved peptides to %s", out_path)
        return out_path

    def close_driver(self) -> None:
        if self.driver:
            self.driver.quit()
            logger.debug("Browser closed")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run BIOPEP UWM batch analysis and extract peptides from enzyme action results.")
    p.add_argument("--fasta", type=Path, required=True, help="Path to input FASTA file")
    p.add_argument("--out-dir", type=Path, default=Path("results/enzymes_action"), help="Directory to write outputs")
    p.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT, help="Timeout (s) for page loads and waits")
    p.add_argument("--headless", action="store_true", help="Run Chrome in headless mode")
    p.add_argument("--verbose", action="store_true", help="Enable debug logs")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    import sys

    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    cfg = RunConfig(
        fasta_path=args.fasta.expanduser().resolve(),
        out_dir=args.out_dir.expanduser().resolve(),
        headless=args.headless,
        timeout=int(args.timeout),
    )
    try:
        cfg.validate()
    except Exception as e:
        logger.error("Invalid configuration: %s", e)
        return 2

    try:
        with BIOPEPEnzymeAction(cfg) as runner:
            runner.navigate()
            runner.configure()
            runner.run()
            runner.extract()
            out_file = runner.save()
        logger.info("Done. Output: %s", out_file)
        return 0
    except Exception as e:
        logger.exception("Run failed: %s", e)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
