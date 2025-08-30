#!/usr/bin/env python3
"""
BIOPEP UWM — Selenium batch runner (ACE inhibitor extraction)

Purpose
-------
Automates the "Batch processing" workflow on BIOPEP UWM and extracts peptide
sequences annotated with **ACE inhibitor** activity from the resulting tables.

Why this exists
---------------
The original script used interactive prompts and hard‑coded paths. This version
is fully CLI‑driven, portable, and easier to reuse in pipelines and notebooks.

Quick start
-----------
    python biopep_uwm_automation_refactored.py \
        --fasta data/proteome.fasta \
        --out-dir results/ \
        --headless \
        --timeout 600

Inputs/Outputs
--------------
- Input: a FASTA text file passed via --fasta
- Output: a newline‑separated .txt file with unique ACE‑inhibitory sequences
  saved under --out-dir (timestamped filename)

Notes
-----
- Uses webdriver-manager to fetch a matching ChromeDriver.
- Selectors are intentionally conservative and defensive against minor DOM
  changes; if the site updates markup, adjust the CSS/XPath constants below.
"""
from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import (
    TimeoutException,
    NoSuchElementException,
    WebDriverException,
)
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

# -----------------------------------------------------------------------------
# Configuration & constants
# -----------------------------------------------------------------------------

DEFAULT_TIMEOUT = 600  # seconds; tune via --timeout
BIOPEP_BATCH_URL = (
    "https://biochemia.uwm.edu.pl/biopep/batch_processing.php?but_processing.x=61&but_processing.y=16"
)

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("biopep")


@dataclass
class RunConfig:
    fasta_path: Path
    out_dir: Path
    headless: bool = False
    timeout: int = DEFAULT_TIMEOUT
    url: str = BIOPEP_BATCH_URL

    def validate(self) -> None:
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA not found: {self.fasta_path}")
        if not self.fasta_path.is_file():
            raise ValueError(f"Expected a file, got: {self.fasta_path}")
        self.out_dir.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------
# Core automation
# -----------------------------------------------------------------------------

class BIOPEPAnalyzer:
    """Encapsulates browser session and scraping logic for BIOPEP batch runs."""

    def __init__(self, cfg: RunConfig):
        self.cfg = cfg
        self.driver: Optional[webdriver.Chrome] = None
        self.wait: Optional[WebDriverWait] = None
        self.results: List[str] = []

    # --- lifecycle ---------------------------------------------------------
    def __enter__(self):  # context manager for safe teardown
        self.setup_driver()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close_driver()
        # propagate exception if any
        return False

    def setup_driver(self) -> None:
        """Start Chrome with sensible defaults; prefer headless for CI runs."""
        logger.info("Initializing Chrome driver (headless=%s)", self.cfg.headless)
        try:
            options = webdriver.ChromeOptions()
            options.add_argument("--no-sandbox")
            options.add_argument("--disable-dev-shm-usage")
            options.add_argument("--disable-gpu")
            options.add_argument(
                "--user-agent=Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
            )
            if self.cfg.headless:
                # newer Chrome supports headless=new for better rendering parity
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

    # --- high‑level steps --------------------------------------------------
    def navigate(self) -> None:
        logger.info("Opening BIOPEP batch page: %s", self.cfg.url)
        try:
            assert self.driver is not None and self.wait is not None
            self.driver.get(self.cfg.url)
            self.wait.until(EC.presence_of_element_located((By.TAG_NAME, "body")))
            logger.info("Page loaded")
        except TimeoutException:
            logger.error("Timeout while loading BIOPEP page")
            raise

    def read_fasta(self) -> str:
        text = self.cfg.fasta_path.read_text(encoding="utf-8")
        logger.info("Read FASTA: %s (%,d chars)", self.cfg.fasta_path.name, len(text))
        return text

    def configure(self) -> None:
        """Set checkboxes, select activity and enzymes, paste FASTA."""
        logger.info("Configuring analysis settings…")
        fasta_content = self.read_fasta()
        assert self.driver is not None and self.wait is not None

        # Uncheck non‑essential checkboxes if present
        for value in ["profiles", "profilemarker", "calculationsatas", "profileepi"]:
            try:
                el = self.driver.find_element(By.CSS_SELECTOR, f"input[name='wybor[]'][value='{value}']")
                if el.is_selected():
                    el.click()
                    logger.debug("Unchecked %s", value)
            except NoSuchElementException:
                # Site layout may change; silently continue
                pass

        # Ensure required options are checked
        for value in ["calculationsab", "enzymesaction", "searchfragments", "calculationsdvb"]:
            try:
                el = self.driver.find_element(By.CSS_SELECTOR, f"input[name='wybor[]'][value='{value}']")
                if not el.is_selected():
                    el.click()
                logger.debug("Checked %s", value)
            except NoSuchElementException:
                pass

        # Activity: ACE inhibitor
        try:
            Select(self.driver.find_element(By.NAME, "sel_activity")).select_by_value("ACE inhibitor")
            logger.debug("Selected activity: ACE inhibitor")
        except NoSuchElementException:
            logger.warning("Activity select not found — continuing anyway")

        # Paste FASTA content
        try:
            textarea = self.wait.until(EC.presence_of_element_located((By.NAME, "txt_seq")))
            self.driver.execute_script("arguments[0].value = arguments[1];", textarea, fasta_content)
            logger.info("Pasted sequence (%d chars)", len(fasta_content))
        except TimeoutException:
            logger.error("Sequence textarea not found")
            raise

        # Enzymes: pepsin (pH 1.3), trypsin, chymotrypsin (A)
        enzyme_settings = {"enz1": "13", "enz2": "12", "enz3": "11"}
        for name, value in enzyme_settings.items():
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
        """Submit the form and wait for result table to render."""
        assert self.wait is not None
        logger.info("Starting analysis…")
        try:
            btn = self.wait.until(EC.element_to_be_clickable((By.NAME, "but_report")))
            btn.click()
            logger.debug("Clicked report button")
            self.wait.until(
                EC.presence_of_element_located(
                    (
                        By.CSS_SELECTOR,
                        "form[name='fm_action'][action='batch_report_cutting_for_seq'] table.table-out",
                    )
                )
            )
            logger.info("Result tables present")
        except TimeoutException:
            logger.error("Timeout waiting for results")
            raise

    def extract(self) -> List[str]:
        """Scan tables and collect unique sequences with ACE inhibitor activity."""
        assert self.driver is not None
        logger.info("Extracting sequences…")
        tables = self.driver.find_elements(By.CSS_SELECTOR, "form[name='fm_action'] table.table-out")
        if not tables:
            logger.warning("No result tables found")
            return []

        hits: set[str] = set()
        for table in tables:
            rows = table.find_elements(By.CSS_SELECTOR, "tr")
            if not rows:
                continue

            # Heuristic: pick first header row; fall back to first row.
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

            header_cells = header_row.find_elements(By.TAG_NAME, "th") or header_row.find_elements(By.TAG_NAME, "td")
            headers = [" ".join(c.text.split()).strip().lower() for c in header_cells]

            # Locate columns of interest
            seq_idx = next((i for i, h in enumerate(headers) if "sequence" in h), None)
            act_idx = next((i for i, h in enumerate(headers) if "activity" in h), None)
            if seq_idx is None or act_idx is None:
                continue

            for row in rows[header_idx + 1 :]:
                cells = row.find_elements(By.TAG_NAME, "td")
                if not cells or len(cells) <= max(seq_idx, act_idx):
                    continue
                activity = cells[act_idx].text.strip().lower()
                if "ace inhibitor" in activity:
                    seq = cells[seq_idx].text.strip()
                    if seq:
                        hits.add(seq)

        self.results = sorted(hits)
        logger.info("Found %d ACE‑inhibitory sequences", len(self.results))
        return self.results

    def save(self) -> Path:
        if not self.results:
            logger.warning("No results to save")
            # still produce an empty file for traceability
        species = self.cfg.fasta_path.stem
        ts = datetime.now().strftime("%Y%m%d%H%M%S")
        out_path = self.cfg.out_dir / f"biopep_results_{species}_{ts}.txt"
        out_path.write_text("\n".join(self.results) + ("\n" if self.results else ""), encoding="utf-8")
        logger.info("Saved results to %s", out_path)
        return out_path

    def close_driver(self) -> None:
        if self.driver:
            self.driver.quit()
            logger.debug("Browser closed")


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run BIOPEP UWM batch analysis and extract ACE inhibitors.")
    p.add_argument("--fasta", type=Path, required=True, help="Path to input FASTA file")
    p.add_argument("--out-dir", type=Path, default=Path("results"), help="Directory to write outputs")
    p.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT, help="Timeout (s) for page loads and waits")
    p.add_argument("--headless", action="store_true", help="Run Chrome in headless mode")
    p.add_argument("--verbose", action="store_true", help="Enable debug logs")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # logging setup
    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    cfg = RunConfig(
        fasta_path=args.fasta.expanduser().resolve(),
        out_dir=args.out_dir.expanduser().resolve(),
        headless=args.headless,
        timeout=int(args.timeout),
    )
    cfg.validate()

    try:
        with BIOPEPAnalyzer(cfg) as runner:
            runner.navigate()
            runner.configure()
            runner.run()
            runner.extract()
            out_file = runner.save()
        logger.info("Done. Output: %s", out_file)
        return 0
    except Exception as e:
        logger.exception("Failed: %s", e)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
