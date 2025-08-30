#!/usr/bin/env python3
"""
Novel peptide comparison utility

Purpose
-------
Compares a list of peptides (TXT) against a reference set (ODS) and outputs
potentially novel peptides that are not present in the reference.

This version is CLI‑driven, uses Pathlib, and writes results to a user‑defined
output directory.

Quick start
-----------
    python novel_peptides_compare_refactored.py \
        --input peptides.txt \
        --reference data/known_ace_inhibitory_peptides.ods \
        --out-dir results/novel_peptides

Outputs
-------
- `<species>_potentially_novel_peptides.txt` under --out-dir
- One peptide per line, sorted alphabetically.

Notes
-----
- Requires pandas with ODF support (engine="odf").
- Reference ODS: peptides in the **first column**, header row ignored.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
import pandas as pd

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("novel_peptides")


def compare_peptides(input_txt: Path, reference_ods: Path, out_dir: Path) -> Path:
    # Load reference peptides from first column
    try:
        known_df = pd.read_excel(reference_ods, engine="odf")
    except Exception as e:
        logger.error("Failed to read reference ODS: %s", e)
        raise

    known_peptides = set(known_df.iloc[:, 0].dropna().astype(str).str.strip())
    logger.info("Loaded %d known peptides from %s", len(known_peptides), reference_ods)

    # Load input peptides from TXT
    input_peptides = {
        line.strip()
        for line in input_txt.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    logger.info("Loaded %d input peptides from %s", len(input_peptides), input_txt)

    # Compute novel peptides
    novel_peptides = sorted(input_peptides - known_peptides)
    logger.info("Identified %d potentially novel peptides", len(novel_peptides))

    # Prepare output
    out_dir.mkdir(parents=True, exist_ok=True)
    species = input_txt.stem
    out_path = out_dir / f"{species}_potentially_novel_peptides.txt"
    out_path.write_text("\n".join(novel_peptides) + ("\n" if novel_peptides else ""), encoding="utf-8")
    logger.info("Saved results to %s", out_path)
    return out_path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Compare peptides against reference set and list novel entries.")
    p.add_argument("--input", type=Path, required=True, help="TXT file with peptides (one per line)")
    p.add_argument("--reference", type=Path, required=True, help="ODS file with known peptides (first column)")
    p.add_argument("--out-dir", type=Path, default=Path("results/novel_peptides"), help="Output directory")
    p.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        input_txt = args.input.expanduser().resolve()
        reference_ods = args.reference.expanduser().resolve()
        out_dir = args.out_dir.expanduser().resolve()

        if not input_txt.exists():
            logger.error("Input file not found: %s", input_txt)
            return 2
        if not reference_ods.exists():
            logger.error("Reference file not found: %s", reference_ods)
            return 2

        out_file = compare_peptides(input_txt, reference_ods, out_dir)
        logger.info("Done. Output file: %s", out_file)
        return 0
    except Exception as e:
        logger.exception("Comparison failed: %s", e)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
