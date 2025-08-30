#!/usr/bin/env python3
"""
Extract SMILES strings from text file

Purpose
-------
Reads a text file containing colon-separated entries (e.g. "id: SMILES")
and writes only the SMILES strings to a new `.smi` file.

Quick start
-----------
    python extract_smiles.py \
        --input data/raw_smiles.txt \
        --out-dir results/openbabel_output

Outputs
-------
- `<basename>2.smi` saved under --out-dir
- Each line contains one SMILES string.
"""

from __future__ import annotations
import argparse
import logging
from pathlib import Path

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("extract_smiles")


def extract_smiles(inp: Path, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{inp.stem}2.smi"

    n_ok = n_bad = 0
    with inp.open("r", encoding="utf-8") as fin, out_file.open("w", encoding="utf-8") as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue

            # take only the part after colon (SMILES)
            smiles = line.split(":", 1)[1].strip() if ":" in line else line
            if smiles:
                fout.write(smiles + "\n")
                n_ok += 1
            else:
                n_bad += 1

    logger.info("Extracted %d SMILES → %s", n_ok, out_file)
    if n_bad:
        logger.warning("Skipped %d invalid/empty lines", n_bad)
    return out_file


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Extract SMILES from colon-separated text file.")
    p.add_argument("--input", type=Path, required=True, help="Input TXT file with SMILES (optionally colon-separated)")
    p.add_argument("--out-dir", type=Path, default=Path("results/openbabel_output"), help="Directory for .smi output")
    p.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    inp = args.input.expanduser().resolve()
    if not inp.exists():
        logger.error("Input file not found: %s", inp)
        return 2

    out_dir = args.out_dir.expanduser().resolve()
    try:
        out_file = extract_smiles(inp, out_dir)
        logger.info("Done. Output: %s", out_file)
        return 0
    except Exception as e:
        logger.exception("Extraction failed: %s", e)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
