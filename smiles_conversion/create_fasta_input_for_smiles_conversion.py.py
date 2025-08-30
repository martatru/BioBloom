#!/usr/bin/env python3
"""
Convert peptide list (TXT) into FASTA file for OpenBabel input

Purpose
-------
Reads a plain TXT file with one peptide sequence per line and outputs a FASTA
file with headers in the form `> <basename>_pepN`.

Quick start
-----------
    python peptides_to_fasta_refactored.py \
        --input data/peptides.txt \
        --out-dir results/OPENBABEL_INPUT

Outputs
-------
- `<basename>.fasta` saved under --out-dir
- Each sequence prefixed with `> <basename>_pepN`

Example FASTA output
--------------------
    > sample_pep1
    AAAAK
    > sample_pep2
    GGFL

"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("peptides_to_fasta")


def peptides_to_fasta(input_txt: Path, out_dir: Path) -> Path:
    if not input_txt.exists():
        raise FileNotFoundError(f"Input file not found: {input_txt}")

    peptides = [line.strip() for line in input_txt.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not peptides:
        raise ValueError("No peptides found in input file")

    out_dir.mkdir(parents=True, exist_ok=True)
    base_name = input_txt.stem
    out_fasta = out_dir / f"{base_name}.fasta"

    with out_fasta.open("w", encoding="utf-8") as f:
        for i, pep in enumerate(peptides, start=1):
            f.write(f">{base_name}_pep{i}\n")
            f.write(pep + "\n")

    logger.info("Wrote FASTA with %d peptides → %s", len(peptides), out_fasta)
    return out_fasta


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Convert plain peptide TXT file into FASTA format for OpenBabel input.")
    p.add_argument("--input", type=Path, required=True, help="TXT file with peptide sequences (one per line)")
    p.add_argument("--out-dir", type=Path, default=Path("results/OPENBABEL_INPUT"), help="Directory to save FASTA output")
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
        out_dir = args.out_dir.expanduser().resolve()

        out_path = peptides_to_fasta(input_txt, out_dir)
        logger.info("Done. Output: %s", out_path)
        return 0
    except Exception as e:
        logger.exception("Conversion failed: %s", e)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
