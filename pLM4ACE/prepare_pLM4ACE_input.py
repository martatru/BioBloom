#!/usr/bin/env python3
"""
Format peptide list for pLM4ACE input

Purpose
-------
Takes a plain TXT file with peptide sequences (one per line) and reformats it
into a FASTA‑style file where each sequence is prefixed with "> peptideX".

Quick start
-----------
    python format_peptides_for_plm4ace.py \
        --input data/peptides.txt \
        --out-dir results/pLM4ACE_input

Outputs
-------
- `<species_name>_formatted.txt` containing:

      > peptide1
      AAAAK
      > peptide2
      GGFL
      ...

Notes
-----
- Species name is inferred from the input filename (without extension).
- Empty lines in the input are skipped.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("plm4ace.formatter")


def format_peptides(input_file: Path, out_dir: Path) -> Path:
    if not input_file.exists() or not input_file.is_file():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    species_name = input_file.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    output_file = out_dir / f"{species_name}_formatted.txt"

    lines = input_file.read_text(encoding="utf-8").splitlines()
    formatted: list[str] = []
    for idx, line in enumerate(lines, start=1):
        seq = line.strip()
        if not seq:
            continue
        formatted.append(f"> peptide{idx}\n")
        formatted.append(seq + "\n")

    output_file.write_text("".join(formatted), encoding="utf-8")
    logger.info("Saved %d peptides to %s", len(formatted) // 2, output_file)
    return output_file


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Format plain peptide list into FASTA‑style for pLM4ACE input.")
    p.add_argument("--input", type=Path, required=True, help="Path to TXT file with peptide sequences")
    p.add_argument("--out-dir", type=Path, default=Path("results/pLM4ACE_input"), help="Output directory")
    p.add_argument("--verbose", action="store_true", help="Enable debug logs")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        input_file = args.input.expanduser().resolve()
        out_dir = args.out_dir.expanduser().resolve()

        out_path = format_peptides(input_file, out_dir)
        logger.info("Done. Output file: %s", out_path)
        return 0
    except Exception as e:
        logger.exception("Formatting failed: %s", e)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
