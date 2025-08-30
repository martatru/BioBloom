#!/usr/bin/env python3
"""
Split FASTA‑style peptide file into N parts

Purpose
-------
Splits a file containing peptide entries (two lines per entry: header + sequence)
into a given number of parts. Each output file preserves the FASTA structure.

Quick start
-----------
    python split_fasta_refactored.py \
        --input data/peptides_formatted.txt \
        --num-parts 5 \
        --out-dir results/pLM4ACE_input/splits

Outputs
-------
- `split_part_1.txt`, `split_part_2.txt`, … in the target directory
- Each file contains roughly equal number of entries; remainder distributed across first parts.

Assumptions
-----------
- Input file has an even number of lines (header and sequence alternating).
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("split_fasta")


def split_fasta(input_path: Path, output_dir: Path, num_parts: int) -> list[Path]:
    lines = input_path.read_text(encoding="utf-8").splitlines()
    if len(lines) % 2 != 0:
        raise ValueError("Input file must have pairs of lines (header + sequence)")

    # Group entries
    entries = [(lines[i].strip(), lines[i + 1].strip()) for i in range(0, len(lines), 2)]
    total_entries = len(entries)
    if num_parts <= 0:
        raise ValueError("Number of parts must be positive")
    if num_parts > total_entries:
        raise ValueError("Number of parts cannot exceed number of entries")

    part_size = total_entries // num_parts
    remainder = total_entries % num_parts

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    start = 0
    for i in range(num_parts):
        end = start + part_size + (1 if i < remainder else 0)
        part_entries = entries[start:end]
        start = end

        out_path = output_dir / f"split_part_{i+1}.txt"
        with out_path.open("w", encoding="utf-8") as f:
            for header, seq in part_entries:
                f.write(f"{header}\n{seq}\n")
        logger.info("Created %s with %d entries", out_path, len(part_entries))
        outputs.append(out_path)

    return outputs


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Split FASTA‑style peptide file into N parts.")
    p.add_argument("--input", type=Path, required=True, help="Path to peptide TXT/FASTA file")
    p.add_argument("--num-parts", type=int, required=True, help="Number of parts to split into")
    p.add_argument("--out-dir", type=Path, default=Path("results/pLM4ACE_input/splits"), help="Output directory")
    p.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        input_path = args.input.expanduser().resolve()
        out_dir = args.out_dir.expanduser().resolve()

        if not input_path.exists():
            logger.error("Input file not found: %s", input_path)
            return 2

        outputs = split_fasta(input_path, out_dir, args.num_parts)
        logger.info("Done. Created %d files in %s", len(outputs), out_dir)
        return 0
    except Exception as e:
        logger.exception("Splitting failed: %s", e)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
