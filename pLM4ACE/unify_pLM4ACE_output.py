#!/usr/bin/env python3
"""
Summarize high‑activity peptide sequences from Excel files

Purpose
-------
Scans a folder of `.xlsx` result files, extracts sequences with "high activity"
labels, and writes them into a single summary TXT file.

Quick start
-----------
    python summarize_high_activity_sequences.py \
        --input-dir results/excel_outputs \
        --out-dir results/txt_summaries

Outputs
-------
- `<folder_name>.txt` under --out-dir, containing all sequences with "high activity"
  found in the input `.xlsx` files.

Assumptions
-----------
- Each Excel file has at least two columns: `sequence` and `activity`.
- Matching is case‑insensitive (`high activity`).
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
import pandas as pd

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("summarize_high_activity")


def collect_high_activity_sequences(folder: Path) -> list[str]:
    sequences: list[str] = []
    for file_path in folder.glob("*.xlsx"):
        try:
            df = pd.read_excel(file_path)
            if {"sequence", "activity"}.issubset(df.columns):
                hits = df[df["activity"].str.contains("high activity", case=False, na=False)][
                    "sequence"
                ].dropna().astype(str).str.strip()
                sequences.extend(hits.tolist())
                logger.debug("%s → %d hits", file_path.name, len(hits))
            else:
                logger.warning("%s missing required columns", file_path.name)
        except Exception as e:
            logger.error("Error reading %s: %s", file_path.name, e)
    return sequences


def summarize(folder: Path, out_dir: Path) -> Path | None:
    seqs = collect_high_activity_sequences(folder)
    if not seqs:
        logger.warning("No high‑activity sequences found in %s", folder)
        return None

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{folder.name}.txt"
    out_path.write_text("\n".join(seqs) + "\n", encoding="utf-8")
    logger.info("Saved %d sequences → %s", len(seqs), out_path)
    return out_path


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Summarize 'high activity' sequences from Excel files in a folder.")
    p.add_argument("--input-dir", type=Path, required=True, help="Folder containing .xlsx files")
    p.add_argument("--out-dir", type=Path, default=Path("results/txt_summaries"), help="Directory to write the summary file")
    p.add_argument("--verbose", action="store_true", help="Enable debug logs")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    folder = args.input_dir.expanduser().resolve()
    if not folder.is_dir():
        logger.error("Input path is not a folder: %s", folder)
        return 2

    out_dir = args.out_dir.expanduser().resolve()

    out_path = summarize(folder, out_dir)
    if not out_path:
        logger.warning("No output file created.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
