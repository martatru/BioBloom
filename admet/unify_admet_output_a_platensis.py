#!/usr/bin/env python3
"""
Merge multiple CSV files (same header) into one file

Purpose
-------
Concatenate all CSV files in a folder (matching a pattern) into a single
`merged.csv` while writing the header only once.

This version is CLI‑driven (no interactive prompts), uses Pathlib + logging,
handles UTF‑8 with BOM (`utf-8-sig`) gracefully, and skips empty files.

Quick start
-----------
    python merge_csvs_refactored.py \
        --input-dir data/csvs \
        --pattern "*.csv" \
        --output merged.csv

Options
-------
- `--input-dir`   Directory containing source CSVs
- `--pattern`     Glob pattern for selecting files (default: `*.csv`)
- `--output`      Output path (default: `<input-dir>/merged.csv`)
- `--encoding`    Input encoding (default: `utf-8-sig`)
- `--out-encoding` Output encoding (default: `utf-8`)
- `--strict-headers`  Fail if a file's header differs from the first header

Notes
-----
- Files are processed in sorted order by filename.
- The output file is automatically excluded from the inputs if it lives in the
  same directory.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("merge_csvs")


def iter_csv_paths(folder: Path, pattern: str, output_path: Path) -> Iterable[Path]:
    for p in sorted(folder.glob(pattern)):
        # Exclude the output file itself to avoid self-append
        if p.resolve() == output_path.resolve():
            continue
        yield p


def merge_csvs(
    folder: Path,
    pattern: str = "*.csv",
    output: Path | None = None,
    encoding: str = "utf-8-sig",
    out_encoding: str = "utf-8",
    strict_headers: bool = False,
) -> tuple[Path, int, int]:
    """Merge CSVs. Returns (output_path, files_merged, lines_written)."""
    if not folder.is_dir():
        raise NotADirectoryError(f"Not a directory: {folder}")

    out_path = output or (folder / "merged.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    wrote_header = False
    lines_written = 0
    files_merged = 0
    first_header: str | None = None

    csv_paths = list(iter_csv_paths(folder, pattern, out_path))
    if not csv_paths:
        logger.info("No files matched %s in %s", pattern, folder)
        # still create empty file for reproducibility
        out_path.write_text("", encoding=out_encoding)
        return out_path, 0, 0

    with out_path.open("w", encoding=out_encoding, newline="") as out_f:
        for path in csv_paths:
            with path.open("r", encoding=encoding, newline="") as f:
                header = f.readline()
                if header == "":  # empty file
                    logger.debug("Skipping empty file: %s", path.name)
                    continue

                # Normalize header line endings/whitespace minimally
                header_norm = header.strip("\n\r")
                if not wrote_header:
                    out_f.write(header)
                    lines_written += 1
                    wrote_header = True
                    first_header = header_norm
                else:
                    if strict_headers and header_norm != first_header:
                        raise ValueError(
                            f"Header mismatch in {path.name}. Expected '{first_header}', got '{header_norm}'."
                        )

                for line in f:
                    out_f.write(line)
                    lines_written += 1
                files_merged += 1

    if not wrote_header:
        # All files were empty
        out_path.write_text("", encoding=out_encoding)
        logger.warning("All CSV files were empty. Created an empty %s", out_path.name)
    else:
        logger.info(
            "Merged %d file(s) into %s (total lines incl. header: %d)",
            files_merged,
            out_path,
            lines_written,
        )

    return out_path, files_merged, lines_written


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Merge CSV files in a folder into a single CSV with one header.")
    p.add_argument("--input-dir", type=Path, required=True, help="Folder containing CSV files")
    p.add_argument("--pattern", type=str, default="*.csv", help="Glob pattern for CSV selection (default: *.csv)")
    p.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Path to output CSV (default: <input-dir>/merged.csv)",
    )
    p.add_argument("--encoding", type=str, default="utf-8-sig", help="Input file encoding (default: utf-8-sig)")
    p.add_argument("--out-encoding", type=str, default="utf-8", help="Output file encoding (default: utf-8)")
    p.add_argument("--strict-headers", action="store_true", help="Fail if a file has a different header")
    p.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    folder = args.input_dir.expanduser().resolve()
    output = args.output.expanduser().resolve() if args.output else None

    try:
        out_path, files_merged, lines_written = merge_csvs(
            folder=folder,
            pattern=args.pattern,
            output=output,
            encoding=args.encoding,
            out_encoding=args.out_encoding,
            strict_headers=args.strict_headers,
        )
        logger.info(
            "Done. Merged %d file(s) into %s (lines incl. header: %d)",
            files_merged,
            out_path,
            lines_written,
        )
        return 0
    except Exception as e:
        logger.exception("Merge failed: %s", e)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
