#!/usr/bin/env python3
"""
Merge multiple peptide list files and remove duplicates

Purpose
-------
Combines peptide TXT files (one sequence per line), drops duplicates,
optionally sorts entries, and writes a unified list.

Quick start
-----------
    # Option A: pass explicit files
    python unify_peptide_lists_refactored.py \
      --inputs data/part1.txt data/part2.txt data/part3.txt data/part4.txt \
      --out-dir results/aplatensis_unified \
      --output unified_peptides.txt

    # Option B: merge every .txt in a folder
    python unify_peptide_lists_refactored.py \
      --input-dir data/splits --pattern "*.txt" \
      --out-dir results/aplatensis_unified

Outputs
-------
- Unified TXT file under --out-dir (default name: unified_peptides.txt)

Notes
-----
- Duplicates are removed with exact string comparison after trimming whitespace.
- Use --no-sort to preserve discovery order (folder glob order).
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("unify_peptides")


def iter_input_files(inputs: list[Path] | None, input_dir: Path | None, pattern: str) -> Iterable[Path]:
    if inputs:
        for p in inputs:
            yield p
    elif input_dir:
        yield from sorted(input_dir.glob(pattern))
    else:
        raise ValueError("Provide --inputs or --input-dir")


def unify_peptides(
    inputs: list[Path] | None,
    input_dir: Path | None,
    pattern: str,
    out_dir: Path,
    output_name: str,
    sort_output: bool = True,
) -> tuple[Path, int]:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / output_name

    seen: set[str] = set()
    ordered: list[str] = []

    n_files = 0
    for path in iter_input_files(inputs, input_dir, pattern):
        if not path.exists() or not path.is_file():
            logger.warning("Skipping missing file: %s", path)
            continue
        n_files += 1
        for line in path.read_text(encoding="utf-8").splitlines():
            pep = line.strip()
            if not pep:
                continue
            if pep not in seen:
                seen.add(pep)
                ordered.append(pep)

    if sort_output:
        ordered = sorted(ordered)

    out_path.write_text("\n".join(ordered) + ("\n" if ordered else ""), encoding="utf-8")
    logger.info("Unified %d file(s) → %s (unique peptides: %d)", n_files, out_path, len(ordered))
    return out_path, len(ordered)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Merge peptide TXT files, drop duplicates, and write a unified list.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--inputs", type=Path, nargs="+", help="Explicit input files")
    src.add_argument("--input-dir", type=Path, help="Folder containing input files")
    p.add_argument("--pattern", type=str, default="*.txt", help="Glob when using --input-dir (default: *.txt)")
    p.add_argument("--out-dir", type=Path, default=Path("results/unified_peptides"), help="Output directory")
    p.add_argument("--output", type=str, default="unified_peptides.txt", help="Output filename")
    p.add_argument("--no-sort", action="store_true", help="Do not sort output (preserve discovery order)")
    p.add_argument("--verbose", action="store_true", help="Enable debug logs")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        inputs = [p.expanduser().resolve() for p in (args.inputs or [])]
        input_dir = args.input_dir.expanduser().resolve() if args.input_dir else None
        out_dir = args.out_dir.expanduser().resolve()

        out_path, n_unique = unify_peptides(
            inputs=inputs or None,
            input_dir=input_dir,
            pattern=args.pattern,
            out_dir=out_dir,
            output_name=args.output,
            sort_output=not args.no_sort,
        )
        logger.info("Done. Unique peptides: %d", n_unique)
        return 0
    except Exception as e:
        logger.exception("Unification failed: %s", e)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
