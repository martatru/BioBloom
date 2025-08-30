#!/usr/bin/env python3
"""
Generate peptide PDBs from a FASTA (PyRosetta)

Purpose
-------
Reads a (multi)FASTA of peptide sequences and writes a PDB per entry using
PyRosetta's `pose_from_sequence` with the full-atom residue type set.

Quick start
-----------
    python generate_peptide_structures_pyrosetta_refactored.py \
      /path/to/peptides.fasta \
      --outdir pdbs \
      --minlen 2 --maxlen 50 \
      --silent

Notes
-----
- Requires a working PyRosetta installation and license.
- Only 20 canonical L-amino acids are considered valid by default. Non-standard
  letters are flagged; use `--allow-letters` to whitelist extra codes (e.g. `U`).
"""
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Iterable, List, Tuple

from pyrosetta import init, pose_from_sequence

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("pyrosetta.pepgen")

VALID_AA_DEFAULT = set("ACDEFGHIKLMNPQRSTVWY")  # standard L-amino acids (1-letter)


# -----------------------------------------------------------------------------
# FASTA parsing
# -----------------------------------------------------------------------------

def read_fasta(path: Path) -> List[Tuple[str, str]]:
    """Return list of (header, seq). Accepts files with/without headers."""
    entries: List[Tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    text = path.read_text(encoding="utf-8")
    if not text.strip():
        return []

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_parts)))
            header = line[1:].strip() or "seq"
            seq_parts = []
        else:
            seq_parts.append(line)

    if header is not None:
        entries.append((header, "".join(seq_parts)))

    if not entries:
        # no headers: treat whole file as one sequence
        seq = re.sub(r"\s+", "", text)
        entries = [("seq", seq)]
    return entries


def sanitize_seq_id(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s[:60] if s else "seq"


def normalize_seq(seq: str, allow: set[str]) -> tuple[str, list[str]]:
    """Uppercase sequence and report letters outside allow-set."""
    seq_u = seq.upper()
    bad = sorted({c for c in seq_u if c not in allow})
    return seq_u, bad


# -----------------------------------------------------------------------------
# Core
# -----------------------------------------------------------------------------

def generate_pdbs(
    fasta: Path,
    outdir: Path,
    minlen: int,
    maxlen: int,
    silent: bool,
    allow_letters: set[str] | None = None,
    overwrite: bool = True,
) -> int:
    allow = set(allow_letters) if allow_letters else set(VALID_AA_DEFAULT)

    # Initialize PyRosetta
    init_flags = "-mute all" if silent else ""
    logger.info("Initializing PyRosetta (%s)", "silent" if silent else "verbose")
    init(init_flags)

    entries = read_fasta(fasta)
    if not entries:
        logger.error("No sequences found in %s", fasta)
        return 2

    outdir.mkdir(parents=True, exist_ok=True)

    n_ok = 0
    for raw_name, raw_seq in entries:
        name = sanitize_seq_id(raw_name)
        seq, bad = normalize_seq(raw_seq, allow)
        if bad:
            logger.warning("%s: found non-allowed letters %s", name, bad)
        if not seq:
            logger.warning("%s: empty sequence — skipping", name)
            continue
        if not (minlen <= len(seq) <= maxlen):
            logger.info("%s: length %d outside [%d,%d] — skipping", name, len(seq), minlen, maxlen)
            continue

        out_path = outdir / f"{name}.pdb"
        if out_path.exists() and not overwrite:
            logger.info("%s exists — skipping (use --overwrite to replace)", out_path.name)
            continue

        try:
            pose = pose_from_sequence(seq, "fa_standard")
            pose.dump_pdb(str(out_path))
            logger.info("OK %s: %d aa -> %s", name, len(seq), out_path)
            n_ok += 1
        except Exception as e:
            logger.exception("FAIL %s: %s", name, e)

    logger.info("Done. Wrote %d PDB(s) to %s", n_ok, outdir)
    return 0 if n_ok > 0 else 1


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Convert FASTA peptides to PDB using PyRosetta.")
    p.add_argument("fasta", type=Path, help="Input FASTA file (multi-FASTA supported)")
    p.add_argument("-o", "--outdir", type=Path, default=Path("out_pdbs"), help="Output directory")
    p.add_argument("--minlen", type=int, default=2, help="Minimum length to keep")
    p.add_argument("--maxlen", type=int, default=50, help="Maximum length to keep")
    p.add_argument("--allow-letters", type=str, default="", help="Extra allowed letters (e.g., 'UOZ')")
    p.add_argument("--silent", action="store_true", help="Mute Rosetta output (-mute all)")
    p.add_argument("--no-overwrite", action="store_true", help="Do not overwrite existing PDBs")
    p.add_argument("--verbose", action="store_true", help="Enable debug logs")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    fasta = args.fasta.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()

    if not fasta.exists():
        logger.error("FASTA not found: %s", fasta)
        return 2

    allow = set(c for c in args.allow_letters.upper()) if args.allow_letters else None
    rc = generate_pdbs(
        fasta=fasta,
        outdir=outdir,
        minlen=args.minlen,
        maxlen=args.maxlen,
        silent=args.silent,
        allow_letters=allow,
        overwrite=not args.no_overwrite,
    )
    return rc


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
