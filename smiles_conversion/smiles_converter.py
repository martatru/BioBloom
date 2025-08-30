#!/usr/bin/env python3
"""
Convert between FASTA and SMILES (CLI, portable paths)

Purpose
-------
Two-way conversion with pluggable backends:
- FASTA → SMILES via an external CLI (default: `fasta2smi` from p2smi)
- SMILES → FASTA via RDKit heuristics (recover peptide sequences)

This refactor removes interactive `input()` calls, adds `argparse`, logging,
safer path handling, and explicit overwrite control.

Quick start
-----------
# FASTA → SMILES
    python convert_fasta_smiles_refactored.py \
      --mode f2s \
      --input data/peptides.fasta \
      --out-dir results/converted \
      --fasta2smi-cmd fasta2smi

# SMILES → FASTA
    python convert_fasta_smiles_refactored.py \
      --mode s2f \
      --input data/peptides.smi \
      --out-dir results/converted

Notes
-----
- For FASTA→SMILES you must have an external converter installed and available
  on PATH (default tool name: `fasta2smi`). You can override with
  `--fasta2smi-cmd /full/path/to/tool`.
- RDKit is required for SMILES→FASTA. If unavailable, the script exits with a
  clear error.
"""
from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("conv.fasta_smiles")


# -----------------------------------------------------------------------------
# CLI config
# -----------------------------------------------------------------------------
@dataclass
class RunConfig:
    mode: str  # "f2s" or "s2f" (or "auto")
    input_path: Path
    out_dir: Path
    overwrite: bool
    fasta2smi_cmd: str

    def resolve(self) -> None:
        self.input_path = self.input_path.expanduser().resolve()
        self.out_dir = self.out_dir.expanduser().resolve()
        self.out_dir.mkdir(parents=True, exist_ok=True)

    def infer_mode_if_auto(self) -> None:
        if self.mode != "auto":
            return
        ext = self.input_path.suffix.lower()
        if ext in {".fa", ".fasta", ".faa", ".fna", ".fas"}:
            self.mode = "f2s"
        elif ext in {".smi", ".smiles", ".txt"}:
            self.mode = "s2f"
        else:
            raise ValueError(f"Cannot infer mode from extension: {ext}")


# -----------------------------------------------------------------------------
# FASTA → SMILES via external CLI
# -----------------------------------------------------------------------------

def fasta_to_smiles(fasta_path: Path, out_smi: Path, tool: str) -> None:
    if shutil.which(tool) is None:
        raise FileNotFoundError(
            f"Required converter '{tool}' not found on PATH. Provide with --fasta2smi-cmd or install it."
        )
    cmd = [tool, "-i", str(fasta_path), "-o", str(out_smi)]
    logger.info("Running: %s", " ".join(cmd))
    res = subprocess.run(cmd, check=False)
    if res.returncode != 0:
        raise RuntimeError(f"{tool} failed with exit code {res.returncode}")
    logger.info("Saved SMILES → %s", out_smi)


# -----------------------------------------------------------------------------
# SMILES → FASTA via RDKit (heuristic sequence recovery)
# -----------------------------------------------------------------------------

def _load_rdkit():
    try:
        from rdkit import Chem  # type: ignore
        from rdkit import RDLogger  # type: ignore
        RDLogger.DisableLog("rdApp.*")
        return Chem
    except Exception as e:
        raise RuntimeError("RDKit not available. Install rdkit to enable SMILES→FASTA.") from e

# Canonical amino acid residue motifs (simplified heuristics)
AA_SMILES = {
    "ALA": "C[C@H](N)C=O",
    "CYS": "N[C@H](C=O)CS",
    "ASP": "N[C@H](C=O)CC(=O)O",
    "GLU": "N[C@H](C=O)CCC(=O)O",
    "PHE": "N[C@H](C=O)Cc1ccccc1",
    "GLY": "NCC=O",
    "HIS": "N[C@H](C=O)Cc1c[nH]cn1",
    "ILE": "CC[C@H](C)[C@H](N)C=O",
    "LYS": "NCCCC[C@H](N)C=O",
    "LEU": "CC(C)C[C@H](N)C=O",
    "MET": "CSCC[C@H](N)C=O",
    "ASN": "NC(=O)C[C@H](N)C=O",
    "PRO": "O=C[C@@H]1CCCN1",
    "GLN": "NC(=O)CC[C@H](N)C=O",
    "ARG": "N=C(N)NCCC[C@H](N)C=O",
    "SER": "N[C@H](C=O)CO",
    "THR": "C[C@@H](O)[C@H](N)C=O",
    "VAL": "CC(C)[C@H](N)C=O",
    "TRP": "N[C@H](C=O)Cc1c[nH]c2ccccc12",
    "TYR": "N[C@H](C=O)Cc1ccc(O)cc1",
}
AA_ORDER = [
    "GLY",
    "ALA",
    "VAL",
    "CYS",
    "ASP",
    "GLU",
    "PHE",
    "HIS",
    "ILE",
    "LYS",
    "LEU",
    "MET",
    "ASN",
    "PRO",
    "GLN",
    "ARG",
    "SER",
    "THR",
    "TRP",
    "TYR",
]


def _mol_to_sequence_via_rdkit(Chem, mol) -> Optional[str]:
    if mol is None:
        return None

    # Mark C-alpha atoms so Chem.MolToSequence has residue hints
    CA_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[C:0](=[O:1])[C:2][N:3]"))
    for atoms in CA_matches:
        a = mol.GetAtomWithIdx(atoms[2])
        info = Chem.AtomPDBResidueInfo()
        info.SetName(" CA ")
        a.SetMonomerInfo(info)

    # Assign residue names by matching AA templates
    for aa in AA_ORDER:
        tmpl = Chem.MolFromSmiles(AA_SMILES[aa])
        for match in mol.GetSubstructMatches(tmpl):
            for idx in match:
                a = mol.GetAtomWithIdx(idx)
                mi = a.GetMonomerInfo()
                if mi is not None and mi.GetName() == " CA ":
                    info = Chem.AtomPDBResidueInfo()
                    info.SetName(" CA ")
                    info.SetResidueName(aa)
                    a.SetMonomerInfo(info)

    # Try renumbering to get a contiguous backbone, then extract sequence
    try:
        seq = Chem.MolToSequence(mol)
        return seq if seq else None
    except Exception:
        return None


def _parse_smi_line(line: str) -> Tuple[str, str]:
    parts = line.strip().split()
    if not parts:
        return "", ""
    smi = parts[0]
    name = parts[1] if len(parts) > 1 else ""
    return smi, name


def smiles_to_fasta(smi_path: Path, out_fasta: Path) -> None:
    Chem = _load_rdkit()
    ok, fail = 0, 0
    with smi_path.open("r", encoding="utf-8", errors="ignore") as fin, out_fasta.open(
        "w", encoding="utf-8"
    ) as fout:
        for i, line in enumerate(fin, 1):
            if not line.strip():
                continue
            smi, name = _parse_smi_line(line)
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            seq = _mol_to_sequence_via_rdkit(Chem, mol)
            if seq:
                header = name if name else f"mol_{i}"
                fout.write(f">{header}\n{seq}\n")
                ok += 1
            else:
                fail += 1
                logger.warning("line %d: could not derive sequence", i)
    logger.info("SMILES→FASTA: wrote %d sequences (%d failed) → %s", ok, fail, out_fasta)


# -----------------------------------------------------------------------------
# Orchestration
# -----------------------------------------------------------------------------

def run(cfg: RunConfig) -> Path:
    cfg.resolve()
    cfg.infer_mode_if_auto()

    base = cfg.input_path.stem
    if cfg.mode == "f2s":
        out_path = cfg.out_dir / f"{base}.smi"
        if out_path.exists() and not cfg.overwrite:
            raise FileExistsError(f"Output exists: {out_path}. Use --overwrite to replace.")
        fasta_to_smiles(cfg.input_path, out_path, cfg.fasta2smi_cmd)
        return out_path

    if cfg.mode == "s2f":
        out_path = cfg.out_dir / f"{base}.fasta"
        if out_path.exists() and not cfg.overwrite:
            raise FileExistsError(f"Output exists: {out_path}. Use --overwrite to replace.")
        smiles_to_fasta(cfg.input_path, out_path)
        return out_path

    raise ValueError(f"Unknown mode: {cfg.mode}")


# -----------------------------------------------------------------------------
# CLI entrypoint
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Convert between FASTA and SMILES.")
    p.add_argument("--mode", choices=["f2s", "s2f", "auto"], default="auto", help="Conversion direction")
    p.add_argument("--input", type=Path, required=True, help="Path to input file (.fasta/.fa or .smi)")
    p.add_argument("--out-dir", type=Path, default=Path("results/converted"), help="Where to save output")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing output file if present")

    # External tool for FASTA→SMILES
    p.add_argument(
        "--fasta2smi-cmd",
        type=str,
        default="fasta2smi",
        help="Command/binary to use for FASTA→SMILES conversion",
    )

    p.add_argument("--verbose", action="store_true", help="Enable debug logs")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    cfg = RunConfig(
        mode=args.mode,
        input_path=args.input,
        out_dir=args.out_dir,
        overwrite=args.overwrite,
        fasta2smi_cmd=args.fasta2smi_cmd,
    )

    try:
        out_path = run(cfg)
        logger.info("Done. Output: %s", out_path)
        return 0
    except Exception as e:
        logger.exception("Conversion failed: %s", e)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
