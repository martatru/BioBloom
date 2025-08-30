#!/usr/bin/env python3
"""
Repack side chains of a receptor PDB with PyRosetta

Purpose
-------
Loads a receptor PDB, performs side-chain repacking (no mutations) using
standard full-atom scoring, and saves the repacked structure.

Quick start
-----------
    python repack_receptor.py \
      --input /path/to/ACE_clean.pdb \
      --output /path/to/ACE_repacked.pdb \
      --init-flags "-mute all -ex1 -ex2aro -use_input_sc"

Notes
-----
- Requires PyRosetta installed and licensed.
- By default enables Dunbrack extra rotamers (-ex1 -ex2aro) and uses input side chains.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from pyrosetta import init, pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    RestrictToRepacking,
    InitializeFromCommandline,
    IncludeCurrent,
)
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("pyrosetta.repack")


def repack_sidechains(input_pdb: Path, output_pdb: Path) -> None:
    """Repack side chains only; no mutations."""
    logger.info("Loading pose from %s", input_pdb)
    pose = pose_from_pdb(str(input_pdb))

    scorefxn = get_fa_scorefxn()
    logger.debug("Score function: %s", scorefxn.get_name())

    tf = TaskFactory()
    tf.push_back(InitializeFromCommandline())  # honor command-line init flags
    tf.push_back(IncludeCurrent())             # keep current rotamers as options
    tf.push_back(RestrictToRepacking())        # disallow design (no mutations)

    task = tf.create_task_and_apply_taskoperations(pose)
    packer = PackRotamersMover(scorefxn, task)

    logger.info("Applying PackRotamersMover (repacking only)")
    packer.apply(pose)

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    pose.dump_pdb(str(output_pdb))
    logger.info("Saved repacked receptor → %s", output_pdb)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Repack receptor side chains with PyRosetta (no mutations)."
    )
    p.add_argument("--input", type=Path, required=True, help="Path to input PDB")
    p.add_argument("--output", type=Path, required=True, help="Path to output PDB")
    p.add_argument(
        "--init-flags",
        type=str,
        default="-mute all -ex1 -ex2aro -use_input_sc",
        help="Flags passed to PyRosetta init (quoted string)",
    )
    p.add_argument("--verbose", action="store_true", help="Enable debug logs")
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    inp = args.input.expanduser().resolve()
    outp = args.output.expanduser().resolve()

    if not inp.exists():
        logger.error("Input PDB not found: %s", inp)
        return 2

    # Initialize PyRosetta with user-provided flags (e.g., -ex1 -ex2aro)
    logger.info("Initializing PyRosetta with flags: %s", args.init_flags)
    init(args.init_flags)

    try:
        repack_sidechains(inp, outp)
        logger.info("Done.")
        return 0
    except Exception as e:
        logger.exception("Repacking failed: %s", e)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
