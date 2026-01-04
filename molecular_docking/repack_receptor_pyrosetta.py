#!/usr/bin/env python3
"""
Prepare ACE-I like in the paper:
- Remove lisinopril (LPR) and waters (HOH)
- Keep Zn2+ by default (toggle with --remove-zinc)
- Add polar H/charges via Rosetta chemistry
- Rebuild incomplete side chains using Dunbrack (use_input_sc false + dense rotamers)
- Optional short constrained relax
"""

from __future__ import annotations
import argparse, logging
from pathlib import Path

from pyrosetta import init, pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    RestrictToRepacking,
    InitializeFromCommandline,
)
from pyrosetta.rosetta.core.pack.task.operation import ExtraRotamersGeneric
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
logger = logging.getLogger("pyrosetta.repack")

def strip_lpr_water_keep_zn(pose, keep_zinc: bool = True):
    """Remove LPR and HOH. Keep Zn unless --remove-zinc is used. Remove other hetero."""
    idx_to_del = []
    for i in range(1, pose.size()+1):
        r = pose.residue(i)
        n3 = r.name3()
        if n3 == 'HOH' or n3 == 'LPR':
            idx_to_del.append(i)
        elif not r.is_protein():
            if keep_zinc and n3 in ('ZN','ZN2','ZN+'):
                continue
            idx_to_del.append(i)
    for i in reversed(idx_to_del):
        pose.delete_residue_slow(i)
    logger.info("Removed %d residues (LPR/HOH/hetero%s).",
                len(idx_to_del),
                ", kept Zn" if keep_zinc else ", removed Zn")

def likely_missing_sc_heavy_atoms(pose):
    """Crude detector: residues without CB / chi1-heavy atoms."""
    missing = []
    for i in range(1, pose.size()+1):
        rsd = pose.residue(i)
        if not rsd.is_protein():
            continue
        if rsd.name1() == 'G':  # gly: no CB by design
            continue
        if not rsd.has('CB'):
            missing.append(i); continue
        try:
            if rsd.nchi() > 0:
                _ = rsd.chi(1)  # throws if chi atoms absent
        except Exception:
            missing.append(i)
    return missing

def prepare(input_pdb: Path, output_pdb: Path, keep_zinc: bool, do_relax: bool):
    logger.info("Loading pose from %s", input_pdb)
    pose = pose_from_pdb(str(input_pdb))

    # Clean per paper
    strip_lpr_water_keep_zn(pose, keep_zinc=keep_zinc)

    before = likely_missing_sc_heavy_atoms(pose)
    if before:
        logger.warning("Likely missing side-chain heavy atoms BEFORE: %s", before)
    else:
        logger.info("No obvious missing side-chain heavy atoms BEFORE.")

    sfxn = get_fa_scorefxn()

    # Task: repack only (no design) + dense Dunbrack rotamers
    tf = TaskFactory()
    tf.push_back(InitializeFromCommandline())  # will honor init flags
    tf.push_back(RestrictToRepacking())
    ex = ExtraRotamersGeneric()
    ex.ex1(True); ex.ex2(True); ex.ex2aro(True); ex.extrachi_cutoff(0)
    tf.push_back(ex)
    task = tf.create_task_and_apply_taskoperations(pose)

    logger.info("PackRotamersMover: rebuild/repack side chains")
    PackRotamersMover(sfxn, task).apply(pose)

    if do_relax:
        logger.info("FastRelax (coordinate constrained via init flags)")
        fr = FastRelax()
        fr.set_scorefxn(sfxn)
        fr.max_iter(200)
        fr.apply(pose)

    after = likely_missing_sc_heavy_atoms(pose)
    if after:
        logger.warning("Residues STILL look missing AFTER: %s", after)
        logger.warning("Rosetta cannot invent MISSING RESIDUES, only side-chain atoms on present residues.")
    else:
        logger.info("No obvious missing side-chain heavy atoms AFTER.")

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    pose.dump_pdb(str(output_pdb))
    logger.info("Saved → %s", output_pdb)

def build_parser():
    p = argparse.ArgumentParser(description="ACE-I prep: remove LPR/HOH, keep Zn, rebuild side chains (Dunbrack).")
    p.add_argument("--input", type=Path, required=True, help="Input PDB file")
    p.add_argument("--output", type=Path, required=True, help="Output PDB file")
    p.add_argument("--remove-zinc", action="store_true", help="Also remove Zn2+ (paper keeps it)")
    p.add_argument("--relax", action="store_true", help="Run a short constrained FastRelax")
    p.add_argument("--verbose", action="store_true", help="Debug logging")
    p.add_argument(
        "--init-flags",
        type=str,
        # Critical bits:
        #  -use_input_sc false → do NOT keep truncated side chains
        #  -ex1 -ex2aro -extrachi_cutoff 0 → dense rotamers
        #  -in:ignore_zero_occupancy false (and legacy form) → keep atoms with occ=0.00
        #  -out:pdb_write_hydrogens true → include polar H in output PDB
        #  Relax constraints via flags (portable across builds)
        default="-mute all "
                "-use_input_sc false -ex1 -ex2aro -extrachi_cutoff 0 "
                "-in:ignore_zero_occupancy false -ignore_zero_occupancy false "
                "-out:pdb_write_hydrogens true "
                "-relax:constrain_relax_to_start_coords -relax:ramp_constraints false "
                "-relax:coord_cst_stdev 0.5",
        help="Flags passed to PyRosetta init (quoted string)",
    )
    return p

def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    inp = args.input.expanduser().resolve()
    outp = args.output.expanduser().resolve()
    if not inp.exists():
        logger.error("Input not found: %s", inp); return 2

    logger.info("Initializing PyRosetta with flags: %s", args.init_flags)
    init(args.init_flags)

    try:
        prepare(inp, outp, keep_zinc=not args.remove_zinc, do_relax=args.relax)
        logger.info("Done.")
        return 0
    except Exception as e:
        logger.exception("Failed: %s", e)
        return 2

if __name__ == "__main__":
    raise SystemExit(main())
