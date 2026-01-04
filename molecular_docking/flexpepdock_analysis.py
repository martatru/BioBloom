#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FlexPepDock post-processing (PyMOL-driven, interface-restricted):
- H-bonds: chemical H-bonds (mode=2, donor/acceptor, after h_add) ONLY peptide (ligand) ↔ ACE (receptor)
  * Outputs: Number of H-bonds + list of distances (Å)
- Contacts: protein-only residues (polymer, no HOH) in an interface shell
- RMSD: PyMOL (super + rms_cur), heavy atoms, auto-pick best reference chain
- Clustering by peptide pose -> representative (centroid of largest cluster)

CSV per ligand + combined:
Ligand, Species, Number of H-bonds, H-bond distances (Å), Residues in contact, RMSD to native (1O86, heavy atoms)
"""

from __future__ import annotations
import argparse, csv, glob, math, os, re, sys, tempfile, subprocess, shutil
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

# ---------- deps ----------
try:
    import numpy as np
except ImportError:
    print("This script requires numpy. Install via: pip install numpy", file=sys.stderr)
    sys.exit(1)

# ---------- nice ACE residue names (extend if desired) ----------
ACE_KEY_NAMES = {
    281: "Gln281", 353: "His353", 354: "Ala354", 383: "His383", 384: "Glu384",
    387: "His387", 411: "Glu411", 511: "Lys511", 513: "His513",
    520: "Tyr520", 523: "Tyr523", 162: "Glu162",
}

# ---------- dataclasses ----------
@dataclass
class PoseInfo:
    name: str
    path: str
    bb_coords: Optional[np.ndarray] = None
    bb_index: Optional[List[Tuple[int, str]]] = None

@dataclass
class RowOut:
    ligand: str
    species: str
    hb_num: Optional[int]
    hb_dists_str: str
    key_contacts: str
    rmsd_vs_ref_heavy: Optional[float]

# ---------- PDB helpers ----------
def parse_chain_heavy_atom_map(pdb_path: str, chain_id: str) -> Dict[Tuple[int,str], np.ndarray]:
    amap: Dict[Tuple[int,str], np.ndarray] = {}
    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            ch = line[21:22]
            if ch != chain_id:
                continue
            aname = line[12:16].strip()
            if aname.startswith('H'):
                continue
            try:
                resi = int(line[22:26])
            except ValueError:
                continue
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            amap[(resi, aname)] = np.array([x,y,z], dtype=float)
    return amap

# ---------- RMSD/Kabsch (for clustering only) ----------
def kabsch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    C = np.dot(P.T, Q)
    V, S, Wt = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0
    if d:
        V[:, -1] = -V[:, -1]
    U = np.dot(V, Wt)
    return U

def rmsd_super(P: np.ndarray, Q: np.ndarray) -> float:
    if P.shape != Q.shape or P.size == 0:
        return float('inf')
    Pc = P - P.mean(axis=0); Qc = Q - Q.mean(axis=0)
    U = kabsch(Pc, Qc)
    P_rot = np.dot(Pc, U)
    diff = P_rot - Qc
    return math.sqrt((diff * diff).sum() / P.shape[0])

# ---------- clustering ----------
def cluster_by_threshold(pairwise: np.ndarray, cutoff: float) -> List[List[int]]:
    n = pairwise.shape[0]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if pairwise[i,j] <= cutoff:
                adj[i].append(j); adj[j].append(i)
    seen = [False]*n; comps = []
    for i in range(n):
        if seen[i]: continue
        stack = [i]; seen[i] = True; comp = []
        while stack:
            u = stack.pop(); comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True; stack.append(v)
        comps.append(sorted(comp))
    return comps

def pick_top_cluster(clusters: List[List[int]], poses: List[PoseInfo]) -> List[int]:
    return max(clusters, key=lambda comp: len(comp)) if clusters else []

def pick_centroid(indices: List[int], pairwise: np.ndarray) -> int:
    if not indices: return -1
    sub = pairwise[np.ix_(indices, indices)]
    return int(indices[int(np.argmin(sub.mean(axis=1)))])

# ---------- PyMOL helpers ----------
def _pymol_try_import():
    try:
        import pymol2
        return pymol2
    except Exception:
        return None

def _safe_chain_list(cmd, obj):
    try:
        return list(cmd.get_chains(obj))
    except Exception:
        return []

# ---------- PyMOL: H-bonds (interface-only, peptide↔ACE) ----------
def pymol_hbond_pairs(model_pdb: str,
                      receptor_chain: str,  # ACE chain in model
                      peptide_chain: str,   # peptide chain in model
                      cutoff_hb: float = 3.2, angle_hb: float = 55.0,
                      interface_shell: float = 4.0,
                      proxy_cutoff: float = 3.5) -> list[tuple[float, dict, dict]]:
    """
    Return list of (distance Å, info_rec, info_pep) for H-bond-like contacts between
    ACE (receptor_chain) and peptide (peptide_chain).

    1) Try chemical H-bonds (mode=2) after adding hydrogens, restricted to interface shell.
    2) If none found, fall back to proxy: all N/O…N/O pairs within proxy_cutoff Å (also interface-only).

    info_* = {'resi': int|str, 'resn': str, 'name': str, 'chain': str}
    """
    def _format_atom_info(cmd, sel):
        info = {'resi': None, 'resn': None, 'name': None, 'chain': None}
        def _cb(resi,resn,name,chain):
            info['resi']=int(resi) if str(resi).isdigit() else resi
            info['resn']=resn; info['name']=name; info['chain']=chain
        cmd.iterate(sel, "_cb(resi,resn,name,chain)", space={'_cb':_cb})
        return info

    def _pairs_to_list(cmd, pairs):
        out = []
        for (id_rec,_), (id_pep,_) in pairs:
            s1 = f"id {id_rec}"; s2 = f"id {id_pep}"
            try:
                d = float(cmd.get_distance(s1, s2))
            except Exception:
                cmd.distance("tmphb", s1, s2); d = float(cmd.get_distance(s1, s2)); cmd.delete("tmphb")
            out.append((d, _format_atom_info(cmd, s1), _format_atom_info(cmd, s2)))
        return out

    # --- Try in-memory (pymol2) ---
    pymol2 = _pymol_try_import()
    if pymol2 is not None:
        try:
            with pymol2.PyMOL() as pm:
                cmd = pm.cmd
                cmd.feedback('disable', 'all', 'everything')
                cmd.load(model_pdb, 'mob')

                rec_raw = f"(mob and chain {receptor_chain} and polymer and not resn HOH)"
                pep_raw = f"(mob and chain {peptide_chain} and polymer and not resn HOH)"

                rec_if = f"({rec_raw}) within {interface_shell} of ({pep_raw})"
                pep_if = f"({pep_raw}) within {interface_shell} of ({rec_raw})"

                # 1) Chemical H-bonds
                rec_da = f"({rec_if}) and (donor or acceptor)"
                pep_da = f"({pep_if}) and (donor or acceptor)"
                try:
                    cmd.h_add(rec_da); cmd.h_add(pep_da)
                except Exception:
                    pass
                pairs = cmd.find_pairs(rec_da, pep_da, mode=2, cutoff=cutoff_hb, angle=angle_hb)
                out = _pairs_to_list(cmd, pairs)

                # 2) Proxy fallback if needed
                if not out:
                    rec_no = f"({rec_if}) and name N+O"
                    pep_no = f"({pep_if}) and name N+O"
                    pairs2 = cmd.find_pairs(rec_no, pep_no, mode=0, cutoff=proxy_cutoff)
                    out = _pairs_to_list(cmd, pairs2)

                # sort by distance
                out.sort(key=lambda x: x[0])
                return out
        except Exception:
            pass

    # --- Fallback: external pymol (ensure doubled braces for .format) ---
    import tempfile, subprocess, re, os, shutil
    pymol_bin = shutil.which('pymol')
    if pymol_bin is None:
        return []
    script = """from pymol import cmd
cmd.feedback('disable','all','everything')
cmd.load(r'''{model_pdb}''','mob')
rec_raw = "mob and chain {rc} and polymer and not resn HOH"
pep_raw = "mob and chain {pc} and polymer and not resn HOH"
rec_if = "({{}}) within {shell} of ({{}})".format(rec_raw, pep_raw)
pep_if = "({{}}) within {shell} of ({{}})".format(pep_raw, rec_raw)

def pr_atom(id_):
    out = []
    def _cb(resi,resn,name,chain):
        try: ri = int(resi)
        except Exception: ri = resi
        out.append((ri,resn,name,chain))
    cmd.iterate("id {{}}".format(id_), "_cb(resi,resn,name,chain)", space={{'_cb':_cb}})
    return out[0] if out else (None,None,None,None)

# 1) chemical (mode=2)
rec_da = "({{}}) and (donor or acceptor)".format(rec_if)
pep_da = "({{}}) and (donor or acceptor)".format(pep_if)
try:
    cmd.h_add(rec_da); cmd.h_add(pep_da)
except Exception:
    pass
pairs = cmd.find_pairs(rec_da, pep_da, mode=2, cutoff={cut_hb}, angle={ang})
if len(pairs)==0:
    # 2) proxy N/O…N/O
    rec_no = "({{}}) and name N+O".format(rec_if)
    pep_no = "({{}}) and name N+O".format(pep_if)
    pairs = cmd.find_pairs(rec_no, pep_no, mode=0, cutoff={proxy})

for (i1,_), (i2,_) in pairs:
    d = cmd.get_distance("id {{}}".format(i1), "id {{}}".format(i2))
    r1 = pr_atom(i1); r2 = pr_atom(i2)
    print("HB\\t{{:.3f}}\\t{{}}\\t{{}}\\t{{}}\\t{{}}\\t{{}}\\t{{}}\\t{{}}".format(d, r1[0], r1[1], r1[2], r1[3], r2[0], r2[1], r2[2], r2[3]))
""".format(
        model_pdb=model_pdb.replace("\\","\\\\"),
        rc=receptor_chain, pc=peptide_chain,
        shell=interface_shell, cut_hb=cutoff_hb, ang=angle_hb, proxy=proxy_cutoff
    )
    with tempfile.NamedTemporaryFile('w', suffix='.py', delete=False) as tf:
        tf.write(script); tpath = tf.name
    try:
        proc = subprocess.run([pymol_bin, '-cq', tpath], capture_output=True, text=True, timeout=300)
        out_text = (proc.stdout or "") + (proc.stderr or "")
        out = []
        for m in re.finditer(r"^HB\t([0-9.]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)", out_text, flags=re.M):
            dist = float(m.group(1))
            rec = {'resi': int(m.group(2)) if m.group(2).isdigit() else m.group(2),
                   'resn': m.group(3), 'name': m.group(4), 'chain': m.group(5)}
            pep = {'resi': int(m.group(6)) if m.group(6).isdigit() else m.group(6),
                   'resn': m.group(7), 'name': m.group(8), 'chain': None}
            out.append((dist, rec, pep))
        out.sort(key=lambda x: x[0])
        return out
    finally:
        try: os.remove(tpath)
        except Exception: pass

# ---------- PyMOL: RMSD (heavy) with auto ref chain ----------
def pymol_rmsd_heavy(model_pdb: str, chain_model: str, ref_pdb: str, chain_ref: str | None) -> tuple[Optional[float], Optional[str], int, int]:
    def _super_and_rms(cmd, sel_mob, sel_ref):
        n_mob = cmd.count_atoms(sel_mob); n_ref = cmd.count_atoms(sel_ref)
        if n_mob < 4 or n_ref < 4:
            return None, n_mob, n_ref
        try:
            res = cmd.super(sel_mob, sel_ref, cycles=0)
            if isinstance(res,(list,tuple)) and len(res)>0:
                rms_val = float(res[0])
            else:
                rms_val = cmd.rms_cur(sel_mob, sel_ref, matchmaker=4)
            if rms_val is None or math.isnan(rms_val) or math.isinf(rms_val):
                return None, n_mob, n_ref
            return float(rms_val), n_mob, n_ref
        except Exception:
            return None, n_mob, n_ref

    pymol2 = _pymol_try_import()
    if pymol2 is not None:
        try:
            with pymol2.PyMOL() as pm:
                cmd = pm.cmd
                cmd.feedback('disable','all','everything')
                cmd.load(ref_pdb, 'ref'); cmd.load(model_pdb, 'mob')
                sel_mob = f"(mob and chain {chain_model} and not elem H)"
                ref_chains = _safe_chain_list(cmd, 'ref')
                candidates = ([chain_ref] if chain_ref else []) + [c for c in ref_chains if c != chain_ref]
                best = (None, None, 0, 0)
                for c in (candidates if candidates else ref_chains):
                    sel_ref = f"(ref and chain {c} and not elem H)"
                    rms, nm, nr = _super_and_rms(cmd, sel_mob, sel_ref)
                    if rms is not None and (best[0] is None or rms < best[0]):
                        best = (rms, c, nm, nr)
                return best
        except Exception:
            pass

    pymol_bin = shutil.which('pymol')
    if pymol_bin is None:
        return None, None, 0, 0

    script = """from pymol import cmd
import math
cmd.feedback('disable','all','everything')
cmd.load(r'''{ref_pdb}''','ref'); cmd.load(r'''{model_pdb}''','mob')
ref_chains = cmd.get_chains('ref')
candidates = []
pref = {pref}
if pref is not None and pref in ref_chains:
    candidates.append(pref)
for c in ref_chains:
    if c != pref:
        candidates.append(c)
sel_mob = "mob and chain {cm} and not elem H"
best = (None, None, 0, 0)
for c in (candidates if candidates else ref_chains):
    sel_ref = "ref and chain {{}} and not elem H".format(c)
    n_mob = cmd.count_atoms(sel_mob)
    n_ref = cmd.count_atoms(sel_ref)
    rms = None
    if n_mob >= 4 and n_ref >= 4:
        try:
            res = cmd.super(sel_mob, sel_ref, cycles=0)
            if isinstance(res,(list,tuple)) and len(res)>0:
                rms_val = float(res[0])
            else:
                rms_val = cmd.rms_cur(sel_mob, sel_ref, matchmaker=4)
            if rms_val is not None and not math.isnan(rms_val) and not math.isinf(rms_val):
                rms = float(rms_val)
        except Exception:
            rms = None
    print("BEST_TRY c={} mob={} ref={} rms={}".format(c, n_mob, n_ref, rms))
    if rms is not None and (best[0] is None or rms < best[0]):
        best = (rms, c, n_mob, n_ref)
print("BEST c={} rms={} mob={} ref={}".format(best[1], best[0], best[2], best[3]))
""".format(ref_pdb=ref_pdb.replace("\\","\\\\"), model_pdb=model_pdb.replace("\\","\\\\"),
           pref=repr(chain_ref) if chain_ref else "None", cm=chain_model)
    with tempfile.NamedTemporaryFile('w', suffix='.py', delete=False) as tf:
        tf.write(script); tpath = tf.name
    try:
        proc = subprocess.run([pymol_bin, '-cq', tpath], capture_output=True, text=True, timeout=300)
        out = (proc.stdout or "") + (proc.stderr or "")
        m = re.search(r"BEST c=([A-Za-z0-9]) rms=([0-9.+-eE]+|None) mob=([0-9]+) ref=([0-9]+)", out)
        if not m:
            return None, None, 0, 0
        c = m.group(1); rms_s = m.group(2); n_mob = int(m.group(3)); n_ref = int(m.group(4))
        rms = None if rms_s == "None" else float(rms_s)
        return rms, c, n_mob, n_ref
    finally:
        try: os.remove(tpath)
        except Exception: pass

# ---------- PyMOL: contacts (protein-only, interface) ----------
def pymol_contacts(model_pdb: str, receptor_chain: str, peptide_chain: str, cutoff: float = 4.0) -> tuple[list[str], str]:
    """
    Find all receptor residues (polymer only) within cutoff Å of the peptide.
    Output full residue names like 'Glu384', 'Tyr523', etc.
    """
    def _sort_key(name):
        import re
        m = re.search(r'(\d+)', name)
        return int(m.group(1)) if m else 10**9

    pymol2 = _pymol_try_import()
    if pymol2 is not None:
        try:
            with pymol2.PyMOL() as pm:
                cmd = pm.cmd
                cmd.feedback('disable', 'all', 'everything')
                cmd.load(model_pdb, 'mob')

                rec_sel = f"(mob and chain {receptor_chain} and polymer and not resn HOH)"
                pep_sel = f"(mob and chain {peptide_chain} and polymer and not resn HOH)"
                cmd.select('contact_resi', f"({rec_sel}) within {cutoff} of ({pep_sel})")

                names = set()
                def _cb(resn, resi):
                    # always print proper residue name
                    names.add(f"{resn}{resi}")
                cmd.iterate("contact_resi and name CA", "_cb(resn, resi)", space={'_cb': _cb})

                lst = sorted(names, key=_sort_key)
                return lst, (", ".join(lst) if lst else "—")
        except Exception as e:
            print(f"[WARN] pymol2 contacts failed: {e}")
            pass

    # fallback: external pymol process
    pymol_bin = shutil.which('pymol')
    if pymol_bin is None:
        return [], "—"

    script = """from pymol import cmd
cmd.feedback('disable','all','everything')
cmd.load(r'''{model_pdb}''','mob')
rec_sel = "mob and chain {rc} and polymer and not resn HOH"
pep_sel = "mob and chain {pc} and polymer and not resn HOH"
cmd.select('contact_resi', "({{}}) within {cutoff} of ({{}})".format(rec_sel, pep_sel))
def _cb(resn, resi):
    print("CONTACT {}{}".format(resn, resi))
cmd.iterate("contact_resi and name CA", "_cb(resn, resi)", space={{'_cb':_cb}})
""".format(model_pdb=model_pdb.replace("\\","\\\\"), rc=receptor_chain, pc=peptide_chain, cutoff=cutoff)

    import tempfile, subprocess, re, os
    with tempfile.NamedTemporaryFile('w', suffix='.py', delete=False) as tf:
        tf.write(script); tpath = tf.name
    try:
        proc = subprocess.run([pymol_bin, '-cq', tpath], capture_output=True, text=True, timeout=300)
        out = (proc.stdout or "") + (proc.stderr or "")
        names = []
        for match in re.findall(r"CONTACT\s+([A-Za-z0-9]+)", out):
            names.append(match)
        names_sorted = sorted(set(names), key=_sort_key)
        return names_sorted, (", ".join(names_sorted) if names_sorted else "—")
    finally:
        try: os.remove(tpath)
        except Exception: pass

# ---------- format helpers ----------
def fmt(x: Optional[float]) -> str:
    if x is None or (isinstance(x,float) and (math.isnan(x) or math.isinf(x))):
        return ""
    return f"{x:.2f}"

def fmt_int(x: Optional[int]) -> str:
    if x is None:
        return ""
    return str(int(x))

# ---------- MAIN ----------
def main():
    ap = argparse.ArgumentParser(description="Interface H-bonds (count+distances), contacts, RMSD (heavy) via PyMOL.")
    ap.add_argument('--root', required=True, help='Parent folder containing per-ligand folders')
    ap.add_argument('--glob', default='*/output', help='Subfolder glob (default: */output)')
    ap.add_argument('--receptor-chain', default='A', help='ACE chain id in models')
    ap.add_argument('--peptide-chain', default='P', help='Peptide ligand chain id in models')
    ap.add_argument('--rmsd-cutoff', type=float, default=2.0, help='Peptide pose clustering cutoff (Å)')
    ap.add_argument('--contact-cutoff', type=float, default=4.0, help='Interface shell for contacts/H-bonds (Å)')
    ap.add_argument('--ref-pdb', default='/home/marta/Pulpit/proteomic_data/pdb/1O86.pdb')
    ap.add_argument('--ref-peptide-chain', default=None, help='If None, auto-pick reference chain with best RMSD')
    ap.add_argument('--species-map', default=None, help='CSV two columns: ligand_id,species')
    ap.add_argument('--outdir', default='/home/marta/Pulpit/docking_results_analysis')
    ap.add_argument('--out', default='results_min.csv')
    ap.add_argument('--debug', action='store_true')
    args = ap.parse_args()

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)

    ligand_species: Dict[str,str] = {}
    if args.species_map and os.path.isfile(args.species_map):
        with open(args.species_map,'r') as fh:
            rdr = csv.reader(fh)
            for row in rdr:
                if len(row) >= 2:
                    ligand_species[row[0].strip()] = row[1].strip()

    rows: List[RowOut] = []

    folders = sorted(glob.glob(os.path.join(args.root, args.glob)))
    if not folders:
        print("No folders matched. Check --root and --glob.", file=sys.stderr)
        sys.exit(2)

    for fdir in folders:
        parent = os.path.basename(os.path.dirname(os.path.normpath(fdir)))
        ligand_label = parent
        species = ligand_species.get(ligand_label, "")

        model_paths = sorted(glob.glob(os.path.join(fdir, 'model_*.pdb')))
        if not model_paths:
            model_paths = sorted(glob.glob(os.path.join(fdir, 'complex.ppk_*.pdb')))
        if not model_paths:
            print(f"[WARN] No models in {fdir}")
            continue

        # build peptide arrays for clustering (heavy)
        poses: List[PoseInfo] = []
        for mp in model_paths:
            m = PoseInfo(name=os.path.basename(mp), path=mp)
            amap = parse_chain_heavy_atom_map(mp, args.peptide_chain)
            keys = sorted(amap.keys())
            m.bb_index = keys
            m.bb_coords = np.vstack([amap[k] for k in keys]) if keys else np.zeros((0,3), float)
            poses.append(m)

        # choose template
        template_keys = None
        for pz in poses:
            if pz.bb_coords is not None and pz.bb_coords.size > 0:
                template_keys = pz.bb_index; break
        if template_keys is None:
            print(f"[WARN] No peptide atoms found in {fdir} (chain {args.peptide_chain}). Skipping.")
            continue

        def extract_matching(p: PoseInfo) -> Optional[np.ndarray]:
            if not p.bb_index:
                return None
            idx_map = {ra:i for i,ra in enumerate(p.bb_index)}
            rowsB = []
            for ra in template_keys:
                i = idx_map.get(ra)
                if i is None: return None
                rowsB.append(p.bb_coords[i])
            return np.vstack(rowsB)

        aligned = [extract_matching(p) for p in poses]
        valid_idx = [i for i,a in enumerate(aligned) if a is not None and a.size>0]
        if len(valid_idx) < 2:
            print(f"[WARN] Not enough comparable models in {fdir} for clustering.")
            continue

        n = len(poses)
        pair = np.full((n,n), np.inf, dtype=float)
        for i in valid_idx: pair[i,i] = 0.0
        for ii in range(len(valid_idx)):
            for jj in range(ii+1, len(valid_idx)):
                a = valid_idx[ii]; b = valid_idx[jj]
                pair[a,b] = pair[b,a] = rmsd_super(aligned[a], aligned[b])

        comps = cluster_by_threshold(pair, args.rmsd_cutoff)
        if not comps:
            print(f"[WARN] No clusters in {fdir}")
            continue
        top = pick_top_cluster(comps, poses)
        centroid_idx = pick_centroid(top, pair)
        if centroid_idx < 0: centroid_idx = top[0]
        centroid = poses[centroid_idx]

        # ---- PyMOL calculations for representative ----
        # H-bond pairs (peptide↔ACE only), distances list
        hb_pairs = pymol_hbond_pairs(
            centroid.path, receptor_chain=args.receptor_chain, peptide_chain=args.peptide_chain,
            cutoff_hb=3.2, angle_hb=55.0, interface_shell=args.contact_cutoff
        )
        # Sort by distance and produce printable list (2 d.p.)
        hb_dists = sorted([p[0] for p in hb_pairs])
        hb_dists_str = ", ".join(f"{d:.2f}" for d in hb_dists)
        hb_num = len(hb_pairs)

        # Contacts
        _, key_contacts = pymol_contacts(centroid.path, args.receptor_chain, args.peptide_chain, args.contact_cutoff)

        # RMSD heavy vs reference
        rmsd_ref, ref_chain_used, n_mob, n_ref = (None, None, 0, 0)
        if args.ref_pdb:
            rmsd_ref, ref_chain_used, n_mob, n_ref = pymol_rmsd_heavy(
                centroid.path, args.peptide_chain, args.ref_pdb, args.ref_peptide_chain
            )

        rows.append(RowOut(
            ligand=ligand_label, species=species,
            hb_num=hb_num, hb_dists_str=hb_dists_str,
            key_contacts=key_contacts, rmsd_vs_ref_heavy=rmsd_ref
        ))

        # per-ligand CSV
        single_out = os.path.join(args.outdir, f"{ligand_label}.csv")
        with open(single_out, 'w', newline='') as sfh:
            w1 = csv.writer(sfh)
            w1.writerow(["Ligand","Species","Number of H-bonds","H-bond distances (Å)","Residues in contact","RMSD to native (1O86, heavy atoms)"])
            w1.writerow([ligand_label, species, hb_num, hb_dists_str, key_contacts, fmt(rmsd_ref)])

        diag = f" (ref_chain={ref_chain_used}, mob_atoms={n_mob}, ref_atoms={n_ref})" if args.debug else ""
        print(f"{ligand_label}: HBonds={hb_num} [{hb_dists_str}] ; Contacts={key_contacts}; RMSD_heavy_vs_1O86={fmt(rmsd_ref)}{diag}")

    # combined CSV
    with open(args.out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(["Ligand","Species","Number of H-bonds","H-bond distances (Å)","Residues in contact","RMSD to native (1O86, heavy atoms)"])
        for r in rows:
            w.writerow([r.ligand, r.species, r.hb_num, r.hb_dists_str, r.key_contacts, fmt(r.rmsd_vs_ref_heavy)])
    print(f"Wrote {args.out} with {len(rows)} rows.")

if __name__ == '__main__':
    main()
