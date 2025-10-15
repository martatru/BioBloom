#!/usr/bin/env python3
"""
FlexPepDock post-processing for poster-ready tables.

Now includes:
- CLUSTER ALL MODELS by peptide pose (RMSD), not just those with metrics
- Reference validation:
    * --ref-pdb with --rmsd-mode {backbone,heavy} -> RMSD of ligand vs reference
    * optional re-scoring of reference (InterfaceAnalyzer / score_jd2)
- Representative selection: centroid of top cluster; if centroid lacks metrics,
  we try to swap to best-scored member with metrics (but clustering ALWAYS uses all models)

Outputs (per ligand folder + global summary):
Ligand, Gatunek, Energia interfejsu (REU), Energia całkowita (REU), Wiązania wodorowe,
Reszty w kontakcie, Potencjał inhibitorowy,
RMSD vs ref (Å), I_sc ref (REU), ΔI_sc (REU), Total ref (REU), ΔTotal (REU)

Usage example at bottom.
"""
from __future__ import annotations
import argparse, csv, glob, math, os, re, subprocess, sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

try:
    import numpy as np
except ImportError:
    print("This script requires numpy. Install via: pip install numpy", file=sys.stderr)
    sys.exit(1)

# --- ACE residues of interest (edit if needed) ---
ACE_KEY_RESIS = {281, 353, 354, 383, 384, 387, 411, 511, 513, 520, 523, 162}
ACE_KEY_NAMES = {
    281: "Gln281", 353: "His353", 354: "Ala354", 383: "His383", 384: "Glu384",
    387: "His387", 411: "Glu411", 511: "Lys511", 513: "His513",
    520: "Tyr520", 523: "Tyr523", 162: "Glu162",
}

@dataclass
class PoseInfo:
    name: str
    path: str
    score_total: Optional[float] = None
    score_reweighted: Optional[float] = None
    score_I_sc: Optional[float] = None
    bb_coords: Optional[np.ndarray] = None
    bb_index: Optional[List[Tuple[int, str]]] = None

@dataclass
class RowOut:
    ligand: str
    species: str
    iface_energy_reu: Optional[float]
    total_energy_reu: Optional[float]
    hbonds: Optional[int]
    key_contacts: str
    inhibitor_potential: str
    rmsd_vs_ref: Optional[float]
    i_sc_ref: Optional[float]
    di_sc: Optional[float]
    total_ref: Optional[float]
    dtotal: Optional[float]

# -------------------- PDB helpers --------------------
def parse_pdb_backbone_chain(pdb_path: str, chain_id: str, atoms=("N","CA","C","O")) -> Tuple[np.ndarray, List[Tuple[int,str]]]:
    residues: Dict[int, Dict[str, np.ndarray]] = defaultdict(dict)
    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            aname = line[12:16].strip()
            ch = line[21:22]
            try:
                resi = int(line[22:26])
            except ValueError:
                continue
            if ch != chain_id:
                continue
            if aname in atoms:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                residues[resi][aname] = np.array([x,y,z], dtype=float)
    ordered, index = [], []
    for resi in sorted(residues.keys()):
        atoms_here = residues[resi]
        if all(a in atoms_here for a in atoms):
            for a in atoms:
                ordered.append(atoms_here[a])
                index.append((resi, a))
    if not ordered:
        return np.zeros((0,3), float), []
    return np.vstack(ordered), index

def parse_chain_heavy_atom_map(pdb_path: str, chain_id: str) -> Dict[Tuple[int,str], np.ndarray]:
    """Return {(resi, atomname)->xyz} for heavy atoms of given chain, sorted by (resi, atomname)."""
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

def parse_heavy_atoms_list(pdb_path: str, chain_id: str) -> List[np.ndarray]:
    coords = []
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
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            coords.append(np.array([x,y,z], dtype=float))
    return coords

# -------------------- RMSD / Kabsch --------------------
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
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    U = kabsch(Pc, Qc)
    P_rot = np.dot(Pc, U)
    diff = P_rot - Qc
    return math.sqrt((diff * diff).sum() / P.shape[0])

# -------------------- Clustering --------------------
def cluster_by_threshold(pairwise: np.ndarray, cutoff: float) -> List[List[int]]:
    n = pairwise.shape[0]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if pairwise[i,j] <= cutoff:
                adj[i].append(j); adj[j].append(i)
    seen = [False]*n
    comps = []
    for i in range(n):
        if seen[i]:
            continue
        stack = [i]; seen[i] = True
        comp = []
        while stack:
            u = stack.pop(); comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True; stack.append(v)
        comps.append(sorted(comp))
    return comps

# -------------------- Optional Rosetta runners --------------------
def run_interface_analyzer(exe: str, pdb: str, receptor_chain: str, peptide_chain: str) -> Tuple[Optional[float], Optional[int]]:
    try:
        cmd = [exe, '-s', pdb, '-interface', f'{receptor_chain}_{peptide_chain}',
               '-pack_input', 'false', '-pack_separated', 'false',
               '-out:file:scorefile', 'ia_tmp.sc']
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        iface = None; hb = None; header = None
        with open('ia_tmp.sc', 'r') as fh:
            for line in fh:
                if line.startswith('SCORE:'):
                    parts = line.strip().split()
                    if header is None:
                        header = parts; continue
                    try: idx_dg = header.index('dG_separated')
                    except ValueError: idx_dg = None
                    try: idx_hb = header.index('hbonds_int')
                    except ValueError: idx_hb = None
                    if idx_dg is not None: iface = float(parts[idx_dg])
                    if idx_hb is not None: hb = int(float(parts[idx_hb]))
        try: os.remove('ia_tmp.sc')
        except OSError: pass
        return iface, hb
    except Exception:
        return None, None

def run_score_jd2(exe: str, pdb: str) -> Optional[float]:
    try:
        cmd = [exe, '-s', pdb, '-out:file:scorefile', 'total_tmp.sc']
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        total = None; header = None
        with open('total_tmp.sc', 'r') as fh:
            for line in fh:
                if line.startswith('SCORE:'):
                    parts = line.strip().split()
                    if header is None:
                        header = parts; continue
                    try:
                        idx = header.index('total_score'); total = float(parts[idx])
                    except ValueError:
                        pass
        try: os.remove('total_tmp.sc')
        except OSError: pass
        return total
    except Exception:
        return None

# -------------------- Tables & selection --------------------
def read_metrics_table(path: str) -> Dict[str, Dict[str, float]]:
    """Read CSV/TSV and return dict: model_name -> column dict.
       Adds aliases: 'index' → 'model_%02d.pdb' for ROSIE top10 tables.
    """
    if not os.path.isfile(path):
        return {}
    sep = ',' if path.endswith('.csv') else '\t'
    rows: Dict[str, Dict[str,float]] = {}
    with open(path, 'r') as fh:
        header = None
        row_counter = 0
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(sep)
            if header is None:
                header = parts
                continue
            d_raw = {header[i]: parts[i] for i in range(min(len(header), len(parts)))}
            # cast numerics
            d = {}
            for k, v in d_raw.items():
                try: d[k] = float(v)
                except Exception: d[k] = v
            # store under raw key if present
            key = (d_raw.get('description') or d_raw.get('model') or
                   d_raw.get('decoy') or d_raw.get('name') or d_raw.get('Model'))
            if key:
                rows[str(key)] = d
            # alias from 'index' -> model_%02d.pdb
            idx_val = None
            if 'index' in d_raw:
                try: idx_val = int(float(d_raw['index']))
                except Exception: idx_val = None
            else:
                idx_val = row_counter  # fallback sequential
            if idx_val is not None:
                alias = f"model_{idx_val+1:02d}.pdb"
                rows[alias] = d
            # also store any value that looks like a pdb/model id
            for v in d_raw.values():
                if isinstance(v, str) and (v.endswith('.pdb') or v.startswith('model_') or v.startswith('complex.ppk_')):
                    rows[str(v)] = d
            row_counter += 1
    return rows

def pick_top_cluster(clusters: List[List[int]], poses: List[PoseInfo]) -> List[int]:
    best_idx = None; best_key = (float('inf'), -len(poses))
    for i, comp in enumerate(clusters):
        scores = []
        for idx in comp:
            sc = poses[idx].score_reweighted if poses[idx].score_reweighted is not None else poses[idx].score_total
            if sc is not None: scores.append(sc)
        med = np.median(scores) if scores else float('inf')
        key = (med, -len(comp))
        if key < best_key:
            best_key = key; best_idx = i
    return clusters[best_idx] if best_idx is not None else (clusters[0] if clusters else [])

def pick_centroid(indices: List[int], pairwise: np.ndarray) -> int:
    if not indices: return -1
    sub = pairwise[np.ix_(indices, indices)]
    return int(indices[int(np.argmin(sub.mean(axis=1)))])

# -------------------- Contacts & H-bonds --------------------
def format_key_contacts(pdb_path: str, receptor_chain: str, peptide_chain: str, cutoff: float=4.0) -> str:
    rec_atoms = parse_residue_heavy_atoms(pdb_path, receptor_chain)
    pep_heavy = parse_heavy_atoms_list(pdb_path, peptide_chain)
    if not pep_heavy or not rec_atoms:
        return "—"
    pep = np.vstack(pep_heavy)
    hits: List[int] = []
    for resi, atoms in rec_atoms.items():
        A = np.vstack(atoms)
        dmin = np.linalg.norm(pep[:,None,:] - A[None,:,:], axis=2).min()
        if dmin <= cutoff:
            hits.append(resi)
    key = [ACE_KEY_NAMES[r] for r in hits if r in ACE_KEY_RESIS]
    if key:
        return ", ".join(sorted(key, key=lambda x: int(re.findall(r"(\d+)", x)[0])))
    # fallback: nearest 3
    dist_map = []
    for resi, atoms in rec_atoms.items():
        A = np.vstack(atoms)
        dmin = np.linalg.norm(pep[:,None,:] - A[None,:,:], axis=2).min()
        dist_map.append((resi, dmin))
    dist_map.sort(key=lambda x: x[1])
    return ", ".join(ACE_KEY_NAMES.get(r, f"Res{r}") for r,_ in dist_map[:3]) if dist_map else "—"

def parse_residue_heavy_atoms(pdb_path: str, chain_id: str) -> Dict[int, List[np.ndarray]]:
    res_atoms: Dict[int, List[np.ndarray]] = defaultdict(list)
    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not line.startswith("ATOM"):
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
            res_atoms[resi].append(np.array([x,y,z], dtype=float))
    return res_atoms

def count_hbonds_simple(pdb_path: str, receptor_chain: str, peptide_chain: str, dist_cut: float=3.5) -> Optional[int]:
    """Proxy without H atoms: count N/O…N/O pairs across interface ≤ dist_cut (Å)."""
    rec, pep = [], []
    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not (line.startswith('ATOM') or line.startswith('HETATM')): continue
            ch = line[21:22]; aname = line[12:16].strip()
            if aname.startswith('H'): continue
            if aname[0] not in ('N','O'): continue
            try: x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError: continue
            if ch == receptor_chain: rec.append(np.array([x,y,z], float))
            elif ch == peptide_chain: pep.append(np.array([x,y,z], float))
    if not rec or not pep: return None
    R = np.vstack(rec); P = np.vstack(pep)
    D = np.linalg.norm(P[:,None,:] - R[None,:,:], axis=2)
    return int((D <= dist_cut).sum())

def decide_potential(iface: Optional[float], hbonds: Optional[int], key_contacts: str) -> str:
    if iface is None:
        return "(brak danych)"
    key_ok = any(k in key_contacts for k in ["Glu384", "Tyr523", "His383", "His513"])
    hb_ok = (hbonds is not None and hbonds >= 4)
    return "wysoki" if (iface <= -10.0 and hb_ok and key_ok) else "niski"

# -------------------- RMSD vs reference --------------------
def build_atom_array_for_rmsd(pdb_path: str, chain_id: str, mode: str) -> Tuple[np.ndarray, List[Tuple[int,str]]]:
    if mode == 'backbone':
        arr, idx = parse_pdb_backbone_chain(pdb_path, chain_id)
        return arr, idx
    # heavy: build sorted map -> array
    amap = parse_chain_heavy_atom_map(pdb_path, chain_id)
    if not amap:
        return np.zeros((0,3), float), []
    keys = sorted(amap.keys())
    arr = np.vstack([amap[k] for k in keys])
    return arr, keys

def rmsd_between_chains(pdb_model: str, chain_model: str, pdb_ref: str, chain_ref: str, mode: str='heavy') -> Optional[float]:
    """Compute RMSD between ligand atoms in model and reference.
       For 'heavy': intersect common (resi,atom) keys; for 'backbone': backbone sequence alignment via indices.
    """
    if mode == 'backbone':
        A, idxA = parse_pdb_backbone_chain(pdb_model, chain_model)
        B, idxB = parse_pdb_backbone_chain(pdb_ref, chain_ref)
        if A.size==0 or B.size==0:
            return None
        # intersect by (resi,atom)
        setA = {k:i for i,k in enumerate(idxA)}
        selA, selB = [], []
        for j,k in enumerate(idxB):
            if k in setA:
                selA.append(setA[k]); selB.append(j)
        if not selA:
            return None
        return rmsd_super(A[np.array(selA)], B[np.array(selB)])
    else:
        amapA = parse_chain_heavy_atom_map(pdb_model, chain_model)
        amapB = parse_chain_heavy_atom_map(pdb_ref, chain_ref)
        common = sorted(set(amapA.keys()).intersection(amapB.keys()))
        if len(common) < 4:
            return None
        A = np.vstack([amapA[k] for k in common])
        B = np.vstack([amapB[k] for k in common])
        return rmsd_super(A, B)

# -------------------- Main --------------------
def main():
    ap = argparse.ArgumentParser(description="Cluster FlexPepDock outputs and export poster-ready table.")
    ap.add_argument('--root', required=True, help='Parent directory that contains per-peptide folders')
    ap.add_argument('--glob', default='*/output', help='Subfolder glob (default: */output)')
    ap.add_argument('--receptor-chain', default='A')
    ap.add_argument('--peptide-chain', default='P')
    ap.add_argument('--rmsd-cutoff', type=float, default=2.0)
    ap.add_argument('--contact-cutoff', type=float, default=4.0)
    ap.add_argument('--rmsd-mode', choices=['backbone','heavy'], default='backbone', help='Atoms used for clustering and (if no ref) internal RMSDs')
    # Reference validation
    ap.add_argument('--ref-pdb', default=None, help='Path to reference PDB with known ligand pose')
    ap.add_argument('--ref-receptor-chain', default=None, help='Override receptor chain in ref (default: same as --receptor-chain)')
    ap.add_argument('--ref-peptide-chain', default=None, help='Override ligand chain in ref (default: same as --peptide-chain)')
    ap.add_argument('--score-ref', action='store_true', help='If set and IA/score_jd2 are available, compute energies for reference structure')
    # Optional Rosetta binaries
    ap.add_argument('--interfaceanalyzer', default=None)
    ap.add_argument('--score_jd2', default=None)
    # Metadata & output
    ap.add_argument('--species-map', default=None, help='CSV two columns: ligand_id,species list')
    ap.add_argument('--outdir', default='/home/marta/Pulpit/docking_results_analysis')
    ap.add_argument('--out', default='results.csv')
    ap.add_argument('--debug', action='store_true')
    args = ap.parse_args()

    ref_rec_chain = args.ref_receptor_chain or args.receptor_chain
    ref_pep_chain = args.ref_peptide_chain or args.peptide_chain

    ligand_species: Dict[str,str] = {}
    if args.species_map and os.path.isfile(args.species_map):
        with open(args.species_map,'r') as fh:
            rdr = csv.reader(fh)
            for row in rdr:
                if len(row) >= 2:
                    ligand_species[row[0].strip()] = row[1].strip()

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)

    rows: List[RowOut] = []

    folders = sorted(glob.glob(os.path.join(args.root, args.glob)))
    if not folders:
        print("No folders matched. Check --root and --glob.", file=sys.stderr)
        sys.exit(2)

    # Precompute reference energies if requested
    ref_i_sc = None; ref_total = None
    if args.ref_pdb and args.score_ref:
        if args.interfaceanalyzer:
            ref_i_sc, _ = run_interface_analyzer(args.interfaceanalyzer, args.ref_pdb, ref_rec_chain, ref_pep_chain)
        if args.score_jd2:
            ref_total = run_score_jd2(args.score_jd2, args.ref_pdb)

    for fdir in folders:
        # Ligand = parent folder (one level above 'output')
        parent = os.path.basename(os.path.dirname(os.path.normpath(fdir)))
        ligand_label = parent
        species = ligand_species.get(ligand_label, "")

        # Models (ALL models; cluster all)
        model_paths = sorted(glob.glob(os.path.join(fdir, 'model_*.pdb')))
        if not model_paths:
            model_paths = sorted(glob.glob(os.path.join(fdir, 'complex.ppk_*.pdb')))
        if not model_paths:
            print(f"[WARN] No models in {fdir}")
            continue

        # Load metrics tables for mapping I_sc/total_score to models (if present)
        metrics = {}
        for cand in ('top10_scores.tsv','filtered_scores.tsv','metrics.csv'):
            p = os.path.join(fdir, cand)
            if os.path.isfile(p):
                metrics.update(read_metrics_table(p))

        # Build pose objects (scores if available)
        poses: List[PoseInfo] = []
        for mp in model_paths:
            name = os.path.basename(mp)
            m = PoseInfo(name=name, path=mp)
            md = metrics.get(name) or metrics.get(name.replace('.pdb','')) or {}
            for key,val in md.items():
                if not isinstance(val,(int,float)): continue
                low = key.lower()
                if 'reweighted' in low: m.score_reweighted = float(val)
                if ('total_score' in low or low=='score'): m.score_total = float(val)
                if low in ('i_sc','isc','interface_delta','dg_separated'): m.score_I_sc = float(val)
            # atoms for clustering
            if args.rmsd_mode == 'backbone':
                bb, idx = parse_pdb_backbone_chain(mp, args.peptide_chain)
                m.bb_coords = bb; m.bb_index = idx
            else:  # heavy
                amap = parse_chain_heavy_atom_map(mp, args.peptide_chain)
                keys = sorted(amap.keys())
                m.bb_index = keys
                m.bb_coords = np.vstack([amap[k] for k in keys]) if keys else np.zeros((0,3), float)
            poses.append(m)

        # comparable set: create a common key list for intersection (by bb_index)
        # choose first non-empty as template
        template_keys = None
        for pz in poses:
            if pz.bb_coords is not None and pz.bb_coords.size>0:
                template_keys = pz.bb_index; break
        if template_keys is None:
            print(f"[WARN] No peptide atoms found in {fdir} (chain {args.peptide_chain}). Skipping.")
            continue

        # Align arrays to same ordering as template; skip models that can't match
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
        top = pick_top_cluster(comps, poses)
        centroid_idx = pick_centroid(top, pair)
        if centroid_idx < 0: centroid_idx = top[0]

        # If centroid has no metrics, try best-scored WITH metrics (but clustering did use ALL models)
        def has_metrics(i: int) -> bool:
            return poses[i].score_I_sc is not None or poses[i].score_total is not None or poses[i].score_reweighted is not None
        if not has_metrics(centroid_idx):
            candidates = [i for i in top if has_metrics(i)]
            if candidates:
                def keyfun(i):
                    rw = poses[i].score_reweighted
                    ts = poses[i].score_total
                    return (rw if rw is not None else float('inf'),
                            ts if ts is not None else float('inf'))
                centroid_idx = min(candidates, key=keyfun)

        centroid = poses[centroid_idx]

        # Energies for representative
        iface_reu = centroid.score_I_sc
        total_reu = centroid.score_total
        if args.interfaceanalyzer and iface_reu is None:
            iface_reu, hb_rosetta = run_interface_analyzer(args.interfaceanalyzer, centroid.path, args.receptor_chain, args.peptide_chain)
        else:
            hb_rosetta = None
        if args.score_jd2 and total_reu is None:
            total_reu = run_score_jd2(args.score_jd2, centroid.path)

        # H-bonds (prefer IA value if było; w przeciwnym razie geometria)
        hbonds = hb_rosetta if hb_rosetta is not None else count_hbonds_simple(centroid.path, args.receptor_chain, args.peptide_chain)

        key_contacts = format_key_contacts(centroid.path, args.receptor_chain, args.peptide_chain, args.contact_cutoff)
        potential = decide_potential(iface_reu, hbonds, key_contacts)

        # RMSD vs reference (if provided)
        rmsd_ref = None
        if args.ref_pdb:
            rmsd_ref = rmsd_between_chains(centroid.path, args.peptide_chain, args.ref_pdb, ref_pep_chain, mode=('heavy' if args.rmsd_mode=='heavy' else 'backbone'))

        # Reference energies (optional)
        i_sc_ref = ref_i_sc
        total_ref = ref_total
        # If not precomputed, compute per-ligand on demand
        if args.ref_pdb and args.score_ref and (i_sc_ref is None or total_ref is None):
            if args.interfaceanalyzer and i_sc_ref is None:
                i_sc_ref, _ = run_interface_analyzer(args.interfaceanalyzer, args.ref_pdb, ref_rec_chain, ref_pep_chain)
            if args.score_jd2 and total_ref is None:
                total_ref = run_score_jd2(args.score_jd2, args.ref_pdb)

        di_sc = (iface_reu - i_sc_ref) if (iface_reu is not None and i_sc_ref is not None) else None
        dtotal = (total_reu - total_ref) if (total_reu is not None and total_ref is not None) else None

        rows.append(RowOut(
            ligand=ligand_label, species=species,
            iface_energy_reu=iface_reu, total_energy_reu=total_reu,
            hbonds=hbonds, key_contacts=key_contacts, inhibitor_potential=potential,
            rmsd_vs_ref=rmsd_ref, i_sc_ref=i_sc_ref, di_sc=di_sc,
            total_ref=total_ref, dtotal=dtotal
        ))

        # Per-folder CSV
        single_out = os.path.join(args.outdir, f"{ligand_label}.csv")
        with open(single_out, 'w', newline='') as sfh:
            w1 = csv.writer(sfh)
            w1.writerow(["Ligand","Gatunek","Energia interfejsu (REU)","Energia całkowita (REU)","Wiązania wodorowe","Reszty w kontakcie","Potencjał inhibitorowy","RMSD vs ref (Å)","I_sc ref (REU)","ΔI_sc (REU)","Total ref (REU)","ΔTotal (REU)"])
            w1.writerow([ligand_label, species, fmt(iface_reu), fmt(total_reu), fmt_int(hbonds), key_contacts, potential, fmt(rmsd_ref), fmt(i_sc_ref), fmt(di_sc), fmt(total_ref), fmt(dtotal)])

        if args.debug:
            print(f"[DBG] {ligand_label}: models={len(model_paths)}, clusters={len(comps)}, top_size={len(top)}, centroid={poses[centroid_idx].name}, I_sc={iface_reu}, total={total_reu}, RMSD_ref={rmsd_ref}")

    # Summary CSV
    with open(args.out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(["Ligand","Gatunek","Energia interfejsu (REU)","Energia całkowita (REU)","Wiązania wodorowe","Reszty w kontakcie","Potencjał inhibitorowy","RMSD vs ref (Å)","I_sc ref (REU)","ΔI_sc (REU)","Total ref (REU)","ΔTotal (REU)"])
        for r in rows:
            w.writerow([r.ligand, r.species, fmt(r.iface_energy_reu), fmt(r.total_energy_reu), fmt_int(r.hbonds), r.key_contacts, r.inhibitor_potential, fmt(r.rmsd_vs_ref), fmt(r.i_sc_ref), fmt(r.di_sc), fmt(r.total_ref), fmt(r.dtotal)])
    print(f"Wrote {args.out} with {len(rows)} rows.")

def fmt(x: Optional[float]) -> str:
    if x is None or (isinstance(x,float) and (math.isnan(x) or math.isinf(x))):
        return ""
    return f"{x:.2f}"

def fmt_int(x: Optional[int]) -> str:
    if x is None:
        return ""
    return str(int(x))

if __name__ == '__main__':
    main()


"""
uruchamianie z ref:


python /home/marta/Pulpit/biobloom_scripts/molecular_docking/flexpepdock_analysis.py \
  --root /home/marta/Pulpit/docking_output \
  --glob "rosie-*/output" \
  --receptor-chain A \
  --peptide-chain P \
  --rmsd-cutoff 2.0 \
  --contact-cutoff 4.0 \
  --rmsd-mode backbone \
  --ref-pdb /home/marta/Pulpit/docking_input/ACE_structures/1O86.pdb \
  --ref-receptor-chain A \
  --ref-peptide-chain A \
  --score-ref \
  --interfaceanalyzer "$(which InterfaceAnalyzer.linuxgccrelease)" \
  --score_jd2 "$(which score_jd2.linuxgccrelease)" \
  --outdir /home/marta/Pulpit/docking_results_analysis \
  --out /home/marta/Pulpit/docking_results_analysis/summary.csv \
  --debug

"""