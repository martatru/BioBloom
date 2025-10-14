#!/usr/bin/env python3
"""
FlexPepDock post-processing: cluster poses -> pick representative -> compute/interface scores -> export poster-ready table.

- Klasteryzacja po RMSD (backbone peptydu; domyślnie 2.0 Å)
- Wybór top klastra (najniższa mediana reweighted_sc; tie: większy klaster)
- Reprezentant = centroid (min. średni RMSD). Jeśli nie ma metryk, zamiana na najlepiej punktowany model z metrykami.
- Energia interfejsu: I_sc (albo dG_separated z InterfaceAnalyzer, jeśli podasz binarkę)
- Energia całkowita: total_score (albo score_jd2, jeśli podasz binarkę)
- Wiązania wodorowe: prosta metryka geometryczna N/O…N/O ≤ 3.5 Å (gdy brak IA)
- Kontakty: reszty ACE w odległości ≤ 4.0 Å od atomów ciężkich peptydu
- Zapis: CSV zbiorczy + osobne CSV na folder (nazwane jak folder nadrzędny)

Przykład uruchomienia:
python flexpepdock_analysis.py \
  --root /home/marta/Pulpit/docking_output \
  --glob "rosie-*/output" \
  --receptor-chain A --peptide-chain P \
  --rmsd-cutoff 2.0 --contact-cutoff 4.0 \
  --outdir /home/marta/Pulpit/docking_results_analysis \
  --out /home/marta/Pulpit/docking_results_analysis/summary.csv \
  --debug
"""
from __future__ import annotations
import argparse, csv, glob, math, os, re, subprocess, sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

try:
    import numpy as np
except ImportError:
    print("This script requires numpy. Please install it (e.g., pip install numpy)", file=sys.stderr)
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

def parse_heavy_atoms(pdb_path: str, chain_id: str) -> List[np.ndarray]:
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
        cmd = [exe, '-s', pdb, '-interface', f'{receptor_chain}_{peptide_chain}', '-pack_input', 'false', '-pack_separated', 'false', '-out:file:scorefile', 'ia_tmp.sc']
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
    pep_heavy = parse_heavy_atoms(pdb_path, peptide_chain)
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
    # fallback: 3 nearest residues
    dist_map = []
    for resi, atoms in rec_atoms.items():
        A = np.vstack(atoms)
        dmin = np.linalg.norm(pep[:,None,:] - A[None,:,:], axis=2).min()
        dist_map.append((resi, dmin))
    dist_map.sort(key=lambda x: x[1])
    return ", ".join(ACE_KEY_NAMES.get(r, f"Res{r}") for r,_ in dist_map[:3]) if dist_map else "—"

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

# -------------------- Main --------------------
def main():
    ap = argparse.ArgumentParser(description="Cluster FlexPepDock outputs and export poster-ready table.")
    ap.add_argument('--root', required=True, help='Parent directory that contains per-peptide folders')
    ap.add_argument('--glob', default='*/', help='Subfolder glob under root to scan (default: */)')
    ap.add_argument('--receptor-chain', default='A')
    ap.add_argument('--peptide-chain', default='P')
    ap.add_argument('--rmsd-cutoff', type=float, default=2.0)
    ap.add_argument('--contact-cutoff', type=float, default=4.0)
    ap.add_argument('--interfaceanalyzer', default=None, help='Path to InterfaceAnalyzer.linuxgccrelease')
    ap.add_argument('--score_jd2', default=None, help='Path to score_jd2.linuxgccrelease')
    ap.add_argument('--species-map', default=None, help='CSV two columns: ligand_id,species list')
    ap.add_argument('--outdir', default='/home/marta/Pulpit/docking_results_analysis')
    ap.add_argument('--out', default='results.csv')
    ap.add_argument('--debug', action='store_true')
    args = ap.parse_args()

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

    for fdir in folders:
        # Ligand = parent folder (one level above 'output')
        parent = os.path.basename(os.path.dirname(os.path.normpath(fdir)))
        ligand_label = parent
        species = ligand_species.get(ligand_label, "")

        # Wczytaj tabele metryk
        metrics = {}
        for cand in ('top10_scores.tsv','filtered_scores.tsv','metrics.csv'):
            p = os.path.join(fdir, cand)
            if os.path.isfile(p):
                metrics.update(read_metrics_table(p))

        # Wszystkie modele
        all_models = sorted(glob.glob(os.path.join(fdir, 'model_*.pdb')))
        if not all_models:
            all_models = sorted(glob.glob(os.path.join(fdir, 'complex.ppk_*.pdb')))
        if not all_models:
            print(f"[WARN] No models in {fdir}")
            continue

        # --- wybór modeli z metrykami ---
        metric_pdb_names = {k for k in metrics.keys() if isinstance(k, str) and k.endswith('.pdb')}

        max_idx_seen = -1
        for v in metrics.values():
            if isinstance(v, dict) and 'index' in v:
                try:
                    max_idx_seen = max(max_idx_seen, int(float(v['index'])))
                except Exception:
                    pass
        n_by_index = (max_idx_seen + 1) if max_idx_seen >= 0 else 0

        models_with_metrics = [p for p in all_models if os.path.basename(p) in metric_pdb_names]
        if (len(models_with_metrics) < 2) and (n_by_index >= 2):
            models_with_metrics = all_models[:n_by_index]

        if len(models_with_metrics) >= 2:
            model_paths = models_with_metrics
        else:
            model_paths = all_models
            if len(models_with_metrics) == 0:
                print(f"[WARN] No per-model metrics matched in {fdir}; energies may be empty.")

        if args.debug:
            print(f"[DBG] {parent}: all={len(all_models)}, with_metrics_by_name={len([p for p in all_models if os.path.basename(p) in metric_pdb_names])}, n_by_index={n_by_index}, clustered={len(model_paths)}")

        # Pose objects
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
            bb, idx = parse_pdb_backbone_chain(mp, args.peptide_chain)
            m.bb_coords = bb; m.bb_index = idx
            poses.append(m)

        # comparable backbone set
        template_idx = None
        for p in poses:
            if p.bb_coords is not None and p.bb_coords.size>0:
                template_idx = p.bb_index; break
        if template_idx is None:
            print(f"[WARN] No peptide backbone atoms found in {fdir} (chain {args.peptide_chain}). Skipping.")
            continue

        def extract_matching(p: PoseInfo) -> Optional[np.ndarray]:
            if not p.bb_index:
                return None
            idx_map = {ra:i for i,ra in enumerate(p.bb_index)}
            rowsB = []
            for ra in template_idx:
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

        # Jeśli centroid bez metryk -> wybierz najlepiej punktowany z metrykami
        def has_metrics(i: int) -> bool:
            return poses[i].score_I_sc is not None or poses[i].score_total is not None or poses[i].score_reweighted is not None

        if not has_metrics(centroid_idx):
            candidates = [i for i in top if has_metrics(i)]
            if not candidates:
                candidates = [i for i in range(n) if has_metrics(i)]
            if candidates:
                def keyfun(i):
                    rw = poses[i].score_reweighted
                    ts = poses[i].score_total
                    return (rw if rw is not None else float('inf'),
                            ts if ts is not None else float('inf'))
                centroid_idx = min(candidates, key=keyfun)

        centroid = poses[centroid_idx]

        iface_reu = None; hbonds = None; total_reu = None
        if args.interfaceanalyzer:
            iface_reu, hbonds = run_interface_analyzer(args.interfaceanalyzer, centroid.path, args.receptor_chain, args.peptide_chain)
        if args.score_jd2:
            total_reu = run_score_jd2(args.score_jd2, centroid.path)
        if iface_reu is None and centroid.score_I_sc is not None:
            iface_reu = centroid.score_I_sc
        if total_reu is None and centroid.score_total is not None:
            total_reu = centroid.score_total
        if hbonds is None:
            hbonds = count_hbonds_simple(centroid.path, args.receptor_chain, args.peptide_chain)

        key_contacts = format_key_contacts(centroid.path, args.receptor_chain, args.peptide_chain, args.contact_cutoff)
        potential = decide_potential(iface_reu, hbonds, key_contacts)

        rows.append(RowOut(ligand=ligand_label, species=species,
                           iface_energy_reu=iface_reu, total_energy_reu=total_reu,
                           hbonds=hbonds, key_contacts=key_contacts, inhibitor_potential=potential))

        # per-folder CSV (nazwane jak folder nadrzędny)
        try:
            if args.outdir:
                single_out = os.path.join(args.outdir, f"{ligand_label}.csv")
                with open(single_out, 'w', newline='') as sfh:
                    w1 = csv.writer(sfh)
                    w1.writerow(["Ligand", "Gatunek", "Energia interfejsu (REU)", "Energia całkowita (REU)", "Wiązania wodorowe", "Reszty w kontakcie", "Potencjał inhibitorowy"])
                    w1.writerow([ligand_label, species, fmt(iface_reu), fmt(total_reu), fmt_int(hbonds), key_contacts, potential])
        except Exception as e:
            print(f"[WARN] Could not write per-folder CSV for {ligand_label}: {e}")

    # summary CSV
    with open(args.out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(["Ligand", "Gatunek", "Energia interfejsu (REU)", "Energia całkowita (REU)", "Wiązania wodorowe", "Reszty w kontakcie", "Potencjał inhibitorowy"])
        for r in rows:
            w.writerow([r.ligand, r.species, fmt(r.iface_energy_reu), fmt(r.total_energy_reu), fmt_int(r.hbonds), r.key_contacts, r.inhibitor_potential])
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
