"""
how to run:

python select_top_peptides_for_molecular_docking.py --input /home/marta/Desktop/admet_output/t_suecica/t_suecica.csv
python select_top_peptides_for_molecular_docking.py --input /home/marta/Desktop/admet_output/m_salina/m_salina_merged.csv
python select_top_peptides_for_molecular_docking.py --input /home/marta/Desktop/admet_output/a_platensis/a_platensis_merged.csv

"""

import argparse, math, os, sys, csv, re
import pandas as pd
import numpy as np

DEFAULT_OUTROOT = "/home/marta/Desktop/peptides_chosen_for_docking/smiles_lists"

ALL_OUT   = "peptides_ranked_all.csv"
PASS_OUT  = "peptides_ranked_pass.csv"
TOP50_OUT = "top50_peptides.csv"
STATS_OUT = "criterion_pass_rates.txt"

# ---------- IO helpers ----------
def detect_sep(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        sample = f.read(1024*64)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t|")
        return dialect.delimiter
    except Exception:
        first = sample.splitlines()[0] if sample else ""
        for sep in (",",";","\t","|"):
            if first.count(sep) >= 2:
                return sep
    return ","

def read_csv_robust(path):
    sep = detect_sep(path)
    try:
        return pd.read_csv(path, sep=sep, engine="c", low_memory=False)
    except Exception:
        try:
            csv.field_size_limit(sys.maxsize)
        except OverflowError:
            csv.field_size_limit(2**31 - 1)
        return pd.read_csv(path, sep=sep, engine="python", low_memory=False)

# ---------- parsing / normalization ----------
def smart_float(x):
    if pd.isna(x): return np.nan
    if isinstance(x, (int,float,np.integer,np.floating)): return float(x)
    s = str(x).strip()
    if not s or s == "-": return np.nan
    s = s.replace(">","").replace("<","")
    if "," in s and "." in s:
        s = s.replace(",","")     # 1,234.56 -> 1234.56
    elif "," in s and "." not in s:
        s = s.replace(",",".")    # 3,5 -> 3.5
    m = re.search(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', s)
    return float(m.group(0)) if m else np.nan

def to_pct(x):
    v = smart_float(x)
    if pd.isna(v): return np.nan
    return v*100.0 if 0.0 <= v <= 1.0 else v

def is_negative_flag(x, thr=0.5):
    """True = 'bezpieczny/negatywny'. Obsługuje napisy i wartości liczbowe (p<0.5)."""
    if pd.isna(x): return False
    if isinstance(x, (int,float,np.integer,np.floating)):
        return float(x) < thr
    s = str(x).strip().lower()
    if s in ("neg","negative","no","false","0","low","non-blocker","nonblocker"): return True
    if s in ("pos","positive","yes","true","1","blocker","high"): return False
    # numeric in string
    try:
        v = float(s.replace(",","."))
        return v < thr
    except:  # unknown token
        return False

def is_positive_flag(x, thr=0.5):
    """True = 'problem/pozytywny' (np. PAINS = TAK)."""
    if pd.isna(x): return False
    if isinstance(x, (int,float,np.integer,np.floating)):
        return float(x) >= thr
    s = str(x).strip().lower()
    if s in ("pos","positive","yes","true","1","blocker","high"): return True
    if s in ("neg","negative","no","false","0","low","non-blocker","nonblocker"): return False
    try:
        v = float(s.replace(",","."))
        return v >= thr
    except:
        return False

def normalize_01(x, lo, hi):
    if pd.isna(x): return 0.0
    x = min(max(x, lo), hi)
    return 0.0 if hi==lo else (x - lo) / (hi - lo)

def normalize_logrange(x, lo, hi):
    if pd.isna(x) or x<=0 or lo<=0 or hi<=0: return 0.0
    lx, llo, lhi = math.log10(x), math.log10(lo), math.log10(hi)
    lx = min(max(lx, llo), lhi)
    return (lx - llo) / (lhi - llo)

def pick_ba_prob(row):
    # prefer F30, potem F50, F20 – wszystko w % (0–100)
    for c in ("f30","f50","f20"):
        if c in row.index:
            v = to_pct(row[c])
            if not pd.isna(v):
                return v
    return np.nan

# ---------- core build ----------
def build_with_derivatives(df, ba_prob_cut=50.0):
    out = df.copy()

    # fizykochemia
    out["MW_"]   = df["MW"].map(smart_float)
    out["nHD_"]  = df["nHD"].map(smart_float)
    out["nHA_"]  = df["nHA"].map(smart_float)
    out["logP_"] = df["logP"].map(smart_float)

    # Lipiński 3/4
    conds = [
        out["MW_"]   < 500.0,
        out["nHD_"]  < 5.0,
        out["nHA_"]  < 10.0,
        out["logP_"] < 5.0
    ]
    out["Lipinski_passes_"] = sum([c.fillna(False) for c in conds])

    # Absorpcja (HIA, BA)
    out["HIA_pct_"] = df["hia"].map(to_pct)
    out["Fxx_pct_"] = df.apply(pick_ba_prob, axis=1)

    # Dystrybucja (Vd z logVDss)
    out["logVDss_"] = df["logVDss"].map(smart_float)
    out["Vd_Lkg_"]  = np.power(10.0, out["logVDss_"])

    # Eliminacja
    out["t_half_h_"]  = df["t0.5"].map(smart_float)

    # LD50 (jeśli dostępne)
    out["LD50_mgkg_"] = df["LD50_oral"].map(smart_float)
    out["LD50_available_"] = out["LD50_mgkg_"].notna()

    # Bezpieczeństwo (panel WoE)
    out["Ames_neg_"] = df["Ames"].apply(is_negative_flag) if "Ames" in df.columns else False
    out["DILI_neg_"] = df["DILI"].apply(is_negative_flag) if "DILI" in df.columns else False
    # prefer hERG-10um jeśli jest; inaczej hERG
    if "hERG-10um" in df.columns:
        out["hERG_neg_"] = df["hERG-10um"].apply(is_negative_flag)
    elif "hERG" in df.columns:
        out["hERG_neg_"] = df["hERG"].apply(is_negative_flag)
    else:
        out["hERG_neg_"] = False

    # Interferencje (wszystkie muszą być negatywne)
    inter_cols = [c for c in ("PAINS","Reactive","Aggregators","Promiscuous") if c in df.columns]
    if inter_cols:
        inter_any = np.zeros(len(df), dtype=bool)
        for c in inter_cols:
            inter_any |= df[c].apply(is_positive_flag).values
        out["Interference_any_"] = inter_any
    else:
        out["Interference_any_"] = False

    # --------- Kryteria metodologii ---------
    out["C1_Lipinski_3of4"] = out["Lipinski_passes_"] >= 3
    out["C2_HIA_gt30"]      = out["HIA_pct_"] > 30.0
    out["C3_BA_gt30_prob"]  = out["Fxx_pct_"] >= ba_prob_cut
    out["C4_Vd_range"]      = (out["Vd_Lkg_"] >= 0.04) & (out["Vd_Lkg_"] <= 20.0)
    out["C5_t12_ge_0p5h"]   = out["t_half_h_"] >= 0.5

    # (6) Bezpieczeństwo: Ames−, DILI−, hERG− oraz brak interferencji
    out["C6_safety_panel"]  = out["Ames_neg_"] & out["DILI_neg_"] & out["hERG_neg_"] & (~out["Interference_any_"])

    # LD50 logika WoE:
    # Jeśli LD50 jest dostępne → wymagamy LD50>500 ORAZ C6_safety_panel
    # Jeśli LD50 brak           → wymagamy C6_safety_panel (bez warunku LD50)
    out["C6_final"] = np.where(out["LD50_available_"],
                               (out["LD50_mgkg_"] > 500.0) & out["C6_safety_panel"],
                               out["C6_safety_panel"])

    crit = ["C1_Lipinski_3of4","C2_HIA_gt30","C3_BA_gt30_prob","C4_Vd_range","C5_t12_ge_0p5h","C6_final"]
    out["Pass_Count"] = sum([out[c].fillna(False) for c in crit])
    out["Pass_WZ"]    = out["Pass_Count"] == 6

    # --------- Ranking (kompozyt) ---------
    # Bez LD50 w składowej (stosujemy panel bezpieczeństwa)
    safety_avg = (out["Ames_neg_"].astype(float) + out["DILI_neg_"].astype(float) + out["hERG_neg_"].astype(float)) / 3.0
    interf_pen = out["Interference_any_"].astype(float)  # 1 jeśli problem

    def score_row(r):
        lip = (r["Lipinski_passes_"] / 4.0)
        hia = normalize_01(r["HIA_pct_"], 0.0, 100.0)
        ba  = normalize_01(r["Fxx_pct_"], 50.0, 100.0)
        vd  = normalize_logrange(r["Vd_Lkg_"], 0.04, 20.0)
        t12 = normalize_01(r["t_half_h_"], 0.5, 24.0)
        saf = r["safety_avg_"]
        pen = r["interf_pen_"]
        return (1.0*lip + 1.0*hia + 1.0*ba + 0.5*vd + 0.5*t12 + 1.0*saf - 0.5*pen)

    out["safety_avg_"] = safety_avg
    out["interf_pen_"] = interf_pen
    out["ADMET_score_"] = out.apply(score_row, axis=1)

    # Tie helper
    out["Hsum_"] = (out["nHD_"].fillna(99) + out["nHA_"].fillna(99))
    return out

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Rank ALL peptides per new methodology (ADMETlab 3.0).")
    ap.add_argument("--input", required=True, help="Path to species CSV (single header line).")
    ap.add_argument("--ba_prob_cut", type=float, default=50.0,
                    help="Probability cutoff (%) that oral BA > 30% (default 50).")
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Root folder where a subfolder named after the input file will be created.")
    args = ap.parse_args()

    df = read_csv_robust(args.input)
    ranked = build_with_derivatives(df, ba_prob_cut=args.ba_prob_cut)

    # Sort: passers first, then by score, then tie-breakers
    ranked = ranked.sort_values(
        by=["Pass_WZ","ADMET_score_","Lipinski_passes_","MW_","Hsum_"],
        ascending=[False, False, True, True, True]
    )

    base = os.path.splitext(os.path.basename(args.input))[0]
    outdir = os.path.join(args.outroot, base)
    os.makedirs(outdir, exist_ok=True)

    all_path  = os.path.join(outdir, ALL_OUT)
    pass_path = os.path.join(outdir, PASS_OUT)
    top_path  = os.path.join(outdir, TOP50_OUT)
    stats_path= os.path.join(outdir, STATS_OUT)

    ranked.to_csv(all_path, index=False)
    strict_pass = ranked[ranked["Pass_WZ"]]
    strict_pass.to_csv(pass_path, index=False)

    # Always produce a Top-50:
    if len(strict_pass) > 0:
        top = strict_pass.head(min(50, len(strict_pass))); source = "strict (6/6)"
    else:
        relaxed = ranked[ranked["Pass_Count"] >= 5]
        if len(relaxed) > 0:
            top = relaxed.head(min(50, len(relaxed))); source = "relaxed (>=5/6)"
        else:
            top = ranked.head(min(50, len(ranked)));   source = "overall"
    top.to_csv(top_path, index=False)

    # Stats & audit
    n = len(ranked)
    with open(stats_path, "w") as f:
        f.write(f"Total rows: {n}\n")
        f.write(f"Strict pass (all 6): {int(ranked['Pass_WZ'].sum())}\n")
        f.write(f"C1 Lipinski ≥3/4: {int((ranked['Lipinski_passes_'] >= 3).sum())}\n")
        f.write(f"C2 HIA>30%: {int((ranked['HIA_pct_'] > 30).sum())}\n")
        f.write(f"C3 BA>30% prob ≥{args.ba_prob_cut}%: {int((ranked['Fxx_pct_'] >= args.ba_prob_cut).sum())}\n")
        f.write(f"C4 Vd in [0.04,20] L/kg: {int(((ranked['Vd_Lkg_'] >= 0.04) & (ranked['Vd_Lkg_'] <= 20)).sum())}\n")
        f.write(f"C5 t1/2 ≥0.5 h: {int((ranked['t_half_h_'] >= 0.5).sum())}\n")
        f.write(f"C6 safety panel (Ames-, DILI-, hERG- & no PAINS/Reactive/Aggregators/Promiscuous): {int(ranked['C6_safety_panel'].sum())}\n")
        if "LD50_mgkg_" in ranked.columns:
            ld = ranked["LD50_mgkg_"].dropna()
            f.write(f"LD50 available: {len(ld)} / {n}\n")
            f.write(f"LD50>500 among available: {int((ld > 500).sum())}\n")
        f.write("\nInterference columns present: ")
        f.write(", ".join([c for c in ("PAINS","Reactive","Aggregators","Promiscuous") if c in df.columns]) or "none")
        f.write("\n")

    print(f"[OK] Out folder: {outdir}")
    print(f"[OK] Ranked ALL → {all_path}")
    print(f"[OK] Strict PASS → {pass_path}  (rows: {len(strict_pass)})")
    print(f"[OK] Top-50 ({source}) → {top_path}")
    print(f"[OK] Stats → {stats_path}")

if __name__ == "__main__":
    main()
