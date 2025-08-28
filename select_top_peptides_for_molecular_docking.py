"""
how to run:

python select_top_peptides_for_molecular_docking.py --input (path)

paths:

/home/marta/Desktop/admet_output/t_suecica/t_suecica.csv
/home/marta/Desktop/admet_output/m_salina/m_salina_merged.csv
/home/marta/Desktop/admet_output/a_platensis/a_platensis_merged.csv

"""
#!/usr/bin/env python3
import argparse, math, os, sys, csv
import pandas as pd
import numpy as np

OUTDIR  = "/home/marta/Desktop/peptides_chosen_for_docking/smiles_lists"
ALL_OUT = "peptides_ranked_all.csv"
PASS_OUT= "peptides_ranked_pass.csv"

def detect_sep(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        header = f.readline()
    for sep in (",",";","\t","|"):
        if header.count(sep) >= 2:
            return sep
    return ","  # sensible default for ADMETlab

def read_csv_robust(path):
    sep = detect_sep(path)
    # 1) try C engine (fast, no small-field cap)
    try:
        return pd.read_csv(path, sep=sep, engine="c", low_memory=False)
    except Exception:
        # 2) raise Python csv field-size cap and retry with Python engine
        try:
            csv.field_size_limit(sys.maxsize)
        except OverflowError:
            csv.field_size_limit(2**31 - 1)
        return pd.read_csv(path, sep=sep, engine="python", low_memory=False)

def to_float(x):
    if pd.isna(x): return np.nan
    if isinstance(x, str):
        s = x.strip().replace('%','').replace(',','.')
        v = s.lower()
        if v in ('high','h','yes','y','true'):  return 80.0
        if v in ('low','l','no','n','false'):   return 20.0
        try: return float(s)
        except: return np.nan
    return float(x)

def to_pct(x):
    v = to_float(x)
    if pd.isna(v): return np.nan
    return v*100.0 if 0.0 <= v <= 1.0 else v

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
    for c in ("f30","f50","f20"):
        if c in row.index:
            v = to_pct(row[c])
            if not pd.isna(v):
                return v
    return np.nan

def build_with_derivatives(df):
    out = df.copy()

    out["MW_"]    = df["MW"].map(to_float)
    out["nHD_"]   = df["nHD"].map(to_float)
    out["nHA_"]   = df["nHA"].map(to_float)
    out["logP_"]  = df["logP"].map(to_float)

    conds = [
        out["MW_"]    < 500.0,
        out["nHD_"]   < 5.0,
        out["nHA_"]   < 10.0,
        out["logP_"]  < 5.0
    ]
    out["Lipinski_passes_"] = sum([c.fillna(False) for c in conds])

    out["HIA_pct_"] = df["hia"].map(to_pct)
    out["Fxx_pct_"] = df.apply(pick_ba_prob, axis=1)

    out["logVDss_"] = df["logVDss"].map(to_float)
    out["Vd_Lkg_"]  = np.power(10.0, out["logVDss_"])

    out["t_half_h_"]  = df["t0.5"].map(to_float)
    out["LD50_mgkg_"] = df["LD50_oral"].map(to_float)

    out["C1_Lipinski_3of4"] = out["Lipinski_passes_"] >= 3
    out["C2_HIA_gt30"]      = out["HIA_pct_"] > 30.0
    out["C3_BA_gt30_prob"]  = out["Fxx_pct_"] >= 50.0
    out["C4_Vd_range"]      = (out["Vd_Lkg_"] >= 0.04) & (out["Vd_Lkg_"] <= 20.0)
    out["C5_t12_ge_0p5h"]   = out["t_half_h_"] >= 0.5
    out["C6_LD50_gt500"]    = out["LD50_mgkg_"] > 500.0

    out["Pass_WZ"] = (
        out["C1_Lipinski_3of4"] &
        out["C2_HIA_gt30"] &
        out["C3_BA_gt30_prob"] &
        out["C4_Vd_range"] &
        out["C5_t12_ge_0p5h"] &
        out["C6_LD50_gt500"]
    )

    def score_row(r):
        lip = (r["Lipinski_passes_"] / 4.0)
        hia = normalize_01(r["HIA_pct_"], 0.0, 100.0)
        ba  = normalize_01(r["Fxx_pct_"], 50.0, 100.0)
        vd  = normalize_logrange(r["Vd_Lkg_"], 0.04, 20.0)
        t12 = normalize_01(r["t_half_h_"], 0.5, 24.0)
        tox = normalize_01(r["LD50_mgkg_"], 500.0, 5000.0)
        return (1.0*lip + 1.0*hia + 1.0*ba + 0.5*vd + 0.5*t12 + 1.0*tox)

    out["ADMET_score_"] = out.apply(score_row, axis=1)
    out["Hsum_"] = (out["nHD_"].fillna(99) + out["nHA_"].fillna(99))

    def failure_str(r):
        fails = []
        if not r["C1_Lipinski_3of4"]: fails.append("Lipinski<3/4")
        if not r["C2_HIA_gt30"]:      fails.append("HIA<=30%")
        if not r["C3_BA_gt30_prob"]:  fails.append("BA>30% prob<50%")
        if not r["C4_Vd_range"]:      fails.append("Vd out of 0.04–20")
        if not r["C5_t12_ge_0p5h"]:   fails.append("t1/2<0.5h")
        if not r["C6_LD50_gt500"]:    fails.append("LD50<=500")
        return "; ".join(fails)
    out["Fail_Reasons"] = out.apply(failure_str, axis=1)

    return out

def main():
    ap = argparse.ArgumentParser(description="Rank ALL peptides (passers first) from an ADMETlab 3.0 CSV.")
    ap.add_argument("--input", required=True, help="Path to species CSV (single header line).")
    args = ap.parse_args()

    df = read_csv_robust(args.input)
    ranked = build_with_derivatives(df)

    ranked = ranked.sort_values(
        by=["Pass_WZ","ADMET_score_","Lipinski_passes_","MW_","Hsum_"],
        ascending=[False, False, True, True, True]
    )

    os.makedirs(OUTDIR, exist_ok=True)
    all_path  = os.path.join(OUTDIR, ALL_OUT)
    pass_path = os.path.join(OUTDIR, PASS_OUT)

    ranked.to_csv(all_path, index=False)
    ranked[ranked["Pass_WZ"]].to_csv(pass_path, index=False)

    print(f"[OK] Ranked ALL peptides → {all_path}")
    print(f"[OK] Ranked PASS-only   → {pass_path}")
    print("Sorting: Pass_WZ (True first) → ADMET_score_ (desc) → Lipinski_passes_ → MW_ → (nHD_+nHA_)")

if __name__ == "__main__":
    main()
