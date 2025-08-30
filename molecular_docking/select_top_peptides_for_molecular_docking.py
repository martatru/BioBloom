#!/usr/bin/env python3
"""
Select peptides with top ADMET properties from AdmetLab 3.0 batch CSV
and save ranked/filtered outputs.

How to run:
    python select_top_peptides_for_molecular_docking.py \
        --input /path/to/species_admetlab.csv \
        --ba_prob_cut 50 \
        --outroot /path/to/output_root
"""

from __future__ import annotations

import argparse, math, os, sys, csv, re
from pathlib import Path
import pandas as pd
import numpy as np

# --- Defaults & output names --------------------------------------------------

DEFAULT_OUTROOT = "/home/marta/Desktop/peptides_chosen_for_docking/smiles_lists"

ALL_OUT        = "peptides_ranked_all.csv"
PASS_OUT       = "peptides_ranked_pass.csv"
PASS_PLUS5_OUT = "peptides_pass_plus_5of6_ranked.csv"   # main deliverable
STATS_OUT      = "criterion_pass_rates.txt"

# --- IO helpers ---------------------------------------------------------------

def detect_sep(path: Path) -> str:
    sample = path.read_text(encoding="utf-8", errors="ignore")[: 1024 * 64]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t|")
        return dialect.delimiter
    except Exception:
        first = sample.splitlines()[0] if sample else ""
        for sep in (",", ";", "\t", "|"):
            if first.count(sep) >= 2:
                return sep
    return ","

def read_csv_robust(path: Path) -> pd.DataFrame:
    sep = detect_sep(path)
    try:
        return pd.read_csv(path, sep=sep, engine="c", low_memory=False)
    except Exception:
        try:
            csv.field_size_limit(sys.maxsize)
        except OverflowError:
            csv.field_size_limit(2**31 - 1)
        return pd.read_csv(path, sep=sep, engine="python", low_memory=False)

# --- Parsing / normalization --------------------------------------------------

def smart_float(x):
    if pd.isna(x): return np.nan
    if isinstance(x, (int, float, np.integer, np.floating)): return float(x)
    s = str(x).strip()
    if not s or s == "-": return np.nan
    s = s.replace(">", "").replace("<", "")
    if "," in s and "." in s:
        s = s.replace(",", "")       # 1,234.56 -> 1234.56
    elif "," in s and "." not in s:
        s = s.replace(",", ".")      # 3,5 -> 3.5
    m = re.search(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', s)
    return float(m.group(0)) if m else np.nan

def to_pct(x):
    v = smart_float(x)
    if pd.isna(v): return np.nan
    return v * 100.0 if 0.0 <= v <= 1.0 else v

def is_negative_flag(x, thr=0.5):
    """True = safe/negative (e.g. Ames-, hERG-). Supports strings or numeric p-values (< thr)."""
    if pd.isna(x): return False
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x) < thr
    s = str(x).strip().lower()
    if s in ("neg", "negative", "no", "false", "0", "low", "non-blocker", "nonblocker"): return True
    if s in ("pos", "positive", "yes", "true", "1", "blocker", "high"): return False
    try:
        v = float(s.replace(",", "."))
        return v < thr
    except Exception:
        return False

def is_positive_flag(x, thr=0.5):
    """True = problematic/positive (e.g., PAINS = yes)."""
    if pd.isna(x): return False
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x) >= thr
    s = str(x).strip().lower()
    if s in ("pos", "positive", "yes", "true", "1", "blocker", "high"): return True
    if s in ("neg", "negative", "no", "false", "0", "low", "non-blocker", "nonblocker"): return False
    try:
        v = float(s.replace(",", "."))
        return v >= thr
    except Exception:
        return False

def normalize_01(x, lo, hi):
    if pd.isna(x): return 0.0
    x = min(max(x, lo), hi)
    return 0.0 if hi == lo else (x - lo) / (hi - lo)

def normalize_logrange(x, lo, hi):
    if pd.isna(x) or x <= 0 or lo <= 0 or hi <= 0: return 0.0
    lx, llo, lhi = math.log10(x), math.log10(lo), math.log10(hi)
    lx = min(max(lx, llo), lhi)
    return (lx - llo) / (lhi - llo)

def pick_ba_prob(row):
    # prefer F30, then F50, F20 — all as percentages (0–100)
    for c in ("f30", "f50", "f20"):
        if c in row.index:
            v = to_pct(row[c])
            if not pd.isna(v):
                return v
    return np.nan

# --- Core build ---------------------------------------------------------------

def build_with_derivatives(df: pd.DataFrame, ba_prob_cut: float = 50.0) -> pd.DataFrame:
    out = df.copy()

    # Physicochemical
    out["MW_"]   = df["MW"].map(smart_float)   if "MW"   in df.columns else np.nan
    out["nHD_"]  = df["nHD"].map(smart_float)  if "nHD"  in df.columns else np.nan
    out["nHA_"]  = df["nHA"].map(smart_float)  if "nHA"  in df.columns else np.nan
    out["logP_"] = df["logP"].map(smart_float) if "logP" in df.columns else np.nan

    # Lipinski 3/4
    conds = [
        out.get("MW_",   pd.Series(np.nan, index=out.index))   < 500.0,
        out.get("nHD_",  pd.Series(np.nan, index=out.index))   < 5.0,
        out.get("nHA_",  pd.Series(np.nan, index=out.index))   < 10.0,
        out.get("logP_", pd.Series(np.nan, index=out.index))   < 5.0,
    ]
    out["Lipinski_passes_"] = sum([c.fillna(False) for c in conds])

    # Absorption (HIA, BA)
    out["HIA_pct_"] = df["hia"].map(to_pct) if "hia" in df.columns else np.nan
    out["Fxx_pct_"] = df.apply(pick_ba_prob, axis=1)

    # Distribution (Vd from logVDss)
    out["logVDss_"] = df["logVDss"].map(smart_float) if "logVDss" in df.columns else np.nan
    out["Vd_Lkg_"]  = np.power(10.0, out["logVDss_"]) if "logVDss" in df.columns else np.nan

    # Elimination
    out["t_half_h_"] = df["t0.5"].map(smart_float) if "t0.5" in df.columns else np.nan

    # LD50 (if available)
    out["LD50_mgkg_"]     = df["LD50_oral"].map(smart_float) if "LD50_oral" in df.columns else np.nan
    out["LD50_available_"] = out["LD50_mgkg_"].notna()

    # Safety (WoE panel)
    out["Ames_neg_"] = df["Ames"].apply(is_negative_flag) if "Ames" in df.columns else False
    out["DILI_neg_"] = df["DILI"].apply(is_negative_flag) if "DILI" in df.columns else False
    if "hERG-10um" in df.columns:
        out["hERG_neg_"] = df["hERG-10um"].apply(is_negative_flag)
    elif "hERG" in df.columns:
        out["hERG_neg_"] = df["hERG"].apply(is_negative_flag)
    else:
        out["hERG_neg_"] = False

    # Interference (all must be negative)
    inter_cols = [c for c in ("PAINS","Reactive","Aggregators","Promiscuous") if c in df.columns]
    if inter_cols:
        inter_any = np.zeros(len(df), dtype=bool)
        for c in inter_cols:
            inter_any |= df[c].apply(is_positive_flag).values
        out["Interference_any_"] = inter_any
    else:
        out["Interference_any_"] = False

    # Criteria
    out["C1_Lipinski_3of4"] = out["Lipinski_passes_"] >= 3
    out["C2_HIA_gt30"]      = out["HIA_pct_"] > 30.0
    out["C3_BA_gt30_prob"]  = out["Fxx_pct_"] >= ba_prob_cut
    out["C4_Vd_range"]      = (out["Vd_Lkg_"] >= 0.04) & (out["Vd_Lkg_"] <= 20.0)
    out["C5_t12_ge_0p5h"]   = out["t_half_h_"] >= 0.5
    out["C6_safety_panel"]  = out["Ames_neg_"] & out["DILI_neg_"] & out["hERG_neg_"] & (~out["Interference_any_"])
    out["C6_final"] = np.where(
        out["LD50_available_"],
        (out["LD50_mgkg_"] > 500.0) & out["C6_safety_panel"],
        out["C6_safety_panel"],
    )

    crit = ["C1_Lipinski_3of4","C2_HIA_gt30","C3_BA_gt30_prob","C4_Vd_range","C5_t12_ge_0p5h","C6_final"]
    out["Pass_Count"] = sum([out[c].fillna(False) for c in crit])
    out["Pass_WZ"]    = out["Pass_Count"] == 6

    # Composite score
    safety_avg = (out["Ames_neg_"].astype(float) + out["DILI_neg_"].astype(float) + out["hERG_neg_"].astype(float)) / 3.0
    interf_pen = out["Interference_any_"].astype(float)  # 1 if problematic

    def score_row(r):
        lip = (r["Lipinski_passes_"] / 4.0)
        hia = normalize_01(r["HIA_pct_"], 0.0, 100.0)
        ba  = normalize_01(r["Fxx_pct_"], 50.0, 100.0)
        vd  = normalize_logrange(r["Vd_Lkg_"], 0.04, 20.0)
        t12 = normalize_01(r["t_half_h_"], 0.5, 24.0)
        saf = r["safety_avg_"]
        pen = r["interf_pen_"]
        return (1.0*lip + 1.0*hia + 1.0*ba + 0.5*vd + 0.5*t12 + 1.0*saf - 0.5*pen)

    out["safety_avg_"]  = safety_avg
    out["interf_pen_"]  = interf_pen
    out["ADMET_score_"] = out.apply(score_row, axis=1)

    out["Hsum_"] = (out["nHD_"].fillna(99) + out["nHA_"].fillna(99))
    return out

# --- Main ---------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Rank peptides per ADMET methodology (AdmetLab 3.0 batch CSV).")
    ap.add_argument("--input", required=True, help="Path to species CSV (single header line).")
    ap.add_argument("--ba_prob_cut", type=float, default=50.0,
                    help="Probability cutoff (%) that oral BA > 30%% (default 50).")
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT,
                    help="Root folder where a subfolder named after the input file will be created.")
    args = ap.parse_args(argv)

    inp = Path(args.input).expanduser().resolve()
    outroot = Path(args.outroot).expanduser().resolve()

    df = read_csv_robust(inp)
    ranked = build_with_derivatives(df, ba_prob_cut=args.ba_prob_cut)

    # Global sort: strict pass first, then score, then tie-breakers
    ranked = ranked.sort_values(
        by=["Pass_WZ","ADMET_score_","Lipinski_passes_","MW_","Hsum_"],
        ascending=[False, False, True, True, True]
    )

    base = inp.stem
    outdir = outroot / base
    outdir.mkdir(parents=True, exist_ok=True)

    all_path        = outdir / ALL_OUT
    pass_path       = outdir / PASS_OUT
    pass_plus5_path = outdir / PASS_PLUS5_OUT
    stats_path      = outdir / STATS_OUT

    # Full ranking + strict pass (6/6)
    ranked.to_csv(all_path, index=False)
    strict_pass = ranked[ranked["Pass_WZ"]]
    strict_pass.to_csv(pass_path, index=False)

    # Main deliverable: 6/6 + (>=5/6 but !=6/6), sorted as above
    relaxed = ranked[(ranked["Pass_Count"] >= 5) & (~ranked["Pass_WZ"])]
    out_sel = pd.concat([
        strict_pass.assign(Criteria_Tier="6/6"),
        relaxed.assign(Criteria_Tier=">=5/6"),
    ], axis=0).sort_values(
        by=["Pass_WZ","ADMET_score_","Lipinski_passes_","MW_","Hsum_"],
        ascending=[False, False, True, True, True]
    )
    out_sel.to_csv(pass_plus5_path, index=False)

    # Stats / audit
    n = len(ranked)
    with stats_path.open("w", encoding="utf-8") as f:
        f.write(f"Total rows: {n}\n")
        f.write(f"Strict pass (all 6): {int(ranked['Pass_WZ'].sum())}\n")
        f.write(f"Relaxed (>=5/6): {int((ranked['Pass_Count'] >= 5).sum())}\n")
        f.write(f"C1 Lipinski ≥3/4: {int((ranked['Lipinski_passes_'] >= 3).sum())}\n")
        f.write(f"C2 HIA>30%: {int((ranked['HIA_pct_'] > 30).sum())}\n")
        f.write(f"C3 BA>30% prob ≥{args.ba_prob_cut}%: {int((ranked['Fxx_pct_'] >= args.ba_prob_cut).sum())}\n")
        f.write(f"C4 Vd in [0.04,20] L/kg: {int(((ranked['Vd_Lkg_'] >= 0.04) & (ranked['Vd_Lkg_'] <= 20)).sum())}\n")
        f.write(f"C5 t1/2 ≥0.5 h: {int((ranked['t_half_h_'] >= 0.5).sum())}\n")
        f.write("C6 safety panel (Ames-, DILI-, hERG- & no PAINS/Reactive/Aggregators/Promiscuous): "
                f"{int(ranked['C6_safety_panel'].sum())}\n")
        if "LD50_mgkg_" in ranked.columns:
            ld = ranked["LD50_mgkg_"].dropna()
            f.write(f"LD50 available: {len(ld)} / {n}\n")
            f.write(f"LD50>500 among available: {int((ld > 500).sum())}\n")
        inter_cols_present = [c for c in ("PAINS","Reactive","Aggregators","Promiscuous") if c in df.columns]
        f.write("\nInterference columns present: ")
        f.write(", ".join(inter_cols_present) if inter_cols_present else "none")
        f.write("\n")

    print(f"[OK] Out folder: {outdir}")
    print(f"[OK] Ranked ALL → {all_path}")
    print(f"[OK] Strict PASS (6/6) → {pass_path}  (rows: {len(strict_pass)})")
    print(f"[OK] PASS + >=5/6 ranked → {pass_plus5_path}  (rows: {len(out_sel)})")
    print(f"[OK] Stats → {stats_path}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())