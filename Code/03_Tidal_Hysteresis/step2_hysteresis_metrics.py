# -*- coding: utf-8 -*-
"""Step 2: matched-water-level hysteresis metrics for the revised Fig.5 (upgraded).

Upgrades relative to v1
-----------------------
1) Use fewer matched-water bins (default 14) for smoother site-level curves.
2) Export relative water level (0-1) for cleaner plotting.
3) For the spatial trend panel, use ONE representative Coast anchor at x=0,
   instead of many Coast blocks all stacked at 0 m.
4) Export delta(rising-falling) versus relative water level, for the new Fig.5d.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

# 公共计算模块随代码目录一起保留。
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))

from fig5_new_common import (
    SITES, SITE_SHORT, SITE_ORDER,
    choose_root_folder, discover_das_files, read_das_matrix,
    preprocess_das_series, build_exact_length_layout,
)

OUTDIR = "15_NEW_hysteresis_metrics_output"
PHASE_OUTDIR = "12_NEW_actual_tidal_stage_output"
PHASE_FILE = "actual_tidal_stage_DAS.xlsx"
N_LEVEL_BINS = 14
LEVEL_Q_LOW = 2.0
LEVEL_Q_HIGH = 98.0
MIN_POINTS_PER_BRANCH_BIN = 3
MIN_MATCHED_BINS = 5
MIN_CHANNELS_PER_BLOCK = 3
N_BOOT = 1000
RANDOM_SEED = 20260914


def _load_phase_time_table(root: Path):
    f = root / PHASE_OUTDIR / PHASE_FILE
    if not f.exists():
        raise FileNotFoundError(
            "找不到 Step 1 输出：%s\n请先运行 step1_actual_tidal_stage.py 或 RUN_NEW_FIG5_HYSTERESIS.py" % f
        )
    df = pd.read_excel(f, sheet_name="time_phase")
    time_col = df.columns[0]
    df[time_col] = pd.to_datetime(df[time_col], errors="coerce")
    df = df.dropna(subset=[time_col]).copy()
    df = df.rename(columns={time_col: "time"})
    df = df.set_index("time").sort_index()
    return df


def build_level_edges(level_series, n_bins=N_LEVEL_BINS, q_low=LEVEL_Q_LOW, q_high=LEVEL_Q_HIGH):
    lev = np.asarray(level_series, float)
    lev = lev[np.isfinite(lev)]
    if len(lev) < 50:
        raise RuntimeError("有效潮位点太少，无法建立 matched-water bins。")
    lo, hi = np.nanpercentile(lev, [q_low, q_high])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        raise RuntimeError("潮位分箱失败，请检查 tide_observed。")
    edges = np.linspace(lo, hi, int(n_bins) + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    rel = (centers - lo) / (hi - lo)
    return edges, centers, rel, lo, hi


def matched_bin_statistics(level, metric, halfcycle, edges, lo_level, hi_level):
    level = np.asarray(level, float)
    metric = np.asarray(metric, float)
    halfcycle = np.asarray(halfcycle, object)
    rows = []
    for i in range(len(edges) - 1):
        left = edges[i]
        right = edges[i + 1]
        if i < len(edges) - 2:
            in_bin = (level >= left) & (level < right)
        else:
            in_bin = (level >= left) & (level <= right)
        mr = in_bin & (halfcycle == "Rising") & np.isfinite(metric)
        mf = in_bin & (halfcycle == "Falling") & np.isfinite(metric)
        nr = int(np.sum(mr))
        nf = int(np.sum(mf))
        rising_med = float(np.nanmedian(metric[mr])) if nr >= MIN_POINTS_PER_BRANCH_BIN else np.nan
        falling_med = float(np.nanmedian(metric[mf])) if nf >= MIN_POINTS_PER_BRANCH_BIN else np.nan
        xc = float(0.5 * (left + right))
        xr = float((xc - lo_level) / (hi_level - lo_level)) if hi_level > lo_level else np.nan
        rows.append({
            "bin_left": float(left),
            "bin_right": float(right),
            "level_center": xc,
            "relative_level": xr,
            "n_rising": nr,
            "n_falling": nf,
            "rising_median": rising_med,
            "falling_median": falling_med,
            "delta_rf": float(rising_med - falling_med) if np.isfinite(rising_med) and np.isfinite(falling_med) else np.nan,
            "matched_bin": bool(np.isfinite(rising_med) and np.isfinite(falling_med)),
        })
    return pd.DataFrame(rows)


def hysteresis_index_from_bins(bins_df):
    d = bins_df.loc[bins_df["matched_bin"], "delta_rf"].to_numpy(float)
    d = d[np.isfinite(d)]
    if len(d) < MIN_MATCHED_BINS:
        return np.nan, int(len(d))
    return float(np.nanmean(d)), int(len(d))


def bootstrap_median(values, n_boot=N_BOOT, seed=RANDOM_SEED):
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    center = float(np.nanmedian(x))
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot, float)
    for i in range(n_boot):
        idx = rng.integers(0, len(x), size=len(x))
        boots[i] = np.nanmedian(x[idx])
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return center, float(lo), float(hi)


def bootstrap_median_matrix(matrix, n_boot=N_BOOT, seed=RANDOM_SEED):
    x = np.asarray(matrix, float)
    keep = np.isfinite(x).any(axis=1)
    x = x[keep]
    if len(x) == 0:
        nan = np.full(matrix.shape[1], np.nan)
        return nan, nan.copy(), nan.copy()
    center = np.nanmedian(x, axis=0)
    rng = np.random.default_rng(seed)
    boots = np.empty((n_boot, x.shape[1]), float)
    for i in range(n_boot):
        idx = rng.integers(0, len(x), size=len(x))
        boots[i] = np.nanmedian(x[idx], axis=0)
    lo, hi = np.nanpercentile(boots, [2.5, 97.5], axis=0)
    return center, lo, hi


def fit_linear(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if len(x) < 5:
        return None
    X = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    a, b = float(beta[0]), float(beta[1])
    yhat = a + b * x
    sst = np.sum((y - np.mean(y)) ** 2)
    ssr = np.sum((y - yhat) ** 2)
    r2 = 1.0 - ssr / sst if sst > 0 else np.nan
    return {"intercept": a, "slope_per_m": b, "slope_per_100m": b * 100.0, "R2": float(r2), "n_blocks": int(len(x))}


def _aggregate_site_curves(channel_bin_df):
    rows = []
    for isite, site in enumerate(SITES):
        sub = channel_bin_df[channel_bin_df["site"] == site].copy()
        for branch in ["rising_median", "falling_median", "delta_rf"]:
            wide = sub.pivot_table(index="channel", columns="relative_level", values=branch, aggfunc="first")
            wide = wide.reindex(sorted(wide.columns), axis=1)
            med, lo, hi = bootstrap_median_matrix(wide.to_numpy(float), seed=RANDOM_SEED + 100 * isite + (0 if branch == "rising_median" else 1 if branch == "falling_median" else 2))
            for x, y, yl, yh in zip(wide.columns.to_numpy(float), med, lo, hi):
                rows.append({
                    "site": site,
                    "site_short": SITE_SHORT[site],
                    "curve_type": branch,
                    "relative_level": float(x),
                    "median_value": float(y) if np.isfinite(y) else np.nan,
                    "lower95": float(yl) if np.isfinite(yl) else np.nan,
                    "upper95": float(yh) if np.isfinite(yh) else np.nan,
                })
    return pd.DataFrame(rows)




def run(root=None):
    root = Path(root) if root is not None else choose_root_folder("选择包含四个 DAS 1-min RMS Excel 的总文件夹")
    phase_df = _load_phase_time_table(root)
    grid = pd.DatetimeIndex(phase_df.index)

    das_files = discover_das_files(root)
    raw_by_site = {site: read_das_matrix(das_files[site]) for site in SITES}
    z_by_site = {}
    for site in SITES:
        raw = raw_by_site[site]
        z_by_site[site] = {int(ch): preprocess_das_series(raw.loc[ch], grid).to_numpy(float) for ch in sorted(raw.index.astype(int))}

    level = phase_df["tide_observed"].to_numpy(float)
    halfcycle = phase_df["halfcycle"].to_numpy(object)
    edges, centers, rel, lo_level, hi_level = build_level_edges(level)

    bb_channel_rows = []
    bb_bin_rows = []
    for site in SITES:
        for ch, z in z_by_site[site].items():
            bins = matched_bin_statistics(level, z, halfcycle, edges, lo_level, hi_level)
            H, nmatch = hysteresis_index_from_bins(bins)
            bb_channel_rows.append({"site": site, "site_short": SITE_SHORT[site], "channel": int(ch), "broadband_hysteresis_index": H, "n_matched_bins": int(nmatch)})
            tmp = bins.copy()
            tmp.insert(0, "channel", int(ch))
            tmp.insert(0, "site_short", SITE_SHORT[site])
            tmp.insert(0, "site", site)
            bb_bin_rows.append(tmp)

    broadband_channel = pd.DataFrame(bb_channel_rows)
    broadband_bins = pd.concat(bb_bin_rows, ignore_index=True) if bb_bin_rows else pd.DataFrame()
    site_curves = _aggregate_site_curves(broadband_bins)

    site_summary_rows = []
    for i, site in enumerate(SITES):
        vals = broadband_channel.loc[broadband_channel.site == site, "broadband_hysteresis_index"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        med, lo, hi = bootstrap_median(vals, seed=RANDOM_SEED + i)
        site_summary_rows.append({"site": site, "site_short": SITE_SHORT[site], "median_hysteresis_index": med, "lower95": lo, "upper95": hi, "n_channels": int(len(vals))})
    site_summary = pd.DataFrame(site_summary_rows)

    geo = build_exact_length_layout()
    geo = geo[geo.site.isin(SITES)][["site", "channel", "block_index", "distance_to_coast_polyline_m"]].copy()
    merged = broadband_channel.merge(geo, on=["site", "channel"], how="left")

    block_rows = []
    coast_vals = merged.loc[merged.site == "Coast", "broadband_hysteresis_index"].to_numpy(float)
    coast_vals = coast_vals[np.isfinite(coast_vals)]
    if len(coast_vals):
        block_rows.append({
            "site": "Coast", "site_short": "Coast", "block_index": -1,
            "distance_to_coast_polyline_m": 0.0,
            "block_hysteresis_index": float(np.nanmedian(coast_vals)),
            "n_channels": int(len(coast_vals)),
            "anchor_type": "site_anchor",
        })

    for (site, block_index), s in merged[merged.site != "Coast"].groupby(["site", "block_index"], dropna=True):
        vals = s["broadband_hysteresis_index"].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        if len(vals) < MIN_CHANNELS_PER_BLOCK:
            continue
        block_rows.append({
            "site": site,
            "site_short": SITE_SHORT[site],
            "block_index": int(block_index),
            "distance_to_coast_polyline_m": float(np.nanmedian(s.distance_to_coast_polyline_m)),
            "block_hysteresis_index": float(np.nanmedian(vals)),
            "n_channels": int(len(vals)),
            "anchor_type": "lagoon_block",
        })
    blocks = pd.DataFrame(block_rows)

    fit = fit_linear(blocks["distance_to_coast_polyline_m"], blocks["block_hysteresis_index"])
    reg_rows = [{"metric": "Broadband", **fit}] if fit is not None else []
    rng = np.random.default_rng(RANDOM_SEED)
    boot_rows = []
    lagoon_blocks = blocks[blocks.site != "Coast"].copy()
    coast_channel_vals = merged.loc[merged.site == "Coast", "broadband_hysteresis_index"].to_numpy(float)
    coast_channel_vals = coast_channel_vals[np.isfinite(coast_channel_vals)]
    for ib in range(N_BOOT):
        parts = []
        if len(coast_channel_vals):
            idx = rng.integers(0, len(coast_channel_vals), size=len(coast_channel_vals))
            parts.append(pd.DataFrame({"distance_to_coast_polyline_m": [0.0], "block_hysteresis_index": [float(np.nanmedian(coast_channel_vals[idx]))]}))
        for site in ["Lagoon3", "Lagoon2", "Lagoon1"]:
            s = lagoon_blocks[lagoon_blocks.site == site]
            if len(s):
                idx = rng.integers(0, len(s), size=len(s))
                parts.append(s.iloc[idx][["distance_to_coast_polyline_m", "block_hysteresis_index"]].reset_index(drop=True))
        if not parts:
            continue
        bdf = pd.concat(parts, ignore_index=True)
        bfit = fit_linear(bdf["distance_to_coast_polyline_m"], bdf["block_hysteresis_index"])
        if bfit is not None:
            boot_rows.append({"iteration": ib, "metric": "Broadband", **bfit})
    regression_boot = pd.DataFrame(boot_rows)
    if len(regression_boot):
        slo, smed, shi = np.nanpercentile(regression_boot["slope_per_100m"].to_numpy(float), [2.5, 50, 97.5])
    else:
        slo = smed = shi = np.nan
    regression_summary = pd.DataFrame([{
        "metric": "Broadband",
        "slope_per_100m_q025": slo,
        "slope_per_100m_median": smed,
        "slope_per_100m_q975": shi,
        "n_boot": int(len(regression_boot)),
        "trend_positive_supported": bool(np.isfinite(slo) and slo > 0),
        "trend_negative_supported": bool(np.isfinite(shi) and shi < 0),
    }])

    outdir = root / OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(outdir / "hysteresis_metrics_v2_upgraded.xlsx", engine="openpyxl") as w:
        pd.DataFrame({"bin_left": edges[:-1], "bin_right": edges[1:], "level_center": centers, "relative_level": rel}).to_excel(w, sheet_name="level_bins", index=False)
        broadband_channel.to_excel(w, sheet_name="channel_broadband_H", index=False)
        broadband_bins.to_excel(w, sheet_name="broadband_binned_curves", index=False)
        site_curves.to_excel(w, sheet_name="site_loop_curves", index=False)
        site_summary.to_excel(w, sheet_name="site_broadband_summary", index=False)
        blocks.to_excel(w, sheet_name="broadband_spatial_points", index=False)
        pd.DataFrame(reg_rows).to_excel(w, sheet_name="regressions", index=False)
        regression_summary.to_excel(w, sheet_name="regression_summary", index=False)
        regression_boot.to_excel(w, sheet_name="regression_bootstrap", index=False)
        pd.DataFrame({
            "parameter": [
                "broadband_metric", "hysteresis_definition", "positive_sign",
                "tide_for_matching", "water_level_bins", "water_level_bin_percentile_range",
                "minimum_points_per_branch_bin", "minimum_matched_bins_for_channel_H",
                "coast_spatial_treatment", "bootstrap_iterations"
            ],
            "value": [
                "Broadband log10(RMS), winsorized, detrended, within-channel z score",
                "H = mean_k [ median(I_rising, h_k) - median(I_falling, h_k) ]",
                "H > 0 means rising branch stronger at matched water level; H < 0 means falling branch stronger",
                "Observed tide (same uniform time grid as Step 1)",
                N_LEVEL_BINS, f"{LEVEL_Q_LOW}-{LEVEL_Q_HIGH} percentile",
                MIN_POINTS_PER_BRANCH_BIN, MIN_MATCHED_BINS,
                "single Coast site anchor at 0 m; lagoon represented by 20 m blocks",
                N_BOOT,
            ]
        }).to_excel(w, sheet_name="metadata", index=False)

    print("Step 2 完成：", outdir)
    print("请查看 hysteresis_metrics_v2_upgraded.xlsx。")
    return outdir / "hysteresis_metrics_v2_upgraded.xlsx"


if __name__ == "__main__":
    run()
