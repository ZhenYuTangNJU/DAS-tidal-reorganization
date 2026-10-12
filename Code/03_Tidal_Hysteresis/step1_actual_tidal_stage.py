# -*- coding: utf-8 -*-
"""Step 1: observed high/low anchored tidal stage versus broadband DAS intensity.

Outputs:
12_NEW_actual_tidal_stage_output/
  actual_tidal_stage_DAS.xlsx
  潮位转折点数据保存在Excel中，本处理版不生成诊断图

IMPORTANT
---------
The plotted DAS metric is broadband standardized log10(RMS). It is NOT a D/SD
zero-mean band-passed component. Therefore the stage curve answers: at which
OBSERVED tidal stage is broadband DAS vibration relatively stronger?
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter

# 公共计算模块随代码目录一起保留。
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))

from fig5_new_common import (
    SITES, SITE_SHORT, SITE_ORDER, UNIFORM_DT_MIN,
    choose_root_folder, choose_excel_file, discover_das_files, read_das_matrix,
    read_tide_excel, make_uniform_time_grid, resample_series_to_grid,
    preprocess_das_series,
)

OUTDIR = "12_NEW_actual_tidal_stage_output"
PHASE_GRID = np.arange(0.0, 360.0 + 0.1, 5.0)
N_BOOT = 1000
RANDOM_SEED = 20260914

# Turning-point detection controls. These are deliberately conservative and
# are exported to metadata for traceability.
SMOOTH_WINDOW_MIN = 70       # ~1.2 h, suppresses small high-frequency wiggles
MIN_EXTREMA_DISTANCE_H = 3.0 # neighbouring same-type extrema cannot be too close
PROMINENCE_FRACTION = 0.08   # fraction of q95-q05 tide range

STAGES = ["Low", "Rising", "High", "Falling"]


def _odd_window(n):
    n = max(5, int(n))
    return n if n % 2 == 1 else n + 1


def detect_turning_points(tide_u: pd.Series):
    y = tide_u.to_numpy(float)
    dt = UNIFORM_DT_MIN
    win = _odd_window(round(SMOOTH_WINDOW_MIN / dt))
    if win >= len(y):
        win = _odd_window(max(5, len(y)//5))
    if win >= len(y):
        win = len(y)-1 if len(y)%2==0 else len(y)
    smooth = savgol_filter(y, window_length=win, polyorder=2, mode="interp")

    q05, q95 = np.nanpercentile(smooth, [5, 95])
    prominence = max((q95-q05)*PROMINENCE_FRACTION, np.finfo(float).eps)
    distance = max(1, int(round(MIN_EXTREMA_DISTANCE_H*60.0/dt)))

    hi, _ = find_peaks(smooth, distance=distance, prominence=prominence)
    lo, _ = find_peaks(-smooth, distance=distance, prominence=prominence)

    events = []
    for i in hi:
        events.append({"idx": int(i), "type": "High", "level": float(smooth[i])})
    for i in lo:
        events.append({"idx": int(i), "type": "Low", "level": float(smooth[i])})
    events.sort(key=lambda r: r["idx"])

    # Enforce alternation. If two neighbouring events have the same type,
    # keep the more extreme one.
    clean = []
    for e in events:
        if not clean or clean[-1]["type"] != e["type"]:
            clean.append(e)
        else:
            prev = clean[-1]
            take_new = (e["level"] > prev["level"]) if e["type"] == "High" else (e["level"] < prev["level"])
            if take_new:
                clean[-1] = e

    # Remove pathological adjacent events < 1.5 h apart after alternation.
    min_halfcycle_steps = int(round(1.5*60.0/dt))
    changed = True
    while changed and len(clean) >= 3:
        changed = False
        for j in range(len(clean)-1):
            if clean[j+1]["idx"] - clean[j]["idx"] < min_halfcycle_steps:
                # Keep the event that gives the larger contrast to its other neighbour.
                if j == 0:
                    del clean[j]
                elif j+1 == len(clean)-1:
                    del clean[j+1]
                else:
                    left_contrast = abs(clean[j]["level"] - clean[j-1]["level"])
                    right_contrast = abs(clean[j+2]["level"] - clean[j+1]["level"])
                    if left_contrast >= right_contrast:
                        del clean[j+1]
                    else:
                        del clean[j]
                changed = True
                break

    ev = pd.DataFrame(clean)
    if len(ev) < 4:
        raise RuntimeError("识别到的高潮/低潮转折点太少，请检查潮位文件或调整检测参数。")
    ev["time"] = tide_u.index[ev["idx"].to_numpy(int)]
    ev["observed_level"] = tide_u.iloc[ev["idx"].to_numpy(int)].to_numpy(float)
    ev["prominence_threshold"] = prominence
    ev["smooth_window_min"] = SMOOTH_WINDOW_MIN
    ev["min_extrema_distance_h"] = MIN_EXTREMA_DISTANCE_H
    return pd.Series(smooth, index=tide_u.index), ev


def build_actual_tidal_phase(index: pd.DatetimeIndex, events: pd.DataFrame):
    phase = np.full(len(index), np.nan, float)
    cycle_id = np.full(len(index), np.nan, float)
    halfcycle = np.full(len(index), "", dtype=object)

    ev = events.sort_values("idx").reset_index(drop=True)
    cid = 0
    for k in range(len(ev)-1):
        a = ev.iloc[k]; b = ev.iloc[k+1]
        i0, i1 = int(a.idx), int(b.idx)
        if i1 <= i0:
            continue
        if a.type == "Low" and b.type == "High":
            vals = np.linspace(0.0, 180.0, i1-i0+1)
            label = "Rising"
            cid += 1
        elif a.type == "High" and b.type == "Low":
            vals = np.linspace(180.0, 360.0, i1-i0+1)
            label = "Falling"
        else:
            continue
        phase[i0:i1+1] = vals
        cycle_id[i0:i1+1] = cid
        halfcycle[i0:i1+1] = label

    # 360° and 0° are the same low-water phase; wrap only finite values
    # to [0, 360). Values outside complete turning-point intervals remain NaN.
    finite = np.isfinite(phase)
    phase[finite] = np.mod(phase[finite], 360.0)
    return pd.DataFrame({"phase_deg": phase, "cycle_id": cycle_id, "halfcycle": halfcycle}, index=index)


def stage_from_phase(phi):
    phi = np.asarray(phi, float).copy()
    finite = np.isfinite(phi)
    phi[finite] = np.mod(phi[finite], 360.0)
    out = np.full(phi.shape, "", dtype=object)
    out[(phi < 45) | (phi >= 315)] = "Low"
    out[(phi >= 45) & (phi < 135)] = "Rising"
    out[(phi >= 135) & (phi < 225)] = "High"
    out[(phi >= 225) & (phi < 315)] = "Falling"
    return out


def harmonic_prediction(phase_deg, y, phase_grid=PHASE_GRID):
    p = np.asarray(phase_deg, float)
    yy = np.asarray(y, float)
    m = np.isfinite(p) & np.isfinite(yy)
    if m.sum() < 30:
        return np.full(len(phase_grid), np.nan), np.nan
    th = np.deg2rad(p[m])
    X = np.column_stack([
        np.ones(m.sum()), np.cos(th), np.sin(th), np.cos(2*th), np.sin(2*th)
    ])
    beta, *_ = np.linalg.lstsq(X, yy[m], rcond=None)
    pred_obs = X @ beta
    sst = np.sum((yy[m]-np.mean(yy[m]))**2)
    ssr = np.sum((yy[m]-pred_obs)**2)
    r2 = 1.0 - ssr/sst if sst > 0 else np.nan
    tg = np.deg2rad(phase_grid)
    Xg = np.column_stack([
        np.ones(len(tg)), np.cos(tg), np.sin(tg), np.cos(2*tg), np.sin(2*tg)
    ])
    return Xg @ beta, float(r2)


def bootstrap_median_matrix(matrix, n_boot=N_BOOT, seed=RANDOM_SEED):
    x = np.asarray(matrix, float)
    x = x[np.isfinite(x).any(axis=1)]
    if len(x) == 0:
        ncol = matrix.shape[1]
        nan = np.full(ncol, np.nan)
        return nan, nan.copy(), nan.copy()
    med = np.nanmedian(x, axis=0)
    rng = np.random.default_rng(seed)
    boots = np.empty((n_boot, x.shape[1]), float)
    for i in range(n_boot):
        idx = rng.integers(0, len(x), size=len(x))
        boots[i] = np.nanmedian(x[idx], axis=0)
    lo, hi = np.nanpercentile(boots, [2.5, 97.5], axis=0)
    return med, lo, hi


def run(root=None, tide_path=None):
    root = Path(root) if root is not None else choose_root_folder("选择包含四个 DAS 1-min RMS Excel 的总文件夹")
    tide_path = Path(tide_path) if tide_path is not None else choose_excel_file("选择 2025 DAS 同期潮位 Excel", root)

    das_files = discover_das_files(root)
    raw_by_site = {site: read_das_matrix(das_files[site]) for site in SITES}
    tide = read_tide_excel(tide_path)

    starts = [raw.columns.min() for raw in raw_by_site.values()]
    ends = [raw.columns.max() for raw in raw_by_site.values()]
    common_start = max(max(starts), tide.index.min())
    common_end = min(min(ends), tide.index.max())
    grid = make_uniform_time_grid(common_start, common_end, UNIFORM_DT_MIN)
    if len(grid) < 100:
        raise RuntimeError("共同 DAS–tide 时间窗太短。")

    tide_u = resample_series_to_grid(tide, grid)
    tide_smooth, events = detect_turning_points(tide_u)
    phase_df = build_actual_tidal_phase(grid, events)
    phase_df["stage"] = stage_from_phase(phase_df["phase_deg"].to_numpy(float))
    phase_df["tide_observed"] = tide_u.to_numpy(float)
    phase_df["tide_smooth"] = tide_smooth.to_numpy(float)

    pred_rows = []
    channel_stage_rows = []
    channel_qc_rows = []

    for site in SITES:
        raw = raw_by_site[site]
        for ch in sorted(raw.index.astype(int)):
            z = preprocess_das_series(raw.loc[ch], grid).to_numpy(float)
            pred, r2 = harmonic_prediction(phase_df["phase_deg"].to_numpy(float), z)
            if np.isfinite(pred).any():
                pmax = float(PHASE_GRID[np.nanargmax(pred)] % 360.0)
            else:
                pmax = np.nan
            channel_qc_rows.append({
                "site": site, "site_short": SITE_SHORT[site], "channel": int(ch),
                "harmonic_R2": r2, "secondary_fitted_peak_phase_deg": pmax,
            })
            for ph, val in zip(PHASE_GRID, pred):
                pred_rows.append({
                    "site": site, "site_short": SITE_SHORT[site], "channel": int(ch),
                    "phase_deg": float(ph), "fitted_broadband_logRMS_z": float(val) if np.isfinite(val) else np.nan,
                })
            stages = phase_df["stage"].to_numpy(object)
            for stage in STAGES:
                m = (stages == stage) & np.isfinite(z)
                channel_stage_rows.append({
                    "site": site, "site_short": SITE_SHORT[site], "channel": int(ch),
                    "stage": stage, "mean_broadband_logRMS_z": float(np.nanmean(z[m])) if m.sum() else np.nan,
                    "n_time_points": int(m.sum()),
                })

    pred_df = pd.DataFrame(pred_rows)
    ch_stage = pd.DataFrame(channel_stage_rows)
    ch_qc = pd.DataFrame(channel_qc_rows)

    # Site-level curve: median across channels + channel bootstrap CI.
    site_curve_rows = []
    for si, site in enumerate(SITES):
        sub = pred_df[pred_df.site == site]
        wide = sub.pivot(index="channel", columns="phase_deg", values="fitted_broadband_logRMS_z").reindex(columns=PHASE_GRID)
        med, lo, hi = bootstrap_median_matrix(wide.to_numpy(float), seed=RANDOM_SEED+si)
        for ph, y, yl, yh in zip(PHASE_GRID, med, lo, hi):
            site_curve_rows.append({
                "site": site, "site_short": SITE_SHORT[site], "phase_deg": float(ph),
                "median_broadband_logRMS_z": float(y),
                "lower95": float(yl), "upper95": float(yh),
            })
    site_curves = pd.DataFrame(site_curve_rows)

    # Site-stage summary with channel bootstrap CI.
    summary_rows = []
    rng = np.random.default_rng(RANDOM_SEED)
    for site in SITES:
        for stage in STAGES:
            vals = ch_stage[(ch_stage.site==site)&(ch_stage.stage==stage)]["mean_broadband_logRMS_z"].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if len(vals):
                center = float(np.nanmedian(vals))
                boots = np.empty(N_BOOT, float)
                for i in range(N_BOOT):
                    boots[i] = np.nanmedian(vals[rng.integers(0,len(vals),size=len(vals))])
                lo, hi = np.percentile(boots, [2.5,97.5])
            else:
                center=lo=hi=np.nan
            summary_rows.append({
                "site": site, "site_short": SITE_SHORT[site], "stage": stage,
                "median_channel_stage_intensity_z": center,
                "lower95": float(lo), "upper95": float(hi), "n_channels": len(vals),
            })
    site_stage = pd.DataFrame(summary_rows)

    outdir = root / OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)

    events_out = events.copy()
    events_out["time"] = pd.to_datetime(events_out["time"])
    with pd.ExcelWriter(outdir / "actual_tidal_stage_DAS.xlsx", engine="openpyxl") as w:
        events_out.to_excel(w, sheet_name="turning_points", index=False)
        # pandas < 1.5 does not support reset_index(names=...).
        phase_export = phase_df.reset_index()
        phase_export = phase_export.rename(columns={phase_export.columns[0]: "time"})
        phase_export.to_excel(w, sheet_name="time_phase", index=False)
        pred_df.to_excel(w, sheet_name="channel_fitted_curves", index=False)
        site_curves.to_excel(w, sheet_name="site_phase_curves", index=False)
        ch_stage.to_excel(w, sheet_name="channel_stage_means", index=False)
        site_stage.to_excel(w, sheet_name="site_stage_summary", index=False)
        ch_qc.to_excel(w, sheet_name="channel_fit_QC", index=False)
        pd.DataFrame({
            "parameter": [
                "phase_definition", "DAS_metric", "DAS_bandpass_for_stage_curve",
                "harmonic_model", "stage_quadrants", "time_grid_min",
                "smooth_window_min", "min_extrema_distance_h", "prominence_fraction",
                "bootstrap_unit", "bootstrap_iterations",
            ],
            "value": [
                "Observed total tide: Low=0, Rising=90, High=180, Falling=270, next Low=360 deg",
                "Broadband log10(RMS), 0.5-99.5% winsorized, linearly detrended, within-channel z score",
                "NONE",
                "beta0 + beta1 cos(theta)+beta2 sin(theta)+beta3 cos(2theta)+beta4 sin(2theta)",
                "Low:315-45; Rising:45-135; High:135-225; Falling:225-315 deg",
                UNIFORM_DT_MIN, SMOOTH_WINDOW_MIN, MIN_EXTREMA_DISTANCE_H, PROMINENCE_FRACTION,
                "DAS channels", N_BOOT,
            ]
        }).to_excel(w, sheet_name="metadata", index=False)


    print("Step 1 完成：", outdir)
    print("请核对 actual_tidal_stage_DAS.xlsx 中的潮位转折点和阶段判定。")
    return outdir / "actual_tidal_stage_DAS.xlsx"


if __name__ == "__main__":
    run()
