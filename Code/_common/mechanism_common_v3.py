#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mechanism_common_v3.py
======================

公共函数：
- 1-min DAS 读取后的规则化/预处理
- 宽频带 FFT 重建与相关
- 严格长度约束的几何构建
- STI / cosine / OLS / LOSO
- 主导方向扫描
"""

from __future__ import annotations

import sys
import math
from pathlib import Path
from typing import Dict, Tuple, Sequence, Optional

import numpy as np
import pandas as pd
from scipy.signal import detrend, periodogram

from das_rbr_1min_common import (
    SITES, SITE_SHORT, SITE_LABELS, BANDS,
    DAS_WINSOR_PERCENTILES, DAS_USE_LOG10,
    choose_root_folder, discover_das_files, read_das_matrix,
    robust_winsorize,
)
from geometry_config_v3 import (
    CHANNEL_SPACING_M, SITE_CHANNEL_RANGES, PATH_BREAKS,
    EXACT_SEGMENT_LENGTHS_M, TOP_EDGE_AZIMUTH_DEG_FROM_EAST,
    MAP_SEGMENT_AZIMUTHS_DEG, LOWER_VERTEX_GAP_M,
    REPRESENTATIVE_CHANNELS,
)

BAND_NAMES = list(BANDS.keys())
UNIFORM_DT_MIN = 10
DPI = 600
MIN_CHANNEL_COVERAGE = 0.95
MAX_SHORT_GAP_MINUTES = 5






def choose_excel_file(title="选择 Excel 文件", initialdir=None) -> Path:
    try:
        import tkinter as tk
        from tkinter import filedialog
        app = tk.Tk()
        app.withdraw()
        try:
            app.attributes("-topmost", True)
        except Exception:
            pass
        p = filedialog.askopenfilename(
            title=title,
            initialdir=str(initialdir) if initialdir else None,
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        app.destroy()
        if not p:
            print("未选择文件，程序结束。")
            sys.exit(0)
        return Path(p).resolve()
    except Exception as exc:
        raise RuntimeError("无法弹出文件选择窗口，请检查 tkinter。") from exc


# ----------------------------------------------------------------------
# Tide reading
# ----------------------------------------------------------------------
def read_tide_excel(path: Path) -> pd.Series:
    xls = pd.ExcelFile(path)
    candidates = []

    for sheet in xls.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        if df.empty:
            continue

        cols = list(df.columns)
        for tc in cols:
            t = pd.to_datetime(df[tc], errors="coerce")
            if t.notna().sum() < max(10, int(0.5 * len(df))):
                continue
            for vc in cols:
                if vc == tc:
                    continue
                v = pd.to_numeric(df[vc], errors="coerce")
                if v.notna().sum() >= max(10, int(0.5 * len(df))):
                    s = pd.Series(v.to_numpy(float), index=pd.DatetimeIndex(t))
                    s = s.loc[s.index.notna()].dropna().sort_index()
                    s = s.groupby(level=0).mean()
                    candidates.append(s)
                    break

    if not candidates:
        raise RuntimeError(f"无法识别潮位 Excel 的时间列/潮位列：{path.name}")

    return max(candidates, key=len)


# ----------------------------------------------------------------------
# DAS preprocessing for full-channel FFT
# ----------------------------------------------------------------------
def _linear_fill_short_gaps(df_txc: pd.DataFrame, limit=MAX_SHORT_GAP_MINUTES):
    return df_txc.interpolate(
        method="time",
        limit=int(limit),
        limit_direction="both",
        limit_area="inside",
    )


def preprocess_site_matrix(raw_chxt: pd.DataFrame,
                           mode="shape",
                           min_coverage=MIN_CHANNEL_COVERAGE):
    """
    Input: channel x time
    Output: times, channels, matrix(time x channel), QC

    shape:
      positive RMS -> log10 -> winsor -> detrend -> zscale
    absolute:
      positive RMS -> detrend only
    """
    raw = raw_chxt.copy().sort_index(axis=1)
    times = pd.DatetimeIndex(raw.columns)
    vals = raw.to_numpy(float)
    vals[~np.isfinite(vals)] = np.nan
    vals[vals <= 0] = np.nan

    coverage = np.isfinite(vals).mean(axis=1)
    keep = coverage >= min_coverage
    ch = raw.index.to_numpy(int)[keep]
    vals = vals[keep, :]
    if vals.size == 0:
        raise RuntimeError("没有通道达到最低覆盖率要求。")

    df = pd.DataFrame(vals.T, index=times, columns=ch)
    df = _linear_fill_short_gaps(df)
    complete = df.notna().all(axis=0).to_numpy()
    df = df.loc[:, complete]
    ch = df.columns.to_numpy(int)

    if len(ch) == 0:
        raise RuntimeError("短缺口处理后没有完整通道。")

    x = df.to_numpy(float)

    if mode == "shape":
        if DAS_USE_LOG10:
            x = np.log10(x)

        for j in range(x.shape[1]):
            lo, hi = np.nanpercentile(
                x[:, j],
                [DAS_WINSOR_PERCENTILES[0], DAS_WINSOR_PERCENTILES[1]]
            )
            x[:, j] = np.clip(x[:, j], lo, hi)

        x = detrend(x, axis=0, type="linear")
        mu = np.mean(x, axis=0)
        sd = np.std(x, axis=0, ddof=0)
        sd[sd <= 0] = np.nan
        x = (x - mu) / sd

    elif mode == "absolute":
        x = detrend(x, axis=0, type="linear")
    else:
        raise ValueError("mode 必须是 shape 或 absolute")

    qc = pd.DataFrame({
        "channel": raw.index.to_numpy(int),
        "coverage_before_fill": coverage,
        "used_in_fft": np.isin(raw.index.to_numpy(int), ch),
    })

    return times, ch, x, qc


def fft_band_metrics(matrix_txc: np.ndarray, sampling_per_day=1440.0):
    f, p = periodogram(
        matrix_txc,
        fs=float(sampling_per_day),
        window="hann",
        detrend=False,
        scaling="density",
        axis=0,
    )
    keep = (f >= 0.20) & (f <= 4.30)
    f = f[keep]
    p = p[keep, :]

    bp = np.full((matrix_txc.shape[1], len(BAND_NAMES)), np.nan)
    for bi, bn in enumerate(BAND_NAMES):
        f0, f1 = BANDS[bn]
        m = (f >= f0) & (f <= f1)
        if m.sum() >= 2:
            bp[:, bi] = np.trapz(p[m, :], f[m], axis=0)
        elif m.sum() == 1:
            dfreq = np.median(np.diff(f)) if len(f) > 1 else 1.0
            bp[:, bi] = p[m, :][0] * dfreq

    total = np.nansum(bp, axis=1, keepdims=True)
    comp = np.divide(
        bp, total,
        out=np.full_like(bp, np.nan),
        where=total > 0,
    )
    return f, p, bp, comp


# ----------------------------------------------------------------------
# Time grid / band-limited correlation
# ----------------------------------------------------------------------
def make_uniform_time_grid(start, end, dt_min=UNIFORM_DT_MIN):
    start = pd.Timestamp(start).ceil(f"{int(dt_min)}min")
    end = pd.Timestamp(end).floor(f"{int(dt_min)}min")
    return pd.date_range(start, end, freq=f"{int(dt_min)}min")


def resample_series_to_grid(series: pd.Series, grid: pd.DatetimeIndex):
    s = series.copy()
    s = s[~s.index.duplicated(keep="first")].sort_index()
    union_idx = s.index.union(grid)
    temp = s.reindex(union_idx).sort_index()
    temp = temp.interpolate(method="time", limit_direction="both")
    return temp.reindex(grid).astype(float)


def preprocess_signal_series(s: pd.Series,
                             apply_log=True,
                             winsor=True,
                             do_detrend=True,
                             do_zscore=True):
    x = s.astype(float).replace([np.inf, -np.inf], np.nan).copy()
    if apply_log:
        x[x <= 0] = np.nan
        x = np.log10(x)

    x = x.interpolate(method="time", limit_direction="both")

    if winsor:
        x = robust_winsorize(x, *DAS_WINSOR_PERCENTILES)

    arr = x.to_numpy(float)
    if do_detrend:
        arr = detrend(arr, type="linear")

    if do_zscore:
        sd = np.nanstd(arr, ddof=0)
        mu = np.nanmean(arr)
        arr = arr - mu
        if np.isfinite(sd) and sd > 0:
            arr = arr / sd

    return pd.Series(arr, index=x.index)


def fft_bandpass_uniform(y, dt_min, band):
    y = np.asarray(y, dtype=float)
    m = np.isfinite(y)
    if m.sum() < 8:
        return np.full_like(y, np.nan)

    if not m.all():
        idx = np.arange(len(y))
        y = np.interp(idx, idx[m], y[m])

    y = y - np.mean(y)
    n = len(y)
    dt_day = float(dt_min) / 1440.0
    freqs = np.fft.rfftfreq(n, d=dt_day)
    Y = np.fft.rfft(y)
    f0, f1 = band
    mask = (freqs >= f0) & (freqs <= f1)
    Y2 = np.zeros_like(Y)
    Y2[mask] = Y[mask]
    return np.fft.irfft(Y2, n=n)


def band_variance_fraction(y, dt_min, band):
    y = np.asarray(y, dtype=float)
    den = np.nanvar(y)
    if not np.isfinite(den) or den <= 0:
        return np.nan
    yb = fft_bandpass_uniform(y, dt_min, band)
    return float(np.nanvar(yb) / den)


def pearsonr_safe(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 8:
        return np.nan
    x = x[m]
    y = y[m]
    sx = np.std(x, ddof=0)
    sy = np.std(y, ddof=0)
    if sx <= 0 or sy <= 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def max_corr_with_lag(x, y, dt_min, max_lag_hours):
    """
    Positive lag means y lags x.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return np.nan, np.nan

    x = x[m]
    y = y[m]
    max_steps = int(round(max_lag_hours * 60.0 / dt_min))

    best_r = -np.inf
    best_lag = np.nan

    for k in range(-max_steps, max_steps + 1):
        if k < 0:
            xx = x[-k:]
            yy = y[:len(y)+k]
        elif k > 0:
            xx = x[:-k]
            yy = y[k:]
        else:
            xx = x
            yy = y

        if len(xx) < 8:
            continue
        r = pearsonr_safe(xx, yy)
        if np.isfinite(r) and r > best_r:
            best_r = r
            best_lag = k * dt_min / 60.0

    return float(best_r), float(best_lag)


# ----------------------------------------------------------------------
# Exact-length geometry
# ----------------------------------------------------------------------
def rotate_xy(x, y, angle_deg):
    a = np.deg2rad(angle_deg)
    xr = x * np.cos(a) - y * np.sin(a)
    yr = x * np.sin(a) + y * np.cos(a)
    return xr, yr


def point_from_bearing(x0, y0, length, az_deg_from_east):
    a = np.deg2rad(az_deg_from_east)
    return (
        float(x0 + length * np.cos(a)),
        float(y0 + length * np.sin(a))
    )


def triangle_internal_angles_deg(a, b, c):
    """
    sides a,b,c opposite A,B,C
    """
    A = np.degrees(np.arccos(np.clip((b*b + c*c - a*a)/(2*b*c), -1, 1)))
    B = np.degrees(np.arccos(np.clip((a*a + c*c - b*b)/(2*a*c), -1, 1)))
    C = 180.0 - A - B
    return A, B, C


def build_exact_length_layout(top_azimuth_deg=TOP_EDGE_AZIMUTH_DEG_FROM_EAST,
                              segment_azimuths=None,
                              lower_vertex_gap_m=LOWER_VERTEX_GAP_M):
    """
    Build nominal geometry:
    - triangle side lengths exact
    - overall/map bearings approximate
    """
    if segment_azimuths is None:
        segment_azimuths = dict(MAP_SEGMENT_AZIMUTHS_DEG)

    # local triangle before rotation
    # A = CH101, B = CH163, C = lower vertex CH32/257
    AB = EXACT_SEGMENT_LENGTHS_M[(101, 163)]  # 124
    AC = EXACT_SEGMENT_LENGTHS_M[(32, 101)]   # 138
    BC = EXACT_SEGMENT_LENGTHS_M[(163, 257)]  # 188

    A = np.array([0.0, 0.0])
    B = np.array([AB, 0.0])

    cx = (AC*AC + AB*AB - BC*BC) / (2.0 * AB)
    cy = -np.sqrt(max(AC*AC - cx*cx, 0.0))
    C = np.array([cx, cy])

    # rotate to map-like orientation
    def rr(P):
        x, y = rotate_xy(P[0], P[1], top_azimuth_deg)
        return np.array([x, y], dtype=float)

    A = rr(A)
    B = rr(B)
    C = rr(C)

    anchors = {
        101: tuple(A),
        163: tuple(B),
        32: tuple(C),
        257: tuple(C),
    }

    # optional tiny lower-vertex gap, if user later wants it
    if lower_vertex_gap_m != 0:
        # shift CH257 along first connector direction, keeping CH32 at triangle vertex
        az = segment_azimuths[(257, 285)]
        anchors[257] = point_from_bearing(C[0], C[1], lower_vertex_gap_m, az)

    # map-derived connector/coast directions with exact lengths
    for c0, c1 in [(257, 285), (285, 310), (310, 330), (330, 430), (430, 450)]:
        x0, y0 = anchors[c0]
        L = EXACT_SEGMENT_LENGTHS_M[(c0, c1)]
        anchors[c1] = point_from_bearing(x0, y0, L, segment_azimuths[(c0, c1)])

    # interpolate each path segment using exact channel count
    rows = []
    for c0, c1 in zip(PATH_BREAKS[:-1], PATH_BREAKS[1:]):
        x0, y0 = anchors[c0]
        x1, y1 = anchors[c1]
        n = c1 - c0
        seg = []
        for i, ch in enumerate(range(c0, c1 + 1)):
            a = i / float(n) if n > 0 else 0.0
            x = x0 + a * (x1 - x0)
            y = y0 + a * (y1 - y0)
            seg.append((ch, x, y))
        if rows:
            seg = seg[1:]
        rows.extend(seg)

    g = pd.DataFrame(rows, columns=["channel", "x_m", "y_m"])

    def site_of(ch):
        for site, (a, b) in SITE_CHANNEL_RANGES.items():
            if a <= ch <= b:
                return site
        return "Connector"

    g["site"] = g["channel"].map(site_of)
    g["site_label"] = g["site"].map(lambda s: SITE_LABELS.get(s, "Connector"))
    g["is_analysis_site"] = g["site"].isin(SITES)
    g["inside_lagoon"] = g["site"].isin(["Lagoon1", "Lagoon2", "Lagoon3"]).astype(int)

    # local chainage per site
    g["site_chainage_m"] = np.nan
    for site, (c0, c1) in SITE_CHANNEL_RANGES.items():
        m = g["site"] == site
        g.loc[m, "site_chainage_m"] = (g.loc[m, "channel"] - c0) * CHANNEL_SPACING_M

    # exact cable-network distance to CH330 from channel numbering
    g["network_distance_to_CH330_m"] = np.abs(g["channel"] - 330) * CHANNEL_SPACING_M

    # local fiber azimuth from coordinates
    azs = np.full(len(g), np.nan)
    for i in range(len(g)):
        if i == 0:
            p0 = g.loc[i, ["x_m", "y_m"]].to_numpy(float)
            p1 = g.loc[i+1, ["x_m", "y_m"]].to_numpy(float)
        elif i == len(g)-1:
            p0 = g.loc[i-1, ["x_m", "y_m"]].to_numpy(float)
            p1 = g.loc[i, ["x_m", "y_m"]].to_numpy(float)
        else:
            p0 = g.loc[i-1, ["x_m", "y_m"]].to_numpy(float)
            p1 = g.loc[i+1, ["x_m", "y_m"]].to_numpy(float)
        d = p1 - p0
        azs[i] = (np.degrees(np.arctan2(d[1], d[0])) + 360.0) % 180.0
    g["fiber_azimuth_deg"] = azs

    # Coast line segment exact from CH330 to CH430
    pA = np.array(anchors[330], dtype=float)
    pB = np.array(anchors[430], dtype=float)
    v = pB - pA
    vv = float(np.dot(v, v))

    coast_az = (np.degrees(np.arctan2(v[1], v[0])) + 360.0) % 180.0
    coast_normal_az = (coast_az + 90.0) % 180.0

    d_to_line = []
    along = []
    signed_cross = []
    nearest_x = []
    nearest_y = []

    for _, r in g.iterrows():
        p = r[["x_m", "y_m"]].to_numpy(float)
        t = float(np.dot(p - pA, v) / vv)
        tc = min(max(t, 0.0), 1.0)
        q = pA + tc * v
        dvec = p - q
        d_to_line.append(float(np.linalg.norm(dvec)))
        nearest_x.append(float(q[0]))
        nearest_y.append(float(q[1]))
        along.append(float(t * np.sqrt(vv)))
        crossz = v[0] * dvec[1] - v[1] * dvec[0]
        signed_cross.append(float(np.sign(crossz) * np.linalg.norm(dvec)))

    g["distance_to_coast_polyline_m"] = d_to_line
    g["nearest_coast_x_m"] = nearest_x
    g["nearest_coast_y_m"] = nearest_y
    g["coast_alongshore_m"] = along
    g["coast_crossshore_signed_m"] = signed_cross
    g["coast_azimuth_deg"] = coast_az
    g["coast_normal_azimuth_deg"] = coast_normal_az

    # angle to coast axis / coast normal
    def acute_axis_diff(a, b):
        return np.abs((np.asarray(a) - float(b) + 90.0) % 180.0 - 90.0)

    g["acute_angle_to_coast_axis_deg"] = acute_axis_diff(g["fiber_azimuth_deg"], coast_az)
    g["acute_angle_to_coast_normal_deg"] = acute_axis_diff(g["fiber_azimuth_deg"], coast_normal_az)
    g["orientation_cos2_to_coast_axis"] = np.cos(np.deg2rad(g["acute_angle_to_coast_axis_deg"]))**2
    g["orientation_cos2_to_coast_normal"] = np.cos(np.deg2rad(g["acute_angle_to_coast_normal_deg"]))**2

    # distance to representative Coast channel
    ref = g[g["channel"] == REPRESENTATIVE_CHANNELS["Coast"]].iloc[0]
    g["distance_to_Coast413_m"] = np.hypot(g["x_m"] - ref["x_m"], g["y_m"] - ref["y_m"])

    # exact triangle angle metadata
    # angle at CH101 is opposite BC=188
    angle_101 = np.degrees(np.arccos(np.clip((AC*AC + AB*AB - BC*BC)/(2*AC*AB), -1, 1)))
    angle_163 = np.degrees(np.arccos(np.clip((AB*AB + BC*BC - AC*AC)/(2*AB*BC), -1, 1)))
    angle_lower = 180.0 - angle_101 - angle_163

    meta = {
        "triangle_side_CH32_101_m": AC,
        "triangle_side_CH101_163_m": AB,
        "triangle_side_CH163_257_m": BC,
        "triangle_angle_at_CH101_deg": angle_101,
        "triangle_angle_at_CH163_deg": angle_163,
        "triangle_angle_at_lower_vertex_deg": angle_lower,
        "top_edge_azimuth_deg_from_east": top_azimuth_deg,
        "coast_azimuth_deg_from_east": coast_az,
        "coast_normal_azimuth_deg_from_east": coast_normal_az,
        "lower_vertex_gap_m": lower_vertex_gap_m,
    }

    return g, pd.DataFrame({
        "parameter": list(meta.keys()),
        "value": list(meta.values()),
    })


# ----------------------------------------------------------------------
# STI / similarity / models
# ----------------------------------------------------------------------
def spectral_transformation_index(v, ref):
    v = np.asarray(v, float)
    ref = np.asarray(ref, float)
    if not np.isfinite(v).all() or not np.isfinite(ref).all():
        return np.nan
    return float(0.5 * np.sum(np.abs(v - ref)))


def cosine_similarity(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return np.nan
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na <= 0 or nb <= 0:
        return np.nan
    return float(np.dot(a, b) / (na * nb))


def logit(p, eps=1e-6):
    p = np.asarray(p, float)
    p = np.clip(p, eps, 1.0-eps)
    return np.log(p / (1.0-p))


def fisher_z(r, eps=1e-6):
    r = np.asarray(r, float)
    r = np.clip(r, -1.0+eps, 1.0-eps)
    return np.arctanh(r)


def ols_fit(y, X, xnames):
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]

    valid = np.isfinite(y)
    if X.size:
        valid &= np.isfinite(X).all(axis=1)

    y = y[valid]
    X = X[valid] if X.size else np.empty((len(y), 0))

    n = len(y)
    if n < max(8, X.shape[1] + 3):
        return None

    A = np.column_stack([np.ones(n), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ beta
    resid = y - pred
    rss = float(np.sum(resid**2))
    tss = float(np.sum((y-y.mean())**2))
    r2 = 1.0-rss/tss if tss > 0 else np.nan
    k = A.shape[1]
    sigma2 = rss / n if n else np.nan
    aic = n*np.log(sigma2) + 2*k if sigma2 > 0 else np.nan
    aicc = aic + 2*k*(k+1)/(n-k-1) if (n-k-1) > 0 and np.isfinite(aic) else np.nan

    return {
        "n": n,
        "beta": beta,
        "coef_names": ["Intercept"] + list(xnames),
        "r2": float(r2),
        "aicc": float(aicc) if np.isfinite(aicc) else np.nan,
        "rss": rss,
        "pred": pred,
        "resid": resid,
    }


def leave_one_site_out_rmse(df, ycol, xcols, sitecol="site"):
    sq = []
    ntest = 0
    for held in df[sitecol].dropna().unique():
        tr = df[df[sitecol] != held].dropna(subset=[ycol] + list(xcols))
        te = df[df[sitecol] == held].dropna(subset=[ycol] + list(xcols))
        if len(tr) < max(8, len(xcols)+3) or len(te) == 0:
            continue
        fit = ols_fit(tr[ycol].to_numpy(float), tr[list(xcols)].to_numpy(float), list(xcols))
        if fit is None:
            continue
        A = np.column_stack([np.ones(len(te)), te[list(xcols)].to_numpy(float)])
        pred = A @ fit["beta"]
        err = te[ycol].to_numpy(float) - pred
        sq.extend((err**2).tolist())
        ntest += len(te)
    if not sq:
        return np.nan, 0
    return float(np.sqrt(np.mean(sq))), int(ntest)


def direction_factor(fiber_az_deg, phi_deg, eps=1e-4):
    ori = np.cos(np.deg2rad(np.asarray(fiber_az_deg, float)-float(phi_deg)))**2
    return np.log(np.clip(ori, eps, 1.0))
