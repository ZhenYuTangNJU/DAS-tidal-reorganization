# -*- coding: utf-8 -*-
"""Common utilities for the revised Qinhuangdao Fig.5 workflow.

Core methodological changes
---------------------------
1) Actual tidal stage is anchored to observed low/high turning points of TOTAL tide:
   Low=0°, Rising=90°, High=180°, Falling=270°, next Low=360°.
2) DAS stage dependence uses broadband log10(RMS), not the positive peak of a
   zero-mean D/SD band-passed DAS component.
3) D/SD timing is expressed relative to a Coast DAS reference. Positive phase
   offset means the local DAS modulation LAGS the Coast; negative means LEADS.

Frequency bands use the user's final Scheme A.
"""
from __future__ import annotations

import os
import re
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import detrend

# -----------------------------------------------------------------------------
# Core configuration
# -----------------------------------------------------------------------------
SITES = ["Coast", "Lagoon1", "Lagoon2", "Lagoon3"]
SITE_SHORT = {"Coast": "Coast", "Lagoon1": "L1", "Lagoon2": "L2", "Lagoon3": "L3"}
SITE_LABELS = {
    "Coast": "Coast",
    "Lagoon1": "Lagoon near station",
    "Lagoon2": "Lagoon perpendicular",
    "Lagoon3": "Lagoon near coast",
}
SITE_ORDER = ["Coast", "Lagoon3", "Lagoon2", "Lagoon1"]

# Final Scheme A: cycles per day
BANDS = OrderedDict([
    ("Diurnal", (0.65, 1.35)),
    ("Semidiurnal", (1.65, 2.35)),
    ("Terdiurnal", (2.65, 3.35)),
    ("Quarter-diurnal", (3.65, 4.35)),
])

# Preferred representative neighbourhoods already used in the Qinhuangdao workflow.
DAS_NEIGHBOR_WINDOWS = {
    "Coast": (408, 418),
    "Lagoon1": (62, 72),
    "Lagoon2": (125, 135),
    "Lagoon3": (181, 191),
}

UNIFORM_DT_MIN = 10
DAS_WINSOR_PERCENTILES = (0.5, 99.5)
BLOCK_LENGTH_M = 20.0

# -----------------------------------------------------------------------------
# Nominal geometry, copied from the existing v3 geometry workflow
# -----------------------------------------------------------------------------
CHANNEL_SPACING_M = 2.0
SITE_CHANNEL_RANGES = {
    "Lagoon1": (32, 101),
    "Lagoon2": (102, 163),
    "Lagoon3": (164, 257),
    "Coast": (330, 430),
}
PATH_BREAKS = [32, 101, 163, 257, 285, 310, 330, 430, 450]
EXACT_SEGMENT_LENGTHS_M = {
    (32, 101): 138.0,
    (101, 163): 124.0,
    (163, 257): 188.0,
    (257, 285): 56.0,
    (285, 310): 50.0,
    (310, 330): 40.0,
    (330, 430): 200.0,
    (430, 450): 40.0,
}
TOP_EDGE_AZIMUTH_DEG_FROM_EAST = 6.0
MAP_SEGMENT_AZIMUTHS_DEG = {
    (257, 285): -25.0,
    (285, 310): -136.0,
    (310, 330): -48.0,
    (330, 430): 63.0,
    (430, 450): 158.0,
}

# -----------------------------------------------------------------------------
# File selection / discovery
# -----------------------------------------------------------------------------
def choose_root_folder(title="选择秦皇岛分析总文件夹") -> Path:
    env = os.environ.get("QHD_FIG5_ROOT", "").strip()
    if env:
        p = Path(env).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"QHD_FIG5_ROOT 不存在: {p}")
        return p
    try:
        import tkinter as tk
        from tkinter import filedialog
        app = tk.Tk(); app.withdraw()
        try:
            app.attributes("-topmost", True)
        except Exception:
            pass
        folder = filedialog.askdirectory(title=title)
        app.destroy()
        if not folder:
            print("未选择文件夹，程序结束。")
            sys.exit(0)
        return Path(folder).resolve()
    except Exception:
        return Path.cwd().resolve()


def choose_excel_file(title="选择潮位 Excel", initialdir=None) -> Path:
    env = os.environ.get("QHD_TIDE_FILE", "").strip()
    if env:
        p = Path(env).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"QHD_TIDE_FILE 不存在: {p}")
        return p
    try:
        import tkinter as tk
        from tkinter import filedialog
        app = tk.Tk(); app.withdraw()
        try:
            app.attributes("-topmost", True)
        except Exception:
            pass
        p = filedialog.askopenfilename(
            title=title,
            initialdir=str(initialdir) if initialdir else None,
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
        )
        app.destroy()
        if not p:
            print("未选择文件，程序结束。")
            sys.exit(0)
        return Path(p).resolve()
    except Exception as exc:
        raise RuntimeError("无法弹出潮位文件选择窗口。") from exc


def _candidate_das_dirs(root: Path):
    out = []
    for p in [root / "DAS", root / "DAS_1min", root]:
        if p.exists() and p not in out:
            out.append(p)
    return out


def discover_das_files(root: Path):
    files = []
    for d in _candidate_das_dirs(root):
        files.extend([p for p in d.glob("*.xlsx") if not p.name.startswith("~$")])
    files = list(dict.fromkeys(p.resolve() for p in files))

    def score_file(site: str, p: Path):
        t = p.stem.lower().replace(" ", "_")
        score = 0
        if "1min" in t or "1_min" in t:
            score += 20
        if "rms" in t:
            score += 10
        if site == "Coast":
            if "coast" in t and "lagoon" not in t:
                score += 100
            else:
                return -1
        else:
            token = site.lower()
            if token in t or re.search(rf"lagoon[_-]?{site[-1]}", t):
                score += 100
            else:
                return -1
        return score

    found = {}
    for site in SITES:
        ranked = [(score_file(site, p), p) for p in files]
        ranked = [(s, p) for s, p in ranked if s >= 0]
        ranked.sort(key=lambda x: (-x[0], str(x[1])))
        if not ranked:
            raise FileNotFoundError(
                f"没有找到 {site} 的 1-min RMS Excel。建议文件名：{site}_1min_RMS.xlsx"
            )
        best_score = ranked[0][0]
        best = [p for s, p in ranked if s == best_score]
        if len(best) > 1:
            raise RuntimeError(
                f"{site} 匹配到多个同优先级文件，请只保留一个：\n" + "\n".join(str(p) for p in best)
            )
        found[site] = ranked[0][1]
    return found


def read_das_matrix(path: Path) -> pd.DataFrame:
    xls = pd.ExcelFile(path)
    sheet = "RMS_1min" if "RMS_1min" in xls.sheet_names else xls.sheet_names[0]
    raw = pd.read_excel(path, sheet_name=sheet, index_col=0)
    if raw.empty:
        raise ValueError(f"DAS Excel 为空：{path}")
    raw.index = pd.to_numeric(raw.index, errors="coerce")
    raw = raw.loc[raw.index.notna()].copy()
    raw.index = raw.index.astype(int)
    times = pd.to_datetime(raw.columns, errors="coerce")
    good = ~pd.isna(times)
    raw = raw.loc[:, good].copy()
    raw.columns = pd.DatetimeIndex(times[good])
    raw = raw.apply(pd.to_numeric, errors="coerce").sort_index(axis=1)
    return raw.T.groupby(level=0).mean().T


def read_tide_excel(path: Path) -> pd.Series:
    """Identify a datetime column and a numeric tide/water-level column."""
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
                    candidates.append((len(s), str(vc).lower(), s))
    if not candidates:
        raise RuntimeError(f"无法识别潮位 Excel 的时间列/潮位列：{path.name}")
    # Prefer names suggesting tide/water level, then length.
    def pref(item):
        n, name, _ = item
        bonus = 100000 if any(k in name for k in ["tide", "water", "level", "潮", "水位"]) else 0
        return bonus + n
    return max(candidates, key=pref)[2]

# -----------------------------------------------------------------------------
# Time series utilities
# -----------------------------------------------------------------------------
def make_uniform_time_grid(start, end, dt_min=UNIFORM_DT_MIN):
    start = pd.Timestamp(start).ceil(f"{int(dt_min)}min")
    end = pd.Timestamp(end).floor(f"{int(dt_min)}min")
    return pd.date_range(start, end, freq=f"{int(dt_min)}min")


def resample_series_to_grid(series: pd.Series, grid: pd.DatetimeIndex):
    s = series.astype(float).replace([np.inf, -np.inf], np.nan).copy()
    s = s[~s.index.duplicated(keep="first")].sort_index()
    union_idx = s.index.union(grid)
    temp = s.reindex(union_idx).sort_index()
    temp = temp.interpolate(method="time", limit_direction="both")
    return temp.reindex(grid).astype(float)


def robust_winsorize(series: pd.Series, low=0.5, high=99.5):
    s = series.astype(float).copy()
    valid = s.dropna()
    if len(valid) < 10:
        return s
    lo, hi = np.nanpercentile(valid.to_numpy(float), [low, high])
    return s.clip(lower=lo, upper=hi)


def preprocess_das_series(series: pd.Series, grid: pd.DatetimeIndex):
    """Positive RMS -> log10 -> winsor -> linear detrend -> z score.

    This is BROADBAND preprocessing. No D/SD bandpass is applied here.
    """
    s = resample_series_to_grid(series, grid)
    s[s <= 0] = np.nan
    s = np.log10(s)
    s = s.interpolate(method="time", limit_direction="both")
    s = robust_winsorize(s, *DAS_WINSOR_PERCENTILES)
    arr = detrend(s.to_numpy(float), type="linear")
    mu = np.nanmean(arr)
    sd = np.nanstd(arr, ddof=0)
    if np.isfinite(sd) and sd > 0:
        arr = (arr - mu) / sd
    else:
        arr = arr - mu
    return pd.Series(arr, index=grid)


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

# -----------------------------------------------------------------------------
# Circular utilities
# -----------------------------------------------------------------------------
def wrap_deg(x):
    return (np.asarray(x, dtype=float) + 180.0) % 360.0 - 180.0


def circular_mean_deg(phi_deg, weights=None):
    phi = np.asarray(phi_deg, dtype=float)
    m = np.isfinite(phi)
    if weights is None:
        w = np.ones_like(phi)
    else:
        w = np.asarray(weights, dtype=float)
        m &= np.isfinite(w) & (w > 0)
    if m.sum() == 0:
        return np.nan, np.nan
    a = np.deg2rad(phi[m])
    ww = w[m]
    z = np.sum(ww * np.exp(1j * a)) / np.sum(ww)
    return float(np.degrees(np.angle(z))), float(np.abs(z))

# -----------------------------------------------------------------------------
# Nominal geometry
# -----------------------------------------------------------------------------
def _rotate_xy(x, y, angle_deg):
    a = np.deg2rad(angle_deg)
    return x*np.cos(a)-y*np.sin(a), x*np.sin(a)+y*np.cos(a)


def _point_from_bearing(x0, y0, length, az_deg):
    a = np.deg2rad(az_deg)
    return float(x0 + length*np.cos(a)), float(y0 + length*np.sin(a))


def build_exact_length_layout():
    AB = EXACT_SEGMENT_LENGTHS_M[(101, 163)]
    AC = EXACT_SEGMENT_LENGTHS_M[(32, 101)]
    BC = EXACT_SEGMENT_LENGTHS_M[(163, 257)]
    A = np.array([0.0, 0.0]); B = np.array([AB, 0.0])
    cx = (AC*AC + AB*AB - BC*BC)/(2.0*AB)
    cy = -np.sqrt(max(AC*AC-cx*cx, 0.0))
    C = np.array([cx, cy])
    A = np.array(_rotate_xy(*A, TOP_EDGE_AZIMUTH_DEG_FROM_EAST))
    B = np.array(_rotate_xy(*B, TOP_EDGE_AZIMUTH_DEG_FROM_EAST))
    C = np.array(_rotate_xy(*C, TOP_EDGE_AZIMUTH_DEG_FROM_EAST))
    anchors = {101: tuple(A), 163: tuple(B), 32: tuple(C), 257: tuple(C)}
    for c0, c1 in [(257,285),(285,310),(310,330),(330,430),(430,450)]:
        x0, y0 = anchors[c0]
        anchors[c1] = _point_from_bearing(
            x0, y0, EXACT_SEGMENT_LENGTHS_M[(c0,c1)], MAP_SEGMENT_AZIMUTHS_DEG[(c0,c1)]
        )
    rows = []
    for c0, c1 in zip(PATH_BREAKS[:-1], PATH_BREAKS[1:]):
        x0,y0 = anchors[c0]; x1,y1 = anchors[c1]; n=c1-c0
        seg=[]
        for i,ch in enumerate(range(c0,c1+1)):
            a=i/float(n) if n else 0.0
            seg.append((ch, x0+a*(x1-x0), y0+a*(y1-y0)))
        if rows:
            seg=seg[1:]
        rows.extend(seg)
    g=pd.DataFrame(rows, columns=["channel","x_m","y_m"])
    def site_of(ch):
        for site,(a,b) in SITE_CHANNEL_RANGES.items():
            if a <= ch <= b:
                return site
        return "Connector"
    g["site"] = g["channel"].map(site_of)
    g["site_chainage_m"] = np.nan
    for site,(c0,c1) in SITE_CHANNEL_RANGES.items():
        m=g["site"]==site
        g.loc[m,"site_chainage_m"]=(g.loc[m,"channel"]-c0)*CHANNEL_SPACING_M
    pA=np.array(anchors[330],float); pB=np.array(anchors[430],float); v=pB-pA; vv=float(np.dot(v,v))
    dist=[]
    for _,r in g.iterrows():
        p=r[["x_m","y_m"]].to_numpy(float)
        t=float(np.dot(p-pA,v)/vv); tc=min(max(t,0.0),1.0); q=pA+tc*v
        dist.append(float(np.linalg.norm(p-q)))
    g["distance_to_coast_polyline_m"] = dist
    g["network_distance_to_CH330_m"] = np.abs(g["channel"]-330)*CHANNEL_SPACING_M
    g["block_index"] = np.floor(g["site_chainage_m"]/BLOCK_LENGTH_M).astype("Int64")
    return g
