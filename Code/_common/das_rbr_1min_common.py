#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
秦皇岛 DAS–RBR 1-min 分析公共模块
=================================

用途
----
供同目录下的 01–05 脚本共同调用，保证所有图使用完全一致的：
- 新版 DAS 1-min RMS 输入；
- RBR .rsk 读取、气压校正与收放仪异常清理；
- 潮汐频带定义；
- 频谱、频带积分、cosine / MAE 定义。

DAS 输入
--------
由 Qinhuangdao_DAS_raw_to_1min_RMS.py 生成：
    Coast_1min_RMS.xlsx
    Lagoon1_1min_RMS.xlsx
    Lagoon2_1min_RMS.xlsx
    Lagoon3_1min_RMS.xlsx

每个 Excel：
- sheet: RMS_1min
- 第一列 channel
- 后续列为 1-min 时间
- 单元格为由 200-Hz 原始 DAS 计算的 demeaned RMS

重要方法学说明
--------------
1. 1-min RMS 是基础产品，不再回到旧 30-min RMS 文件。
2. 潮汐频带最高仅 4.25 cycles/day。1-min 采样对此严重过采样。
   为降低 Lomb–Scargle 计算量，默认直接使用完整 1-min 序列进行频谱计算，不再强制降到 30 min。
3. 如果 Figure 2 / Figure 4B 的邻近通道 Lomb–Scargle 在你的电脑上运行过慢，
   可把 ANALYSIS_RESAMPLE_MINUTES 从 1 改成 5，仅在频谱计算阶段临时聚合；
   原 1-min Excel 不会被修改。
4. 频带组成分析：DAS 使用 log10(RMS) + 0.5–99.5% winsorization，
   再 detrend + 标准化后计算 Lomb–Scargle；这与旧版“比较谱形”逻辑一致。
5. 空间绝对格局（Figure 4B）不使用逐序列标准化；使用原始 RMS
   （仅 detrend）保留四位置间的相对幅值信息。
"""

from __future__ import annotations

import math
import re
import sqlite3
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.signal import lombscargle

# -----------------------------------------------------------------------------
# 基础配置
# -----------------------------------------------------------------------------
SITES = ["Coast", "Lagoon1", "Lagoon2", "Lagoon3"]

SITE_LABELS = {
    "Coast": "Coast",
    "Lagoon1": "Lagoon near station",
    "Lagoon2": "Lagoon perpendicular",
    "Lagoon3": "Lagoon near coast",
}

SITE_FIBER_LABEL = {
    "Coast": "Coast fiber",
    "Lagoon1": "Lagoon near station fiber",
    "Lagoon2": "Lagoon perpendicular fiber",
    "Lagoon3": "Lagoon near coast fiber",
}

SITE_COMPACT_LABEL = {
    "Coast": "Coast",
    "Lagoon1": "Lagoon near\nstation",
    "Lagoon2": "Lagoon\nperpendicular",
    "Lagoon3": "Lagoon near\ncoast",
}

SITE_SHORT = {
    "Coast": "Coast",
    "Lagoon1": "L1",
    "Lagoon2": "L2",
    "Lagoon3": "L3",
}

SELECTED_DAS_CHANNELS = {
    "Coast": 413,
    "Lagoon1": 67,
    "Lagoon2": 130,
    "Lagoon3": 186,
}

DAS_NEIGHBOR_WINDOWS = {
    "Coast": (408, 418),
    "Lagoon1": (62, 72),
    "Lagoon2": (125, 135),
    "Lagoon3": (181, 191),
}

RBR_SERIAL_TO_SITE = {
    "204920": "Coast",
    "211699": "Lagoon1",
    "204928": "Lagoon2",
    "204923": "Lagoon3",
}

BANDS = OrderedDict([
    ("Diurnal", (0.65, 1.35)),
    ("Semidiurnal", (1.65, 2.35)),
    ("Terdiurnal", (2.65, 3.35)),
    ("Quarter-diurnal", (3.65, 4.35)),
])

CONSTITUENTS = OrderedDict([
    ("O1", 0.9295357067),
    ("K1", 1.0027379067),
    ("M2", 1.9322736133),
    ("S2", 2.0000000000),
    ("MK3", 2.9350115267),
    ("SK3", 3.0027379067),
    ("M4", 3.8645472267),
    ("MS4", 3.9322736133),
])

# 频谱计算阶段的临时时间步长。
# 默认 1 = 完全逐分钟 Lomb–Scargle。
# 若邻近通道计算明显过慢，可改为 5；对 0.7–4.25 cpd 潮汐频带仍足够。
ANALYSIS_RESAMPLE_MINUTES = 1
MIN_BIN_COVERAGE = 0.60

DAS_USE_LOG10 = True
DAS_WINSOR_PERCENTILES = (0.5, 99.5)

# 频率网格。约 5 天记录的 Rayleigh 分辨率约 0.2 cpd，
# 因而 1200 个频率点已经远密于实际可分辨能力。
FREQ_MIN_CPD = 0.20
FREQ_MAX_CPD = 4.30
FREQ_POINTS = 1200

# RBR 水位参数
WATER_DENSITY_KG_M3 = 1025.0
GRAVITY_M_S2 = 9.80665
RBR_WET_PRESSURE_THRESHOLD_DBAR = 0.03

ATM_ROLLING_HOURS = 12
ATM_LOW_QUANTILE = 0.08
ATM_PRESSURE_TOLERANCE_DBAR = 0.020
ATM_MIN_DRY_SEGMENT_MINUTES = 8

# 末尾回收异常
FINAL_TRIM_MINUTES = {
    "Coast": 15,
    "Lagoon1": 15,
    "Lagoon2": 15,
    "Lagoon3": 15,
}
FINAL_SEARCH_LOOKBACK_MINUTES = 120

# Lagoon1 中途取出检查
AUTO_DETECT_INTERNAL_HANDLING = True
INTERNAL_DETECTION_SITES = ["Lagoon1"]
MIN_ABRUPT_JUMP_M = 0.08
JUMP_MAD_FACTOR = 8.0
MIN_EVENT_DROP_M = 0.08
MAX_INTERNAL_EVENT_MINUTES = 180
PRE_EVENT_BASELINE_MINUTES = 20
POST_EVENT_BASELINE_MINUTES = 20
EVENT_PADDING_MINUTES = 2

MANUAL_EXCLUDE_INTERVALS = {
    "Coast": [],
    "Lagoon1": [
        # 需要时可手动追加：
        # ("2026-07-06 12:20", "2026-07-06 12:45", "mid-campaign inspection"),
    ],
    "Lagoon2": [],
    "Lagoon3": [],
}


# -----------------------------------------------------------------------------
# UI / 文件发现
# -----------------------------------------------------------------------------
def choose_root_folder(title: str = "选择 DAS + RBR 分析总文件夹") -> Path:
    """弹窗选择总文件夹。取消时退出。"""
    try:
        import tkinter as tk
        from tkinter import filedialog

        app = tk.Tk()
        app.withdraw()
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
        # 没有 Tk 时退回当前目录
        return Path.cwd().resolve()


def _candidate_das_dirs(root: Path) -> List[Path]:
    out = []
    for p in [root / "DAS", root / "DAS_1min", root]:
        if p.exists() and p not in out:
            out.append(p)
    return out


def discover_das_files(root: Path) -> Dict[str, Path]:
    """优先寻找 *_1min_RMS.xlsx，避免误读旧 30-min merged 文件。"""
    files: List[Path] = []
    for d in _candidate_das_dirs(root):
        files.extend([p for p in d.glob("*.xlsx") if not p.name.startswith("~$")])

    # 去重
    files = list(dict.fromkeys(p.resolve() for p in files))

    def score_file(site: str, p: Path) -> int:
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
            if token in t:
                score += 100
            elif site == "Lagoon1" and re.search(r"lagoon[_-]?1", t):
                score += 100
            elif site == "Lagoon2" and re.search(r"lagoon[_-]?2", t):
                score += 100
            elif site == "Lagoon3" and re.search(r"lagoon[_-]?3", t):
                score += 100
            else:
                return -1
        return score

    found: Dict[str, Path] = {}
    for site in SITES:
        ranked = [(score_file(site, p), p) for p in files]
        ranked = [(s, p) for s, p in ranked if s >= 0]
        ranked.sort(key=lambda x: (-x[0], str(x[1])))
        if not ranked:
            raise FileNotFoundError(
                f"没有找到 {site} 的 1-min RMS Excel。\n"
                f"建议文件名：{site}_1min_RMS.xlsx"
            )
        best_score = ranked[0][0]
        best = [p for s, p in ranked if s == best_score]
        if len(best) > 1:
            raise RuntimeError(
                f"{site} 匹配到多个同优先级文件，请只保留一个：\n" +
                "\n".join(str(p) for p in best)
            )
        found[site] = ranked[0][1]
    return found


def _candidate_rbr_dirs(root: Path) -> List[Path]:
    out = []
    for p in [root / "RBR", root]:
        if p.exists() and p not in out:
            out.append(p)
    return out


def discover_rbr_files(root: Path) -> Dict[str, Path]:
    files: List[Path] = []
    for d in _candidate_rbr_dirs(root):
        files.extend(d.glob("*.rsk"))
    files = list(dict.fromkeys(p.resolve() for p in files))
    if not files:
        raise FileNotFoundError("没有找到任何 .rsk 文件。建议放在 RBR 子文件夹中。")

    found: Dict[str, Path] = {}
    for path in files:
        con = sqlite3.connect(str(path))
        try:
            row = con.execute("SELECT serialID FROM instruments LIMIT 1").fetchone()
        finally:
            con.close()
        if row is None:
            continue
        serial = str(row[0])
        if serial in RBR_SERIAL_TO_SITE:
            found[RBR_SERIAL_TO_SITE[serial]] = path

    missing = [s for s in SITES if s not in found]
    if missing:
        raise FileNotFoundError("缺少以下 RBR 点位：" + ", ".join(missing))
    return found


# -----------------------------------------------------------------------------
# DAS 1-min RMS 读取与处理
# -----------------------------------------------------------------------------
def read_das_matrix(path: Path) -> pd.DataFrame:
    """读取 RMS_1min sheet；若 sheet 不存在则退回第一个 sheet。"""
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
    # 同一时间列若重复，取平均
    raw = raw.T.groupby(level=0).mean().T
    return raw


def robust_winsorize(series: pd.Series, low: float, high: float) -> pd.Series:
    s = series.astype(float).copy()
    valid = s.dropna()
    if len(valid) < 10:
        return s
    lo, hi = np.nanpercentile(valid.to_numpy(float), [low, high])
    return s.clip(lower=lo, upper=hi)


def resample_series(series: pd.Series,
                    minutes: int = ANALYSIS_RESAMPLE_MINUTES,
                    statistic: str = "median") -> pd.Series:
    """
    频谱阶段临时降到 minutes；不修改源 1-min Excel。
    对每个 bin 要求至少 MIN_BIN_COVERAGE 的 1-min 点存在。
    """
    s = series.sort_index().astype(float)
    if minutes <= 1:
        return s.groupby(level=0).mean().sort_index()

    rule = f"{int(minutes)}min"
    count = s.resample(rule).count()
    need = max(1, int(math.ceil(minutes * MIN_BIN_COVERAGE)))
    if statistic == "mean":
        out = s.resample(rule).mean()
    else:
        out = s.resample(rule).median()
    out[count < need] = np.nan
    return out


def process_das_series_for_shape(series: pd.Series,
                                 resample_minutes: int = ANALYSIS_RESAMPLE_MINUTES) -> pd.Series:
    """频带组成 / 谱形：正 RMS -> log10 -> winsor -> 临时重采样。"""
    s = series.astype(float).replace([np.inf, -np.inf], np.nan).copy()
    s[s <= 0] = np.nan
    if DAS_USE_LOG10:
        s = np.log10(s)
    s = robust_winsorize(s, *DAS_WINSOR_PERCENTILES)
    s = resample_series(s, resample_minutes, statistic="median")
    return s


def process_das_series_for_spatial(series: pd.Series,
                                   resample_minutes: int = ANALYSIS_RESAMPLE_MINUTES) -> pd.Series:
    """
    空间绝对格局：保留原始 RMS 尺度，不 log10、不逐序列标准化。
    仅把 1-min RMS 在频谱计算阶段临时做 mean aggregation。
    """
    s = series.astype(float).replace([np.inf, -np.inf], np.nan).copy()
    s[s <= 0] = np.nan
    return resample_series(s, resample_minutes, statistic="mean")


def prepare_das_representatives(das_files: Dict[str, Path],
                                resample_minutes: int = ANALYSIS_RESAMPLE_MINUTES
                                ) -> Dict[str, pd.Series]:
    out = {}
    for site in SITES:
        raw = read_das_matrix(das_files[site])
        ch = SELECTED_DAS_CHANNELS[site]
        if ch not in raw.index:
            raise KeyError(f"{das_files[site].name} 中没有 Ch {ch}")
        out[site] = process_das_series_for_shape(raw.loc[ch], resample_minutes)
    return crop_das_to_common_window(out)


def prepare_das_neighbors(das_files: Dict[str, Path],
                          mode: str = "shape",
                          resample_minutes: int = ANALYSIS_RESAMPLE_MINUTES
                          ) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    for site in SITES:
        raw = read_das_matrix(das_files[site])
        center = SELECTED_DAS_CHANNELS[site]
        c0, c1 = DAS_NEIGHBOR_WINDOWS[site]
        channels = [ch for ch in raw.index if c0 <= ch <= c1]
        if center not in channels:
            raise KeyError(f"{site} 中没有中心通道 Ch {center}")
        if len(channels) == 0:
            raise KeyError(f"{site} 中没有邻域 Ch {c0}-{c1}")

        series_map = {}
        for ch in channels:
            if mode == "spatial":
                s = process_das_series_for_spatial(raw.loc[ch], resample_minutes)
            else:
                s = process_das_series_for_shape(raw.loc[ch], resample_minutes)
            series_map[int(ch)] = s

        out[site] = {
            "center_channel": center,
            "channels": [int(ch) for ch in channels],
            "series": series_map,
        }

    # 统一四个 DAS 位置的名义时间窗
    starts, ends = [], []
    for site in SITES:
        for s in out[site]["series"].values():
            v = s.dropna()
            if len(v):
                starts.append(v.index.min())
                ends.append(v.index.max())
    if starts and ends:
        start, end = max(starts), min(ends)
        for site in SITES:
            for ch in list(out[site]["series"]):
                out[site]["series"][ch] = out[site]["series"][ch].loc[start:end]
        for site in SITES:
            out[site]["analysis_start"] = start
            out[site]["analysis_end"] = end
    return out


def crop_das_to_common_window(series_dict: Dict[str, pd.Series]) -> Dict[str, pd.Series]:
    starts, ends = [], []
    for s in series_dict.values():
        v = s.dropna()
        if len(v):
            starts.append(v.index.min())
            ends.append(v.index.max())
    if not starts or not ends:
        return series_dict
    start, end = max(starts), min(ends)
    return {site: s.loc[start:end] for site, s in series_dict.items()}


# -----------------------------------------------------------------------------
# RBR 1-min 读取与清理（沿用旧代码的规则，但分析不再强制 30-min）
# -----------------------------------------------------------------------------
def aggregate_rsk_one_minute(path: Path, site: str) -> pd.DataFrame:
    print(f"读取并聚合 RBR：{SITE_LABELS[site]} | {path.name}", flush=True)
    con = sqlite3.connect(str(path))
    try:
        query = """
            SELECT
                CAST(tstamp / 60000 AS INTEGER) * 60000 AS tbin,
                AVG(channel01) AS pressure_dbar,
                AVG(channel01 * channel01) AS pressure_sq_mean,
                COUNT(*) AS sample_count
            FROM data
            GROUP BY tbin
            ORDER BY tbin
        """
        df = pd.read_sql_query(query, con)
    finally:
        con.close()

    if df.empty:
        raise RuntimeError(f"{path.name} 中没有 RBR data。")
    df["time"] = pd.to_datetime(df["tbin"], unit="ms")
    var = np.maximum(
        df["pressure_sq_mean"].to_numpy(float)
        - df["pressure_dbar"].to_numpy(float) ** 2,
        0.0,
    )
    df["pressure_std_dbar"] = np.sqrt(var)
    df["site"] = site
    return (
        df.drop(columns=["tbin", "pressure_sq_mean"])
        .set_index("time")
        .sort_index()
    )


def consecutive_true_segments(mask: pd.Series):
    m = mask.fillna(False).astype(bool)
    if m.empty:
        return []
    group = (m != m.shift(fill_value=False)).cumsum()
    out = []
    for _, part in m.groupby(group):
        if bool(part.iloc[0]):
            out.append((part.index[0], part.index[-1], len(part)))
    return out


def estimate_atmospheric_pressure(coast_1min: pd.DataFrame) -> pd.Series:
    p = coast_1min["pressure_dbar"].astype(float).sort_index()
    std = coast_1min["pressure_std_dbar"].astype(float).reindex(p.index)
    coast_end = p.index.max()
    atmosphere_input_mask = p.index < (
        coast_end - pd.Timedelta(minutes=FINAL_TRIM_MINUTES["Coast"])
    )

    window = max(int(round(ATM_ROLLING_HOURS * 60.0)), 60)
    min_periods = max(int(window * 0.25), 60)
    lower_env = p.rolling(
        window=window, center=True, min_periods=min_periods
    ).quantile(ATM_LOW_QUANTILE)
    lower_env = lower_env.interpolate(method="time", limit_direction="both")

    if std.notna().sum():
        std_threshold = float(np.nanpercentile(std.dropna().to_numpy(float), 45))
    else:
        std_threshold = 0.01
    std_threshold = min(max(std_threshold, 0.002), 0.020)

    candidate = (
        (p <= lower_env + ATM_PRESSURE_TOLERANCE_DBAR)
        & (std <= std_threshold)
        & atmosphere_input_mask
    )

    accepted = pd.Series(False, index=p.index)
    for start, end, n in consecutive_true_segments(candidate):
        if n >= ATM_MIN_DRY_SEGMENT_MINUTES:
            accepted.loc[start:end] = True

    if accepted.sum() >= 30:
        anchors = p.where(accepted).rolling(window=15, center=True, min_periods=1).median()
        atmosphere = anchors.interpolate(method="time", limit_direction="both")
        atmosphere = atmosphere.rolling(window=121, center=True, min_periods=1).median()
    else:
        atmosphere = lower_env.rolling(window=121, center=True, min_periods=1).median()
    return atmosphere


def correct_rbr_water_level(one_min: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    atmosphere = estimate_atmospheric_pressure(one_min["Coast"])
    depth_factor = 10000.0 / (WATER_DENSITY_KG_M3 * GRAVITY_M_S2)
    corrected = {}
    for site in SITES:
        df = one_min[site].copy()
        union = df.index.union(atmosphere.index)
        atm = (
            atmosphere.reindex(union)
            .interpolate(method="time", limit_direction="both")
            .reindex(df.index)
        )
        df["estimated_atmosphere_dbar"] = atm
        df["water_pressure_dbar"] = df["pressure_dbar"] - atm
        df["corrected_water_level_m"] = df["water_pressure_dbar"] * depth_factor
        df["display_water_level_m"] = df["corrected_water_level_m"].clip(lower=0.0)
        if site == "Coast":
            df["wet"] = df["water_pressure_dbar"] > RBR_WET_PRESSURE_THRESHOLD_DBAR
        else:
            df["wet"] = df["water_pressure_dbar"] > -0.02
        corrected[site] = df
    return corrected


def robust_mad(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan
    med = float(np.median(arr))
    return float(np.median(np.abs(arr - med)))


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = []
    s0, e0, r0 = intervals[0]
    for s, e, r in intervals[1:]:
        if s <= e0 + pd.Timedelta(minutes=1):
            e0 = max(e0, e)
            if r not in r0:
                r0 += "; " + r
        else:
            merged.append((s0, e0, r0))
            s0, e0, r0 = s, e, r
    merged.append((s0, e0, r0))
    return merged


def detect_internal_handling_intervals(level: pd.Series, site: str):
    if not AUTO_DETECT_INTERNAL_HANDLING or site not in INTERNAL_DETECTION_SITES:
        return []
    s = level.astype(float).sort_index()
    s = s.interpolate(method="time", limit=2, limit_direction="both")
    smooth = s.rolling(window=3, center=True, min_periods=1).median()
    diff = smooth.diff()
    mad_diff = robust_mad(diff.to_numpy(float))
    robust_sigma = 1.4826 * mad_diff if np.isfinite(mad_diff) else 0.0
    threshold = max(MIN_ABRUPT_JUMP_M, JUMP_MAD_FACTOR * robust_sigma)
    downs = list(diff.index[diff <= -threshold])
    ups = list(diff.index[diff >= 0.65 * threshold])
    intervals = []
    used_until = None
    for down in downs:
        if used_until is not None and down <= used_until:
            continue
        minutes_to_end = (s.index.max() - down).total_seconds() / 60.0
        if minutes_to_end <= FINAL_SEARCH_LOOKBACK_MINUTES:
            continue
        search_end = down + pd.Timedelta(minutes=MAX_INTERNAL_EVENT_MINUTES)
        candidate_ups = [t for t in ups if down < t <= search_end]
        if not candidate_ups:
            continue
        up = candidate_ups[0]
        before = smooth.loc[
            down - pd.Timedelta(minutes=PRE_EVENT_BASELINE_MINUTES):
            down - pd.Timedelta(minutes=1)
        ].dropna()
        during = smooth.loc[down:up].dropna()
        after = smooth.loc[
            up + pd.Timedelta(minutes=1):
            up + pd.Timedelta(minutes=POST_EVENT_BASELINE_MINUTES)
        ].dropna()
        if len(before) < 5 or len(during) < 1 or len(after) < 5:
            continue
        drop_before = float(before.median() - during.median())
        drop_after = float(after.median() - during.median())
        duration = (up - down).total_seconds() / 60.0
        if (
            drop_before >= MIN_EVENT_DROP_M
            and drop_after >= MIN_EVENT_DROP_M
            and 1 <= duration <= MAX_INTERNAL_EVENT_MINUTES
        ):
            start = down - pd.Timedelta(minutes=EVENT_PADDING_MINUTES)
            end = up + pd.Timedelta(minutes=EVENT_PADDING_MINUTES)
            intervals.append((start, end, "automatic mid-campaign removal/inspection"))
            used_until = end
    return merge_intervals(intervals)


def detect_final_retrieval_interval(level: pd.Series, site: str):
    s = level.astype(float).sort_index()
    end = s.index.max()
    fixed_start = end - pd.Timedelta(minutes=int(FINAL_TRIM_MINUTES[site]))

    # 所有点位都尝试在末尾窗口寻找“分钟级突降”。
    # Coast 虽然会正常低潮出露，但自然潮位通常不会在 1 min 内骤降 >= 8 cm；
    # 因此这里使用突变阈值而不是把所有低水位都判为回收。
    tail = s.loc[end - pd.Timedelta(minutes=FINAL_SEARCH_LOOKBACK_MINUTES):end]
    tail = tail.interpolate(method="time", limit=2, limit_direction="both")
    smooth = tail.rolling(window=3, center=True, min_periods=1).median()
    diff = smooth.diff()
    mad_diff = robust_mad(diff.to_numpy(float))
    robust_sigma = 1.4826 * mad_diff if np.isfinite(mad_diff) else 0.0
    threshold = max(MIN_ABRUPT_JUMP_M, JUMP_MAD_FACTOR * robust_sigma)
    negative = diff[diff <= -threshold]
    if not negative.empty:
        jump_time = negative.idxmin()
        automatic_start = jump_time - pd.Timedelta(minutes=EVENT_PADDING_MINUTES)
        start = min(fixed_start, automatic_start)
        reason = "automatic final retrieval + fixed tail trim"
    else:
        start = fixed_start
        reason = "fixed final tail trim"
    return start, end, reason


def parse_manual_intervals(site: str):
    out = []
    for item in MANUAL_EXCLUDE_INTERVALS.get(site, []):
        if len(item) != 3:
            raise ValueError("手动区间必须为 (start, end, reason)")
        s, e, r = pd.Timestamp(item[0]), pd.Timestamp(item[1]), str(item[2])
        if e < s:
            raise ValueError(f"无效手动区间：{item}")
        out.append((s, e, r))
    return out


def clean_handling_artifacts(corrected: Dict[str, pd.DataFrame]):
    cleaned = {}
    rows = []
    for site in SITES:
        df = corrected[site].copy()
        level = df["display_water_level_m"]
        intervals = []
        intervals.extend(detect_internal_handling_intervals(level, site))
        intervals.append(detect_final_retrieval_interval(level, site))
        intervals.extend(parse_manual_intervals(site))
        intervals = merge_intervals(intervals)

        df["excluded_handling"] = False
        for start, end, reason in intervals:
            actual_start = max(start, df.index.min())
            actual_end = min(end, df.index.max())
            if actual_end < actual_start:
                continue
            mask = (df.index >= actual_start) & (df.index <= actual_end)
            df.loc[mask, "excluded_handling"] = True
            rows.append({
                "site": site,
                "location": SITE_FIBER_LABEL[site],
                "start": actual_start,
                "end": actual_end,
                "duration_min": (actual_end - actual_start).total_seconds() / 60.0 + 1.0,
                "reason": reason,
            })

        df["clean_water_level_m"] = df["display_water_level_m"].where(~df["excluded_handling"])
        # Coast 正常低潮出露不是“异常”，但出露时不能作为水位用于潮汐频谱。
        if site == "Coast":
            df.loc[~df["wet"], "clean_water_level_m"] = np.nan
        cleaned[site] = df

    intervals_df = pd.DataFrame(rows)
    if not intervals_df.empty:
        intervals_df = intervals_df.sort_values(["site", "start"]).reset_index(drop=True)
    return cleaned, intervals_df


def load_and_clean_rbr(root: Path):
    files = discover_rbr_files(root)
    one_min = {site: aggregate_rsk_one_minute(files[site], site) for site in SITES}
    corrected = correct_rbr_water_level(one_min)
    cleaned, intervals = clean_handling_artifacts(corrected)
    return cleaned, intervals, files


def common_rbr_nominal_window(cleaned: Dict[str, pd.DataFrame]):
    starts, ends = [], []
    for site in SITES:
        s = cleaned[site]["clean_water_level_m"].dropna()
        if len(s):
            starts.append(s.index.min())
            ends.append(s.index.max())
    if not starts or not ends:
        raise RuntimeError("RBR 清理后没有足够有效水位。")
    return max(starts), min(ends)


def linear_detrend_series(series: pd.Series) -> pd.Series:
    s = series.astype(float).replace([np.inf, -np.inf], np.nan)
    v = s.dropna()
    out = pd.Series(np.nan, index=s.index, dtype=float)
    if len(v) < 3:
        return out
    t = (v.index - v.index[0]).total_seconds().to_numpy() / 86400.0
    coef = np.polyfit(t, v.to_numpy(float), 1)
    anomaly = v.to_numpy(float) - np.polyval(coef, t)
    anomaly -= np.nanmean(anomaly)
    out.loc[v.index] = anomaly
    return out


def prepare_rbr_analysis(cleaned: Dict[str, pd.DataFrame],
                         resample_minutes: int = ANALYSIS_RESAMPLE_MINUTES
                         ) -> Tuple[Dict[str, pd.DataFrame], pd.Timestamp, pd.Timestamp]:
    start, end = common_rbr_nominal_window(cleaned)
    output = {}
    for site in SITES:
        s1 = cleaned[site]["clean_water_level_m"].loc[start:end].copy()
        if resample_minutes <= 1:
            s = s1
        else:
            rule = f"{int(resample_minutes)}min"
            count = s1.resample(rule).count()
            s = s1.resample(rule).mean()
            need = max(1, int(math.ceil(resample_minutes * MIN_BIN_COVERAGE)))
            s[count < need] = np.nan
        anomaly = linear_detrend_series(s)
        output[site] = pd.DataFrame({
            "water_level_m": s,
            "level_anomaly_m": anomaly,
        })
    return output, start, end


# -----------------------------------------------------------------------------
# 频谱工具
# -----------------------------------------------------------------------------
def _valid_t_y(index: pd.DatetimeIndex, values: Sequence[float]):
    s = pd.Series(np.asarray(values, dtype=float), index=pd.DatetimeIndex(index))
    s = s.replace([np.inf, -np.inf], np.nan).dropna().sort_index()
    if len(s) < 3:
        return np.array([]), np.array([]), s
    t = (s.index - s.index[0]).total_seconds().to_numpy() / 86400.0
    return t.astype(float), s.to_numpy(float), s


def detrend_array(t_days: np.ndarray, y: np.ndarray, standardize: bool) -> np.ndarray:
    if len(y) < 3:
        return np.full_like(y, np.nan, dtype=float)
    coef = np.polyfit(t_days, y, 1)
    out = y - np.polyval(coef, t_days)
    out -= np.nanmean(out)
    if standardize:
        sd = np.nanstd(out)
        if sd > 0:
            out = out / sd
    return out


def lomb_shape_spectrum(index: pd.DatetimeIndex,
                        values: Sequence[float],
                        freq_cpd: np.ndarray) -> Tuple[np.ndarray, int, float]:
    """用于谱形/频带组成：detrend + z-scale + normalized Lomb–Scargle。"""
    t, y, s = _valid_t_y(index, values)
    if len(y) < 20:
        return np.full_like(freq_cpd, np.nan), len(y), np.nan
    yd = detrend_array(t, y, standardize=True)
    if np.nanstd(yd) == 0:
        return np.full_like(freq_cpd, np.nan), len(y), np.nan
    omega = 2.0 * np.pi * np.asarray(freq_cpd, dtype=float)
    try:
        p = lombscargle(t, yd, omega, precenter=False, normalize=True)
    except TypeError:
        p = lombscargle(t, yd, omega, precenter=False)
        var = np.var(yd)
        if var > 0:
            p = 2.0 * p / (len(yd) * var)
    p = np.maximum(np.asarray(p, dtype=float), 0.0)
    duration = (s.index.max() - s.index.min()).total_seconds() / 86400.0
    return p, len(y), float(duration)


def lomb_absolute_spectrum(index: pd.DatetimeIndex,
                           values: Sequence[float],
                           freq_cpd: np.ndarray) -> Tuple[np.ndarray, int, float]:
    """
    用于 Figure 4B 空间绝对格局：detrend，但不按每条序列标准化。
    raw Lomb power 除以 N，使不同缺测数下的幅值尺度更可比。
    """
    t, y, s = _valid_t_y(index, values)
    if len(y) < 20:
        return np.full_like(freq_cpd, np.nan), len(y), np.nan
    yd = detrend_array(t, y, standardize=False)
    omega = 2.0 * np.pi * np.asarray(freq_cpd, dtype=float)
    p = lombscargle(t, yd, omega, precenter=False, normalize=False)
    p = 2.0 * np.asarray(p, dtype=float) / max(len(yd), 1)
    p = np.maximum(p, 0.0)
    duration = (s.index.max() - s.index.min()).total_seconds() / 86400.0
    return p, len(y), float(duration)


def gaussian_smooth_frequency(power: np.ndarray,
                              freq: np.ndarray,
                              fwhm_cpd: float) -> np.ndarray:
    p = np.asarray(power, dtype=float)
    f = np.asarray(freq, dtype=float)
    if len(p) < 5 or not np.isfinite(fwhm_cpd) or fwhm_cpd <= 0:
        return p.copy()
    df = float(np.nanmedian(np.diff(f)))
    sigma_f = fwhm_cpd / 2.354820045
    sigma_pts = max(sigma_f / df, 0.5)
    radius = int(max(3, math.ceil(4.0 * sigma_pts)))
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (x / sigma_pts) ** 2)
    kernel /= kernel.sum()
    padded = np.pad(p, radius, mode="edge")
    return np.convolve(padded, kernel, mode="same")[radius:-radius]


def trapz_compat(y, x) -> float:
    yy = np.asarray(y, dtype=float)
    xx = np.asarray(x, dtype=float)
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(yy, xx))
    return float(np.trapz(yy, xx))


def band_integrals(freq: np.ndarray, power: np.ndarray) -> Dict[str, float]:
    out = {}
    for band, (f0, f1) in BANDS.items():
        m = (freq >= f0) & (freq <= f1)
        if m.sum() < 2 or np.all(~np.isfinite(power[m])):
            out[band] = np.nan
        else:
            out[band] = trapz_compat(power[m], freq[m])
    return out


def composition(values: Sequence[float]) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    valid = np.isfinite(x) & (x >= 0)
    out = np.full_like(x, np.nan)
    total = np.nansum(x[valid])
    if total > 0:
        out[valid] = x[valid] / total
    return out


def cosine_similarity(x: Sequence[float], y: Sequence[float]) -> float:
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    m = np.isfinite(xx) & np.isfinite(yy)
    if m.sum() < 2:
        return np.nan
    xx, yy = xx[m], yy[m]
    den = np.linalg.norm(xx) * np.linalg.norm(yy)
    if den <= 0:
        return np.nan
    return float(np.dot(xx, yy) / den)


def mae_percentage_points(x: Sequence[float], y: Sequence[float]) -> float:
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    m = np.isfinite(xx) & np.isfinite(yy)
    if not np.any(m):
        return np.nan
    return float(100.0 * np.mean(np.abs(xx[m] - yy[m])))


def determine_common_rayleigh(das_series: Dict[str, pd.Series],
                              rbr_analysis: Dict[str, pd.DataFrame]) -> Tuple[float, float]:
    durations = []
    for s in das_series.values():
        v = s.dropna()
        if len(v) > 1:
            durations.append((v.index.max() - v.index.min()).total_seconds() / 86400.0)
    for site in SITES:
        v = rbr_analysis[site]["level_anomaly_m"].dropna()
        if len(v) > 1:
            durations.append((v.index.max() - v.index.min()).total_seconds() / 86400.0)
    if not durations:
        raise RuntimeError("没有足够数据确定 Rayleigh 分辨率。")
    common_duration = float(min(durations))
    return 1.0 / common_duration, common_duration


def freq_grid(points: int = FREQ_POINTS,
              fmin: float = FREQ_MIN_CPD,
              fmax: float = FREQ_MAX_CPD) -> np.ndarray:
    return np.linspace(fmin, fmax, int(points))


def processing_info_dataframe(das_files: Dict[str, Path],
                              rbr_files: Dict[str, Path],
                              das_series: Dict[str, pd.Series] | None,
                              rbr_start,
                              rbr_end,
                              rayleigh: float | None = None) -> pd.DataFrame:
    rows = [
        ("DAS source", "1-min RMS derived from raw 200-Hz DAS"),
        ("spectral analysis temporary cadence (min)", ANALYSIS_RESAMPLE_MINUTES),
        ("DAS log10 for spectral shape", DAS_USE_LOG10),
        ("DAS winsor percentiles", str(DAS_WINSOR_PERCENTILES)),
        ("RBR analysis", "cleaned water level, no mandatory 30-min aggregation"),
        ("RBR nominal common start", rbr_start),
        ("RBR nominal common end", rbr_end),
        ("Rayleigh smoothing FWHM (cpd)", rayleigh),
    ]
    if das_series:
        starts = [s.dropna().index.min() for s in das_series.values() if len(s.dropna())]
        ends = [s.dropna().index.max() for s in das_series.values() if len(s.dropna())]
        if starts and ends:
            rows.extend([
                ("DAS common start", max(starts)),
                ("DAS common end", min(ends)),
            ])
    for site in SITES:
        rows.append((f"DAS file {site}", str(das_files[site])))
        rows.append((f"RBR file {site}", str(rbr_files[site])))
    return pd.DataFrame(rows, columns=["parameter", "value"])
