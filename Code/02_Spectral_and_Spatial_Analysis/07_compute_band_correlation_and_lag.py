#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
07 | 宽频带相关、峰值相关与高质量滞后（v3）
==========================================

目的：
不再把“单一频率正弦相位”作为传播证据，而是在 D/SD/TD/QD 宽频带上比较：

A. DAS vs 2025 同期实际潮位
- zero-lag correlation
- peak correlation
- lag at peak correlation

B. Lagoon DAS vs Coast DAS center
- zero-lag correlation
- peak correlation
- lag at peak correlation

高质量滞后门槛：
- tide lag:
    peak corr >= 0.40
    且该频带在 2025 潮位中的 variance fraction >= 0.05
- Coast-DAS lag:
    peak corr >= 0.45

这样可避免把弱 TD/QD reference 的相位噪声误解释成传播。
"""

from pathlib import Path
import numpy as np
import pandas as pd

# 公共计算模块随代码目录一起保留。
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))

from das_rbr_1min_common import (
    SITES, SITE_SHORT, SITE_LABELS, BANDS,
    choose_root_folder, discover_das_files, read_das_matrix
)
from geometry_config_v3 import REPRESENTATIVE_CHANNELS
from mechanism_common_v3 import choose_excel_file, read_tide_excel, make_uniform_time_grid, resample_series_to_grid, preprocess_signal_series, fft_bandpass_uniform, band_variance_fraction, pearsonr_safe, max_corr_with_lag, UNIFORM_DT_MIN, BAND_NAMES

OUTDIR = "07_band_limited_transfer_output"


def series_from_raw(raw, channel):
    if int(channel) not in raw.index:
        raise KeyError(f"通道 {channel} 不在 DAS Excel 中。")
    s = raw.loc[int(channel)].astype(float)
    s.index = pd.DatetimeIndex(raw.columns)
    return s.sort_index()


def main():

    root = choose_root_folder("选择包含四个 DAS 1-min RMS Excel 的总文件夹")
    tide_path = choose_excel_file("选择 2025 DAS 同期潮位 Excel", root)

    das_files = discover_das_files(root)
    raw_by_site = {site: read_das_matrix(das_files[site]) for site in SITES}
    tide = read_tide_excel(tide_path)

    starts = [raw.columns.min() for raw in raw_by_site.values()]
    ends = [raw.columns.max() for raw in raw_by_site.values()]
    common_start = max(max(starts), tide.index.min())
    common_end = min(min(ends), tide.index.max())

    grid = make_uniform_time_grid(common_start, common_end, UNIFORM_DT_MIN)
    if len(grid) < 100:
        raise RuntimeError(
            f"共同时间窗过短：{common_start} -> {common_end}。"
            "请检查 DAS 与潮位时间覆盖。"
        )

    # Tide preprocessing
    tide_u = resample_series_to_grid(tide, grid)
    tide_proc = preprocess_signal_series(
        tide_u,
        apply_log=False,
        winsor=False,
        do_detrend=True,
        do_zscore=True
    )

    tide_band = {}
    tide_varfrac = {}
    for bn, band in BANDS.items():
        arr = fft_bandpass_uniform(tide_proc.values, UNIFORM_DT_MIN, band)
        tide_band[bn] = arr
        tide_varfrac[bn] = band_variance_fraction(
            tide_proc.values, UNIFORM_DT_MIN, band
        )

    # Coast center reference
    coast_ch = REPRESENTATIVE_CHANNELS["Coast"]
    coast_s = series_from_raw(raw_by_site["Coast"], coast_ch)
    coast_u = resample_series_to_grid(coast_s, grid)
    coast_p = preprocess_signal_series(
        coast_u, apply_log=True, winsor=True, do_detrend=True, do_zscore=True
    )
    coast_band = {
        bn: fft_bandpass_uniform(coast_p.values, UNIFORM_DT_MIN, band)
        for bn, band in BANDS.items()
    }

    rows = []

    for site in SITES:
        raw = raw_by_site[site]
        for ch in sorted(raw.index.astype(int)):
            s = series_from_raw(raw, ch)
            su = resample_series_to_grid(s, grid)
            sp = preprocess_signal_series(
                su, apply_log=True, winsor=True, do_detrend=True, do_zscore=True
            )

            base = sp.values.astype(float)

            for bn, band in BANDS.items():
                sig = fft_bandpass_uniform(base, UNIFORM_DT_MIN, band)
                ref_tide = tide_band[bn]
                ref_coast = coast_band[bn]

                fcenter = 0.5 * (band[0] + band[1])
                period_h = 24.0 / fcenter
                max_lag_h = 0.5 * period_h

                r0_t = pearsonr_safe(ref_tide, sig)
                rp_t, lag_t = max_corr_with_lag(
                    ref_tide, sig, UNIFORM_DT_MIN, max_lag_h
                )

                r0_c = pearsonr_safe(ref_coast, sig)
                rp_c, lag_c = max_corr_with_lag(
                    ref_coast, sig, UNIFORM_DT_MIN, max_lag_h
                )

                rows.append({
                    "site": site,
                    "site_short": SITE_SHORT[site],
                    "site_label": SITE_LABELS[site],
                    "channel": int(ch),
                    "band": bn,
                    "band_fmin_cpd": band[0],
                    "band_fmax_cpd": band[1],
                    "band_center_cpd": fcenter,
                    "max_lag_hours_searched": max_lag_h,
                    "tide_band_variance_fraction": tide_varfrac[bn],
                    "corr0_vs_tide": r0_t,
                    "corrpeak_vs_tide": rp_t,
                    "lag_at_corrpeak_vs_tide_h": lag_t,
                    "corr0_vs_CoastDAS": r0_c,
                    "corrpeak_vs_CoastDAS": rp_c,
                    "lag_at_corrpeak_vs_CoastDAS_h": lag_c,
                    "signal_band_std": float(np.nanstd(sig, ddof=0)),
                })

    df = pd.DataFrame(rows)

    df["tide_reference_quality"] = np.where(
        df["tide_band_variance_fraction"] >= 0.05,
        "usable",
        "weak"
    )

    df["high_quality_tide_lag"] = (
        (df["corrpeak_vs_tide"] >= 0.40)
        & (df["tide_band_variance_fraction"] >= 0.05)
    )

    df["high_quality_coast_lag"] = (
        df["corrpeak_vs_CoastDAS"] >= 0.45
    )

    summary = (
        df.groupby(["site", "band"])
          .agg(
              n=("channel", "size"),
              tide_corr0_median=("corr0_vs_tide", "median"),
              tide_corrpeak_median=("corrpeak_vs_tide", "median"),
              coast_corr0_median=("corr0_vs_CoastDAS", "median"),
              coast_corrpeak_median=("corrpeak_vs_CoastDAS", "median"),
              tide_lag_median_h=("lag_at_corrpeak_vs_tide_h", "median"),
              coast_lag_median_h=("lag_at_corrpeak_vs_CoastDAS_h", "median"),
              tide_reference_varfrac=("tide_band_variance_fraction", "median"),
              high_quality_tide_lag_fraction=("high_quality_tide_lag", "mean"),
              high_quality_coast_lag_fraction=("high_quality_coast_lag", "mean"),
          )
          .reset_index()
    )

    outdir = root / OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(outdir / "band_limited_transfer_metrics.xlsx", engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="channel_metrics", index=False)
        summary.to_excel(writer, sheet_name="site_band_summary", index=False)
        pd.DataFrame({
            "parameter": [
                "tide_file",
                "common_start",
                "common_end",
                "uniform_dt_min",
                "tide_lag_quality_rule",
                "coast_lag_quality_rule",
            ],
            "value": [
                tide_path.name,
                str(grid.min()),
                str(grid.max()),
                UNIFORM_DT_MIN,
                "peak_corr>=0.40 AND tide_band_variance_fraction>=0.05",
                "peak_corr>=0.45",
            ]
        }).to_excel(writer, sheet_name="metadata", index=False)


    print("=" * 80)
    print("07 v3 完成")
    print("输出：", outdir)
    print("2025 tide band variance fractions:")
    for bn in BAND_NAMES:
        print(f"  {bn:<16s} {tide_varfrac[bn]:.4f}")
    print("=" * 80)


if __name__ == "__main__":
    main()
