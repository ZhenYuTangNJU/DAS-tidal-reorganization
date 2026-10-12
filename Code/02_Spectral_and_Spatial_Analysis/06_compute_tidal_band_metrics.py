#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06 | 全通道四频带空间响应（v3）
==============================

功能：
- 使用 200-Hz raw DAS -> 1-min RMS 的四个 Excel；
- 对全部分析通道计算 D/SD/TD/QD 四频带；
- 同时输出：
    1) within-channel band composition
    2) apparent absolute DAS band power

这是后续 STI / 几何控制 / 方向扫描的基础。
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
from mechanism_common_v3 import BAND_NAMES, preprocess_site_matrix, fft_band_metrics

OUTDIR = "06_all_channel_band_maps_output"


def main():
    root = choose_root_folder("选择包含四个 1-min RMS Excel 的总文件夹")
    outdir = root / OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)

    das_files = discover_das_files(root)

    rows = []
    qc_rows = []

    for site in SITES:
        print(f"处理 {site} ...")
        raw = read_das_matrix(das_files[site])

        # shape/composition
        times_s, ch_s, Xs, qc_s = preprocess_site_matrix(raw, mode="shape")
        f_s, p_s, bp_s, comp_s = fft_band_metrics(Xs)

        # apparent absolute power
        times_a, ch_a, Xa, qc_a = preprocess_site_matrix(raw, mode="absolute")
        f_a, p_a, bp_a, _ = fft_band_metrics(Xa)

        common_ch = sorted(set(ch_s).intersection(set(ch_a)))
        map_s = {int(ch): i for i, ch in enumerate(ch_s)}
        map_a = {int(ch): i for i, ch in enumerate(ch_a)}

        for ch in common_ch:
            i = map_s[ch]
            j = map_a[ch]
            row = {
                "site": site,
                "site_label": SITE_LABELS[site],
                "channel": int(ch),
            }
            for bi, band in enumerate(BAND_NAMES):
                row[f"{band}_composition"] = float(comp_s[i, bi])
                row[f"{band}_shape_band_power"] = float(bp_s[i, bi])
                row[f"{band}_absolute_band_power"] = float(bp_a[j, bi])
                row[f"{band}_log10_absolute_band_power"] = (
                    float(np.log10(bp_a[j, bi])) if bp_a[j, bi] > 0 else np.nan
                )
            rows.append(row)

        qc_s = qc_s.copy()
        qc_s.insert(0, "site", site)
        qc_s["mode"] = "shape"
        qc_rows.append(qc_s)

        qc_a = qc_a.copy()
        qc_a.insert(0, "site", site)
        qc_a["mode"] = "absolute"
        qc_rows.append(qc_a)

    metrics = pd.DataFrame(rows)
    qc = pd.concat(qc_rows, ignore_index=True)

    # nominal chainage in each site
    metrics["site_chainage_m"] = np.nan
    for site in SITES:
        m = metrics["site"] == site
        if m.any():
            cmin = metrics.loc[m, "channel"].min()
            metrics.loc[m, "site_chainage_m"] = (
                metrics.loc[m, "channel"] - cmin
            ) * 2.0


    xlsx = outdir / "all_channel_band_metrics.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as writer:
        metrics.to_excel(writer, sheet_name="all_channel_metrics", index=False)
        qc.to_excel(writer, sheet_name="channel_QC", index=False)
        pd.DataFrame({
            "band": BAND_NAMES,
            "f_low_cpd": [BANDS[b][0] for b in BAND_NAMES],
            "f_high_cpd": [BANDS[b][1] for b in BAND_NAMES],
        }).to_excel(writer, sheet_name="band_definition", index=False)

    print("\n06 完成：", outdir)
    print("关键输出：", xlsx)


if __name__ == "__main__":
    main()
