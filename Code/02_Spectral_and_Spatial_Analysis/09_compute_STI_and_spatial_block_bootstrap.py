#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
09 | STI、距离、方向与“潟湖效应”模型（v3）
==========================================

核心目标：
1. 定量描述 Coast -> Lagoon 后频谱组成被重构的程度；
2. 比较：
   - 欧氏距离到 Coast
   - 沿光纤 network distance
   - 光纤方向
   - 潟湖位置本身
   对空间响应的解释能力。

Spectral Transformation Index:
    STI = 0.5 * sum(|C_i(channel) - C_i(Coast median)|)

STI 越大，说明频带组成相对 Coast reference 的变化越强。

注意：
- “潟湖效应”模型只用于普通 R²/AICc 解释；
- 因为只有一个 Coast 组和三个 Lagoon 组，LOSO 对 binary lagoon indicator
  的外推意义较弱，因此该模型的 LOSO 结果不作为主证据。
"""

from pathlib import Path
import numpy as np
import pandas as pd

# 公共计算模块随代码目录一起保留。
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_common"))

from das_rbr_1min_common import (
    SITES, choose_root_folder
)
from mechanism_common_v3 import build_exact_length_layout, BAND_NAMES, spectral_transformation_index, cosine_similarity, logit, ols_fit, leave_one_site_out_rmse

STI_OUTDIR = "09_STI_geometry_models_v3_output"


def find_metrics(root):
    direct = root / "06_all_channel_band_maps_output" / "all_channel_band_metrics.xlsx"
    if direct.exists():
        return direct
    hits = list(root.rglob("all_channel_band_metrics.xlsx"))
    if not hits:
        raise FileNotFoundError("未找到 06 输出 all_channel_band_metrics.xlsx。")
    return hits[0]


def fit_model(work, ycol, xcols, model_name, do_cv=True):
    cols = [ycol] + list(xcols)
    sub = work.dropna(subset=cols).copy()
    if len(sub) < max(8, len(xcols) + 3):
        return None

    X = sub[xcols].to_numpy(float) if xcols else np.empty((len(sub), 0))
    fit = ols_fit(sub[ycol].to_numpy(float), X, xcols)
    if fit is None:
        return None

    coef = dict(zip(fit["coef_names"], fit["beta"]))

    if do_cv:
        cv_rmse, cv_n = leave_one_site_out_rmse(sub, ycol, xcols)
    else:
        cv_rmse, cv_n = np.nan, 0

    row = {
        "model": model_name,
        "n": fit["n"],
        "R2": fit["r2"],
        "AICc": fit["aicc"],
        "LOSO_RMSE": cv_rmse,
        "LOSO_n": cv_n,
        "predictors": ", ".join(xcols) if xcols else "(intercept only)",
    }
    for k, v in coef.items():
        row[f"coef_{k}"] = v
    return row


def compute_sti(root):
    root = Path(root)

    p06 = find_metrics(root)
    band = pd.read_excel(p06, sheet_name="all_channel_metrics")

    geom, geom_meta = build_exact_length_layout()
    geom_use = geom[geom["is_analysis_site"]].copy()

    df = band.merge(
        geom_use,
        on=["site", "channel"],
        how="left",
        suffixes=("", "_geom")
    )

    # Coast median composition reference
    coast = df[df["site"] == "Coast"].copy()
    ref = np.array(
        [coast[f"{bn}_composition"].median() for bn in BAND_NAMES],
        dtype=float
    )

    sti = []
    cosref = []
    for _, r in df.iterrows():
        v = np.array([r[f"{bn}_composition"] for bn in BAND_NAMES], dtype=float)
        sti.append(spectral_transformation_index(v, ref))
        cosref.append(cosine_similarity(v, ref))

    df["STI_vs_CoastMedian"] = sti
    df["composition_cosine_vs_CoastMedian"] = cosref

    df["log_orientation_to_coast_normal"] = np.log(
        np.clip(df["orientation_cos2_to_coast_normal"].to_numpy(float), 1e-4, 1.0)
    )
    df["log_orientation_to_coast_axis"] = np.log(
        np.clip(df["orientation_cos2_to_coast_axis"].to_numpy(float), 1e-4, 1.0)
    )

    # Model definitions
    model_defs = [
        ("M0_intercept", [], True),
        ("M1_euclidean_distance", ["distance_to_coast_polyline_m"], True),
        ("M2_network_distance", ["network_distance_to_CH330_m"], True),
        ("M3_orientation_normal", ["log_orientation_to_coast_normal"], True),
        ("M4_euclidean_plus_orientation", [
            "distance_to_coast_polyline_m",
            "log_orientation_to_coast_normal"
        ], True),
        ("M5_network_plus_orientation", [
            "network_distance_to_CH330_m",
            "log_orientation_to_coast_normal"
        ], True),
        ("M6_euclidean_orientation_lagoon", [
            "distance_to_coast_polyline_m",
            "log_orientation_to_coast_normal",
            "inside_lagoon"
        ], False),
    ]

    responses = [
        ("STI_vs_CoastMedian", "STI", "raw"),
    ]
    for bn in BAND_NAMES:
        responses.append((f"{bn}_composition", f"{bn}_composition", "logit"))
        responses.append((
            f"{bn}_log10_absolute_band_power",
            f"{bn}_apparent_abs_power",
            "raw"
        ))

    model_rows = []

    for source_col, response_name, scale in responses:
        work = df.copy()
        if scale == "logit":
            work["response_y"] = logit(work[source_col].to_numpy(float))
        else:
            work["response_y"] = pd.to_numeric(work[source_col], errors="coerce")

        for model_name, xcols, do_cv in model_defs:
            row = fit_model(
                work,
                "response_y",
                xcols,
                model_name,
                do_cv=do_cv
            )
            if row is None:
                continue

            row.update({
                "response": response_name,
                "source_column": source_col,
                "response_scale": scale,
                "LOSO_interpretation": (
                    "standard leave-one-site-out"
                    if do_cv
                    else "not used as primary evidence because Lagoon indicator is group-imbalanced"
                )
            })
            model_rows.append(row)

    models = pd.DataFrame(model_rows)

    # best by AICc
    best_rows = []
    for response, sub in models.groupby("response"):
        s = sub.dropna(subset=["AICc"]).sort_values("AICc")
        if len(s):
            best_rows.append(s.iloc[0].to_dict())
    best = pd.DataFrame(best_rows)

    # Site summary
    site_summary = (
        df.groupby("site")
          .agg(
              n=("channel", "size"),
              STI_median=("STI_vs_CoastMedian", "median"),
              STI_q25=("STI_vs_CoastMedian", lambda x: np.nanpercentile(x, 25)),
              STI_q75=("STI_vs_CoastMedian", lambda x: np.nanpercentile(x, 75)),
              cosine_to_Coast_median=("composition_cosine_vs_CoastMedian", "median"),
              euclidean_distance_median_m=("distance_to_coast_polyline_m", "median"),
              network_distance_median_m=("network_distance_to_CH330_m", "median"),
          )
          .reset_index()
    )

    outdir = root / STI_OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(outdir / "STI_and_geometry_model_results_v3.xlsx", engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="channel_metrics", index=False)
        site_summary.to_excel(writer, sheet_name="site_summary", index=False)
        models.to_excel(writer, sheet_name="all_models", index=False)
        best.to_excel(writer, sheet_name="best_models_by_AICc", index=False)
        geom_meta.to_excel(writer, sheet_name="geometry_metadata", index=False)
        pd.DataFrame({
            "parameter": [
                "STI_definition",
                "Coast_reference",
                "primary_distance",
                "alternative_distance",
                "orientation_predictor",
                "lagoon_effect_model_note",
            ],
            "value": [
                "0.5 * sum(|C_i(channel) - C_i(Coast median)|)",
                "median four-band composition vector of all Coast channels",
                "distance_to_coast_polyline_m",
                "network_distance_to_CH330_m",
                "log(cos^2(angle to Coast normal))",
                "M6 includes inside_lagoon but LOSO is not primary evidence",
            ]
        }).to_excel(writer, sheet_name="methods", index=False)


    print("=" * 84)
    print("09 v3 完成")
    print("输出：", outdir)
    print("\nSite summary:")
    print(site_summary.to_string(index=False))
    print("\nBest models by AICc:")
    if len(best):
        print(best[["response", "model", "R2", "AICc", "LOSO_RMSE"]].to_string(index=False))
    print("=" * 84)
    return outdir / "STI_and_geometry_model_results_v3.xlsx"


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
13 | 20-m 空间分块 + Block bootstrap（v4）
=========================================

目的
----
解决一个重要统计问题：
2-m 相邻 DAS 通道具有明显空间相关性，327 个通道不能被当成
327 个完全独立样本。

本脚本把每条分析光纤按 20 m 分块：
    20 m / 2 m = 10 个通道间隔左右

每个 block 用中位数代表，然后：
1. 重新统计 STI / composition / correlation；
2. 在 block 层面重新拟合距离与方向模型；
3. 做 stratified spatial-block bootstrap，给 site median 和模型改进量置信区间。

这样论文中的统计样本从“逐通道伪重复”变成“空间块级证据”。

输出
----
13_spatial_block_analysis_v4_output/
    spatial_block_analysis_v4.xlsx
"""

from pathlib import Path
import numpy as np
import pandas as pd

from das_rbr_1min_common import (
    SITES, SITE_SHORT, choose_root_folder
)
from mechanism_common_v3 import logit, fisher_z, ols_fit, leave_one_site_out_rmse

OUTDIR = "13_spatial_block_analysis_v4_output"
BLOCK_LENGTH_M = 20.0
N_BOOT = 1000
RANDOM_SEED = 20260823


def find_file(root, filename):
    hits = list(root.rglob(filename))
    if not hits:
        raise FileNotFoundError(f"未找到 {filename}")
    return hits[0]


def make_corr_wide(transfer):
    parts = []
    for band in ["Diurnal", "Semidiurnal"]:
        sub = transfer[transfer["band"] == band][[
            "site", "channel",
            "corr0_vs_tide",
            "corrpeak_vs_tide",
            "corr0_vs_CoastDAS",
            "corrpeak_vs_CoastDAS",
            "lag_at_corrpeak_vs_tide_h",
            "lag_at_corrpeak_vs_CoastDAS_h",
            "high_quality_tide_lag",
            "high_quality_coast_lag",
        ]].copy()

        rename = {
            c: f"{band}_{c}"
            for c in sub.columns
            if c not in ["site", "channel"]
        }
        parts.append(sub.rename(columns=rename))

    out = parts[0]
    for p in parts[1:]:
        out = out.merge(p, on=["site", "channel"], how="outer")
    return out


def fit_block_model(df, ycol, xcols, transform="raw"):
    sub = df.copy()

    if transform == "logit":
        sub["response_y"] = logit(sub[ycol].to_numpy(float))
    elif transform == "fisherz":
        sub["response_y"] = fisher_z(sub[ycol].to_numpy(float))
    else:
        sub["response_y"] = pd.to_numeric(sub[ycol], errors="coerce")

    sub = sub.dropna(subset=["response_y"] + list(xcols))
    if len(sub) < max(8, len(xcols) + 3):
        return None

    X = (
        sub[list(xcols)].to_numpy(float)
        if xcols else
        np.empty((len(sub), 0))
    )

    fit = ols_fit(
        sub["response_y"].to_numpy(float),
        X,
        list(xcols)
    )
    if fit is None:
        return None

    coef = dict(zip(fit["coef_names"], fit["beta"]))
    cv_rmse, cv_n = leave_one_site_out_rmse(
        sub,
        "response_y",
        list(xcols)
    )

    return {
        "n_blocks": fit["n"],
        "R2": fit["r2"],
        "AICc": fit["aicc"],
        "LOSO_RMSE": cv_rmse,
        "LOSO_n": cv_n,
        **{f"coef_{k}": v for k, v in coef.items()},
    }


def stratified_block_bootstrap(blocks, rng):
    samples = []
    for site in SITES:
        sub = blocks[blocks["site"] == site]
        if len(sub) == 0:
            continue
        idx = rng.integers(0, len(sub), size=len(sub))
        boot = sub.iloc[idx].copy()
        # Keep origin block ID for traceability but create bootstrap-row ID.
        boot["bootstrap_site"] = site
        samples.append(boot)
    return pd.concat(samples, ignore_index=True)


def percentile_ci(x, q=(2.5, 50, 97.5)):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return (np.nan,) * len(q)
    return tuple(np.percentile(x, q))


def main(root=None):
    root = Path(root) if root is not None else choose_root_folder("选择包含 06 / 07 输出的总文件夹")

    p07 = find_file(root, "band_limited_transfer_metrics.xlsx")
    p09 = compute_sti(root)  # 先计算STI，再使用本次输出做空间分块。

    transfer = pd.read_excel(p07, sheet_name="channel_metrics")
    ch = pd.read_excel(p09, sheet_name="channel_metrics")
    corr = make_corr_wide(transfer)

    df = ch.merge(
        corr,
        on=["site", "channel"],
        how="left"
    )

    if "site_chainage_m" not in df.columns:
        raise RuntimeError(
            "09 输出中没有 site_chainage_m。请确认使用 v3 的 09 输出。"
        )

    # 20-m spatial blocks within each site.
    df["block_index"] = np.floor(
        pd.to_numeric(df["site_chainage_m"], errors="coerce")
        / BLOCK_LENGTH_M
    ).astype("Int64")

    # Columns to aggregate by median.
    numeric_cols = [
        "channel",
        "site_chainage_m",
        "x_m", "y_m",
        "distance_to_coast_polyline_m",
        "network_distance_to_CH330_m",
        "fiber_azimuth_deg",
        "acute_angle_to_coast_axis_deg",
        "acute_angle_to_coast_normal_deg",
        "orientation_cos2_to_coast_axis",
        "orientation_cos2_to_coast_normal",
        "STI_vs_CoastMedian",
        "composition_cosine_vs_CoastMedian",
        "Diurnal_composition",
        "Semidiurnal_composition",
        "Terdiurnal_composition",
        "Quarter-diurnal_composition",
        "Diurnal_log10_absolute_band_power",
        "Semidiurnal_log10_absolute_band_power",
        "Diurnal_corrpeak_vs_tide",
        "Semidiurnal_corrpeak_vs_tide",
        "Diurnal_corrpeak_vs_CoastDAS",
        "Semidiurnal_corrpeak_vs_CoastDAS",
        "Diurnal_lag_at_corrpeak_vs_tide_h",
        "Semidiurnal_lag_at_corrpeak_vs_tide_h",
        "Diurnal_lag_at_corrpeak_vs_CoastDAS_h",
        "Semidiurnal_lag_at_corrpeak_vs_CoastDAS_h",
    ]
    numeric_cols = [c for c in numeric_cols if c in df.columns]

    agg_dict = {
        c: "median" for c in numeric_cols
    }
    agg_dict.update({
        "channel": ["min", "max", "median"],
    })

    grouped = (
        df.groupby(["site", "block_index"], dropna=True)
          .agg(agg_dict)
    )

    # Flatten MultiIndex columns.
    grouped.columns = [
        "_".join([str(x) for x in col if str(x)])
        if isinstance(col, tuple)
        else str(col)
        for col in grouped.columns
    ]
    blocks = grouped.reset_index()

    # Normalize expected names after mixed aggregation.
    rename = {}
    for c in numeric_cols:
        if c == "channel":
            continue
        candidate = f"{c}_median"
        if candidate in blocks.columns:
            rename[candidate] = c

    rename.update({
        "channel_min": "channel_min",
        "channel_max": "channel_max",
        "channel_median": "channel_center",
    })
    blocks = blocks.rename(columns=rename)

    # Channel count per block.
    counts = (
        df.groupby(["site", "block_index"])
          .size()
          .reset_index(name="n_channels")
    )
    blocks = blocks.merge(
        counts,
        on=["site", "block_index"],
        how="left"
    )

    blocks["block_label"] = (
        blocks["site"]
        + "_B"
        + blocks["block_index"].astype(int).astype(str)
    )

    # Log orientation predictor.
    blocks["log_orientation_to_coast_normal"] = np.log(
        np.clip(
            blocks["orientation_cos2_to_coast_normal"].to_numpy(float),
            1e-4,
            1.0
        )
    )

    # Site block summaries.
    site_summary_rows = []
    for site in SITES:
        sub = blocks[blocks["site"] == site]
        if len(sub) == 0:
            continue

        row = {
            "site": site,
            "n_blocks": len(sub),
            "STI_median": np.nanmedian(sub["STI_vs_CoastMedian"]),
            "D_corr_Coast_median": np.nanmedian(sub["Diurnal_corrpeak_vs_CoastDAS"]),
            "SD_corr_Coast_median": np.nanmedian(sub["Semidiurnal_corrpeak_vs_CoastDAS"]),
        }
        site_summary_rows.append(row)

    site_summary = pd.DataFrame(site_summary_rows)

    # Block-level model comparison.
    response_defs = [
        ("STI_vs_CoastMedian", "STI", "raw"),
        ("Diurnal_composition", "Diurnal composition", "logit"),
        ("Semidiurnal_composition", "Semidiurnal composition", "logit"),
        ("Diurnal_corrpeak_vs_CoastDAS", "Diurnal peak corr vs Coast DAS", "fisherz"),
        ("Semidiurnal_corrpeak_vs_CoastDAS", "Semidiurnal peak corr vs Coast DAS", "fisherz"),
    ]

    model_defs = [
        ("B0_intercept", []),
        ("B1_euclidean", ["distance_to_coast_polyline_m"]),
        ("B2_network", ["network_distance_to_CH330_m"]),
        ("B3_orientation", ["log_orientation_to_coast_normal"]),
        ("B4_euclidean_orientation", [
            "distance_to_coast_polyline_m",
            "log_orientation_to_coast_normal"
        ]),
        ("B5_network_orientation", [
            "network_distance_to_CH330_m",
            "log_orientation_to_coast_normal"
        ]),
    ]

    model_rows = []
    for ycol, label, transform in response_defs:
        for model_name, xcols in model_defs:
            res = fit_block_model(
                blocks,
                ycol,
                xcols,
                transform=transform
            )
            if res is None:
                continue
            model_rows.append({
                "response": label,
                "source_column": ycol,
                "transform": transform,
                "model": model_name,
                "predictors": ", ".join(xcols) if xcols else "(intercept only)",
                **res
            })

    model_results = pd.DataFrame(model_rows)

    # Bootstrap.
    rng = np.random.default_rng(RANDOM_SEED)
    boot_rows = []

    for ib in range(N_BOOT):
        boot = stratified_block_bootstrap(blocks, rng)

        # site medians
        for site in SITES:
            sub = boot[boot["site"] == site]
            if len(sub) == 0:
                continue
            boot_rows.append({
                "iteration": ib,
                "result_type": "site_median",
                "site": site,
                "metric": "STI",
                "value": np.nanmedian(sub["STI_vs_CoastMedian"]),
            })
            boot_rows.append({
                "iteration": ib,
                "result_type": "site_median",
                "site": site,
                "metric": "D_corr_Coast",
                "value": np.nanmedian(sub["Diurnal_corrpeak_vs_CoastDAS"]),
            })
            boot_rows.append({
                "iteration": ib,
                "result_type": "site_median",
                "site": site,
                "metric": "SD_corr_Coast",
                "value": np.nanmedian(sub["Semidiurnal_corrpeak_vs_CoastDAS"]),
            })

        # ΔR²: distance+orientation minus distance-only.
        for ycol, label, transform in response_defs:
            r1 = fit_block_model(
                boot,
                ycol,
                ["distance_to_coast_polyline_m"],
                transform=transform
            )
            r2 = fit_block_model(
                boot,
                ycol,
                [
                    "distance_to_coast_polyline_m",
                    "log_orientation_to_coast_normal"
                ],
                transform=transform
            )
            if r1 is not None and r2 is not None:
                boot_rows.append({
                    "iteration": ib,
                    "result_type": "delta_R2",
                    "site": "",
                    "metric": label,
                    "value": r2["R2"] - r1["R2"],
                })

    boot_df = pd.DataFrame(boot_rows)

    # Bootstrap summaries.
    boot_summary_rows = []

    site_boot = boot_df[
        boot_df["result_type"] == "site_median"
    ]
    for (site, metric), sub in site_boot.groupby(["site", "metric"]):
        lo, med, hi = percentile_ci(sub["value"])
        boot_summary_rows.append({
            "result_type": "site_median",
            "site": site,
            "metric": metric,
            "q025": lo,
            "median": med,
            "q975": hi,
        })

    dr_boot = boot_df[
        boot_df["result_type"] == "delta_R2"
    ]
    for metric, sub in dr_boot.groupby("metric"):
        lo, med, hi = percentile_ci(sub["value"])
        boot_summary_rows.append({
            "result_type": "delta_R2",
            "site": "",
            "metric": metric,
            "q025": lo,
            "median": med,
            "q975": hi,
        })

    boot_summary = pd.DataFrame(boot_summary_rows)

    outdir = root / OUTDIR
    outdir.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(
        outdir / "spatial_block_analysis_v4.xlsx",
        engine="openpyxl"
    ) as writer:
        df.to_excel(writer, sheet_name="channel_input", index=False)
        blocks.to_excel(writer, sheet_name="20m_blocks", index=False)
        site_summary.to_excel(writer, sheet_name="site_block_summary", index=False)
        model_results.to_excel(writer, sheet_name="block_models", index=False)
        boot_summary.to_excel(writer, sheet_name="bootstrap_summary", index=False)
        # 1000 boot iterations can be large but still manageable.
        boot_df.to_excel(writer, sheet_name="bootstrap_all", index=False)
        pd.DataFrame({
            "parameter": [
                "block_length_m",
                "channel_spacing_m",
                "bootstrap_iterations",
                "bootstrap_scheme",
                "purpose",
            ],
            "value": [
                BLOCK_LENGTH_M,
                2.0,
                N_BOOT,
                "stratified resampling of 20-m blocks within each site",
                "reduce spatial pseudoreplication and quantify block-level uncertainty",
            ]
        }).to_excel(writer, sheet_name="metadata", index=False)


    print("=" * 84)
    print("13 v4 完成：20-m spatial blocks + block bootstrap")
    print("输出：", outdir)
    print(f"通道级样本数 = {len(df)}")
    print(f"20-m block 数 = {len(blocks)}")
    print("\nBlock-level site summary:")
    print(site_summary.to_string(index=False))
    print("\nBootstrap summary:")
    print(boot_summary.to_string(index=False))
    print("=" * 84)


if __name__ == "__main__":
    main()
