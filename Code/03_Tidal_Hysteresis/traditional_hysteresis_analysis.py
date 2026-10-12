from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2] / "03_Inputs_and_Intermediate_Results"
INPUT = ROOT / "15_NEW_hysteresis_metrics_output" / "hysteresis_metrics_v2_upgraded.xlsx"
OUTDIR = Path(__file__).resolve().parents[2] / "Reprocessed_Outputs" / "Traditional_Hysteresis"
SEED = 20260924
N_BOOT = 10_000

SITE_ORDER = ["Coast", "Lagoon1", "Lagoon2", "Lagoon3"]
SITE_LABEL = {
    "Coast": "Coast",
    "Lagoon1": "Western lagoon",
    "Lagoon2": "Northern lagoon",
    "Lagoon3": "Eastern lagoon",
}


def integrate_adjacent(x: np.ndarray, y: np.ndarray, nominal_dx: float) -> tuple[float, float]:
    """Trapezoidal integral over adjacent valid bins only; returns area and covered width."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    x = x[good]
    y = y[good]
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    area = 0.0
    width = 0.0
    for i in range(len(x) - 1):
        dx = x[i + 1] - x[i]
        if dx <= 0 or dx > 1.5 * nominal_dx:
            continue
        area += 0.5 * (y[i] + y[i + 1]) * dx
        width += dx
    return float(area), float(width)


def channel_loop_metrics(bins: pd.DataFrame) -> pd.DataFrame:
    levels = np.sort(bins["relative_level"].dropna().unique())
    nominal_dx = float(np.nanmedian(np.diff(levels)))
    rows: list[dict] = []
    for (site, channel), g in bins.groupby(["site", "channel"], sort=False):
        g = g.sort_values("relative_level")
        x = g["relative_level"].to_numpy(float)
        delta = g["delta_rf"].to_numpy(float)
        signed, width = integrate_adjacent(x, delta, nominal_dx)
        absolute, width_abs = integrate_adjacent(x, np.abs(delta), nominal_dx)
        positive, _ = integrate_adjacent(x, np.maximum(delta, 0.0), nominal_dx)
        negative, _ = integrate_adjacent(x, np.maximum(-delta, 0.0), nominal_dx)
        if width <= 0 or width_abs <= 0:
            continue
        rows.append(
            {
                "site": site,
                "site_label": SITE_LABEL[site],
                "channel": int(channel),
                "signed_loop_area": signed,
                "absolute_loop_area": absolute,
                "positive_lobe_area": positive,
                "negative_lobe_area": negative,
                "integrated_level_width": width,
                "mean_signed_separation": signed / width,
                "mean_absolute_separation": absolute / width_abs,
                "directionality": signed / absolute if absolute > 0 else np.nan,
                "cancellation_fraction": 1.0 - abs(signed) / absolute if absolute > 0 else np.nan,
                "n_matched_bins": int(np.isfinite(delta).sum()),
            }
        )
    return pd.DataFrame(rows)


def bootstrap_median(values: np.ndarray, seed: int) -> tuple[float, float, float]:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(values), size=(N_BOOT, len(values)))
    boot = np.median(values[idx], axis=1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return float(np.median(values)), float(lo), float(hi)


def summarize(metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    fields = [
        "signed_loop_area",
        "absolute_loop_area",
        "mean_signed_separation",
        "mean_absolute_separation",
        "directionality",
        "cancellation_fraction",
    ]
    for si, site in enumerate(SITE_ORDER):
        g = metrics.loc[metrics["site"] == site]
        row = {"site": site, "site_label": SITE_LABEL[site], "n_channels": int(len(g))}
        for mi, field in enumerate(fields):
            med, lo, hi = bootstrap_median(g[field].to_numpy(float), SEED + 100 * si + mi)
            row[f"{field}_median"] = med
            row[f"{field}_lower95"] = lo
            row[f"{field}_upper95"] = hi
        row["fraction_positive_signed_area"] = float(np.mean(g["signed_loop_area"] > 0))
        rows.append(row)
    return pd.DataFrame(rows)


def site_curves_from_channels(bins: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for si, site in enumerate(SITE_ORDER):
        s = bins.loc[bins["site"] == site]
        for bi, (source_col, branch) in enumerate(
            [("rising_median", "Rising"), ("falling_median", "Falling")]
        ):
            for level, g in s.groupby("relative_level"):
                vals = g[source_col].dropna().to_numpy(float)
                med, lo, hi = bootstrap_median(vals, SEED + 1000 + 100 * si + 10 * bi + int(level * 100))
                rows.append(
                    {
                        "site": site,
                        "site_label": SITE_LABEL[site],
                        "branch": branch,
                        "relative_level": float(level),
                        "median_intensity": med,
                        "lower95": lo,
                        "upper95": hi,
                        "n_channels": int(len(vals)),
                    }
                )
    return pd.DataFrame(rows)














def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    bins = pd.read_excel(INPUT, sheet_name="broadband_binned_curves")
    metrics = channel_loop_metrics(bins)
    curves = site_curves_from_channels(bins)
    summary = summarize(metrics)

    definitions = pd.DataFrame(
        {
            "metric": [
                "signed_loop_area",
                "absolute_loop_area",
                "directionality",
                "cancellation_fraction",
            ],
            "definition": [
                "Integral of rising-minus-falling median DAS intensity over relative water level",
                "Integral of the absolute rising-minus-falling difference over relative water level",
                "signed_loop_area / absolute_loop_area; range -1 to 1",
                "1 - abs(signed_loop_area) / absolute_loop_area; range 0 to 1",
            ],
            "interpretation": [
                "Net branch orientation; positive means rising branch is higher on balance",
                "Total geometric branch separation without cancellation",
                "Consistency of loop orientation across water level",
                "Fraction of geometric area cancelled by oppositely oriented lobes",
            ],
        }
    )

    with pd.ExcelWriter(OUTDIR / "traditional_hysteresis_source_data.xlsx", engine="openpyxl") as writer:
        bins.to_excel(writer, sheet_name="channel_branch_bins", index=False)
        curves.to_excel(writer, sheet_name="site_loop_curves", index=False)
        metrics.to_excel(writer, sheet_name="channel_loop_area", index=False)
        summary.to_excel(writer, sheet_name="site_summary", index=False)
        definitions.to_excel(writer, sheet_name="metric_definitions", index=False)

    print(summary.to_string(index=False))
    print(f"Outputs written to: {OUTDIR}")


if __name__ == "__main__":
    main()
