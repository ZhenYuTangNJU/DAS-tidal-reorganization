#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
DAS 1-min RMS representative channel extraction

Extract:
- center channels:
    Coast Ch413
    Lagoon1 Ch67
    Lagoon2 Ch130
    Lagoon3 Ch186

- neighbor envelope:
    Coast 408-418
    Lagoon1 62-72
    Lagoon2 125-135
    Lagoon3 181-191

Output:
    DAS_curve_data.xlsx
"""

from pathlib import Path
import pandas as pd
import tkinter as tk
from tkinter import filedialog


CENTER_CHANNELS = {
    "Coast": 413,
    "Lagoon1": 67,
    "Lagoon2": 130,
    "Lagoon3": 186,
}

NEIGHBOR_WINDOWS = {
    "Coast": (408, 418),
    "Lagoon1": (62, 72),
    "Lagoon2": (125, 135),
    "Lagoon3": (181, 191),
}

INPUT_FILES = {
    "Coast": "Coast_1min_RMS.xlsx",
    "Lagoon1": "Lagoon1_1min_RMS.xlsx",
    "Lagoon2": "Lagoon2_1min_RMS.xlsx",
    "Lagoon3": "Lagoon3_1min_RMS.xlsx",
}


def select_folder():
    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    folder = filedialog.askdirectory(
        title="选择 DAS 1-min RMS Excel 文件所在文件夹"
    )
    root.destroy()

    if not folder:
        raise RuntimeError("未选择文件夹")

    return Path(folder)


def read_das_excel(filepath):
    df = pd.read_excel(
        filepath,
        sheet_name="RMS_1min",
        index_col=0
    )

    df.index = pd.to_numeric(df.index, errors="coerce").astype(int)
    df.columns = pd.to_datetime(df.columns, errors="coerce")
    df = df.loc[:, ~df.columns.isna()]
    df = df.sort_index(axis=1)

    return df


def extract_site(df, site):

    center_ch = CENTER_CHANNELS[site]
    c0, c1 = NEIGHBOR_WINDOWS[site]

    if center_ch not in df.index:
        raise ValueError(f"{site} 缺少 Ch {center_ch}")

    neighbors = [
        ch for ch in df.index
        if c0 <= ch <= c1
    ]

    center = df.loc[center_ch].astype(float)
    neighbor = df.loc[neighbors].astype(float)

    mean = neighbor.mean(axis=0)
    std = neighbor.std(axis=0)

    return pd.DataFrame({
        "time": df.columns,
        f"{site}_Ch{center_ch}": center.values,
        f"{site}_neighbor_mean": mean.values,
        f"{site}_neighbor_std": std.values,
        f"{site}_upper": (mean + std).values,
        f"{site}_lower": (mean - std).values,
    })


def main():

    folder = select_folder()

    results = []

    for site, filename in INPUT_FILES.items():

        path = folder / filename

        if not path.exists():
            raise FileNotFoundError(path)

        df = read_das_excel(path)

        results.append(
            extract_site(df, site)
        )

    out = results[0]

    for df in results[1:]:
        out = out.merge(
            df,
            on="time",
            how="outer"
        )

    out = out.sort_values("time")

    outfile = folder / "DAS_curve_data.xlsx"

    out.to_excel(
        outfile,
        index=False
    )

    print("完成:", outfile)


if __name__ == "__main__":
    main()
