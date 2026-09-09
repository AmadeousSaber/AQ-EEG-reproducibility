"""Render the revised seven-model summary without individual-level data.

Retains the established Figure 4 composition; intervals now re-estimate the
standardization/PC1 within resamples of the fixed selected feature sets.
"""
from __future__ import annotations

import hashlib
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

import build_figure4_v2 as previous

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PUBLIC_SOURCE = ROOT / "reference_outputs" / "4_corr_v2" / "results_regression_v2_summary.csv"
SOURCE = PUBLIC_SOURCE if PUBLIC_SOURCE.is_file() else ROOT / "4_corr" / "results_regression_v2_summary.csv"
OUTPUT = HERE / "figures_v3"
EXPECTED_SOURCE_SHA256 = "4a6975e1af3cb37903484e532f8860b2d57a02161730be329d8a336c8ff9460b"


def load_units() -> pd.DataFrame:
    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    if digest != EXPECTED_SOURCE_SHA256:
        raise ValueError("Revised summary does not match the reviewed checksum")
    units = pd.read_csv(SOURCE)
    if len(units) != 7:
        raise ValueError("Expected exactly seven selected exploratory models")
    units["unit_code"] = units["unit"].str.extract(r"^(U\d+)", expand=False)
    keys = list(zip(units.cond, units.unit_code, units.aqvar))
    if keys != previous.EXPECTED_KEYS:
        raise ValueError("Unexpected model order or duplicate models")
    if not np.isfinite(units[previous.NUMERIC_COLUMNS].to_numpy(float)).all():
        raise ValueError("Nonfinite plotted data")
    if not np.allclose(units.dq2, units.q2_unit-units.q2_base, atol=1e-12, rtol=0):
        raise ValueError("Invalid leave-one-out arithmetic")
    if not ((units.dR2_boot_lo <= units.dR2_main) & (units.dR2_main <= units.dR2_boot_hi)).all():
        raise ValueError("Interval rendering requires estimates inside the intervals")
    if units.dR2_boot_hi.max() >= .58 or units.dR2_boot_lo.min() < 0:
        raise ValueError("Interval outside displayed axis range")
    if units[["q2_base","q2_unit"]].min().min() <= -.30 or units[["q2_base","q2_unit"]].max().max() >= .047:
        raise ValueError("LOO value outside displayed axis range")
    names = {"detail":"attention to detail", "switching":"attention switching",
             "communication":"communication", "imagination":"imagination"}
    units["label"] = units.cond.map({"open":"EO","closed":"EC"}) + " | " + units.aqvar.map(names)
    units["color"] = units.cond.map({"open":previous.OPEN_COLOR,"closed":previous.CLOSED_COLOR})
    return units.iloc[::-1].reset_index(drop=True)


def main() -> None:
    units = load_units()
    previous.configure_style()
    fig = previous.build_figure(units)
    try:
        fig.axes[0].set_title("Added variance in selected models")
        fig.axes[1].set_xlabel(r"Restricted-scope LOO $q^2$ (vs full-sample mean)")
        OUTPUT.mkdir(parents=True, exist_ok=True)
        stem = OUTPUT / "fig4_post_selection_performance_v3"
        for extension in ("svg","pdf","png","tiff"):
            metadata = ({"Date":None} if extension=="svg" else
                        {"CreationDate":None,"ModDate":None} if extension=="pdf" else None)
            options = {"metadata":metadata} if metadata is not None else {}
            if extension=="tiff":
                with io.BytesIO() as buffer:
                    fig.savefig(buffer, format="tiff", bbox_inches="tight", facecolor="white", dpi=600)
                    buffer.seek(0)
                    with Image.open(buffer) as source:
                        rgb = source.convert("RGB")
                        try:
                            rgb.save(stem.with_suffix(".tiff"),compression="tiff_lzw",dpi=(600,600))
                        finally:
                            rgb.close()
                continue
            fig.savefig(stem.with_suffix("."+extension), bbox_inches="tight", facecolor="white",
                        dpi=600 if extension=="tiff" else 300, **options)
        with Image.open(stem.with_suffix(".tiff")) as check:
            if check.mode != "RGB" or check.info.get("dpi") != (600.,600.):
                raise RuntimeError("TIFF must be RGB at 600 dpi")
    finally:
        plt.close(fig)
    print("Revised Figure 4 saved; source SHA-256:", EXPECTED_SOURCE_SHA256)


if __name__ == "__main__":
    main()
