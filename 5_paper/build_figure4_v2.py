"""Build the revised post-selection regression figure for manuscript v2.

The source table is the immutable regression summary in
``reference_outputs/4_corr``. The figure separates relative improvement over
the demographic baseline (delta q-squared) from absolute predictive
performance (q-squared) so that a positive delta cannot be mistaken for
prediction above the sample-mean benchmark.

For the public Figure 1/2/4 workflow, install requirements-figures.txt and
run scripts/render_public_figures.py from the repository root.
"""

from __future__ import annotations

import hashlib
import io
from pathlib import Path
import sys
import tempfile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import PIL
from PIL import Image


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "reference_outputs" / "4_corr" / "results_regression_summary.csv"
OUTPUT_DIR = HERE / "figures_v2"

OPEN_COLOR = "#4C78A8"
CLOSED_COLOR = "#E28E5B"
BASELINE_COLOR = "#777777"
EXPECTED_SOURCE_SHA256 = "923e7efe71d48575c788dcad541ab01d60e103038e9e73a5a27603b92bfe1858"

EXPECTED_KEYS = [
    ("open", "U1", "detail"),
    ("open", "U2", "switching"),
    ("open", "U3", "communication"),
    ("closed", "U1", "detail"),
    ("closed", "U2", "switching"),
    ("closed", "U3", "communication"),
    ("closed", "U4", "imagination"),
]
REQUIRED_COLUMNS = {
    "cond",
    "unit",
    "family",
    "aqvar",
    "dR2_main",
    "dR2_boot_lo",
    "dR2_boot_hi",
    "q2_base",
    "q2_unit",
    "dq2",
}
NUMERIC_COLUMNS = [
    "dR2_main",
    "dR2_boot_lo",
    "dR2_boot_hi",
    "q2_base",
    "q2_unit",
    "dq2",
]


def configure_style() -> None:
    """Apply a compact, journal-readable style with editable vector text."""

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "aq-manuscript-v2",
        }
    )


def load_confirmatory_units() -> pd.DataFrame:
    """Load and validate the seven confirmatory-labelled post-selection rows."""

    source_sha256 = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    if source_sha256 != EXPECTED_SOURCE_SHA256:
        raise ValueError(
            "Frozen regression summary checksum changed: "
            f"expected {EXPECTED_SOURCE_SHA256}, observed {source_sha256}"
        )

    table = pd.read_csv(SOURCE)
    missing_columns = REQUIRED_COLUMNS.difference(table.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    units = table.loc[table["family"].eq("confirm")].copy()
    if len(units) != 7:
        raise ValueError(f"Expected 7 follow-up units, found {len(units)}")

    units["unit_code"] = units["unit"].str.extract(r"^(U\d+)", expand=False)
    key_columns = ["cond", "unit_code", "aqvar"]
    if units[key_columns].isna().any(axis=None):
        raise ValueError("Condition, unit code, or AQ domain contains missing/unrecognized values")
    if units.duplicated(key_columns).any():
        raise ValueError("Duplicate follow-up unit keys detected")

    observed_keys = set(map(tuple, units[key_columns].itertuples(index=False, name=None)))
    if observed_keys != set(EXPECTED_KEYS):
        raise ValueError(
            "Unexpected follow-up unit set: "
            f"missing={sorted(set(EXPECTED_KEYS) - observed_keys)}, "
            f"extra={sorted(observed_keys - set(EXPECTED_KEYS))}"
        )

    values = units[NUMERIC_COLUMNS].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Non-finite numeric value found in a plotted field")
    if not (
        (units["dR2_boot_lo"] <= units["dR2_main"])
        & (units["dR2_main"] <= units["dR2_boot_hi"])
    ).all():
        raise ValueError("At least one bootstrap interval does not contain its dR2 estimate")
    if not np.allclose(
        units["dq2"],
        units["q2_unit"] - units["q2_base"],
        rtol=1e-12,
        atol=1e-12,
    ):
        raise ValueError("dq2 is inconsistent with q2_unit - q2_base")
    if units["dR2_boot_hi"].max() >= 0.58 or units["dR2_boot_lo"].min() < 0:
        raise ValueError("In-sample interval falls outside the configured panel range")
    q2_values = units[["q2_base", "q2_unit"]].to_numpy(dtype=float)
    if q2_values.min() <= -0.30 or q2_values.max() >= 0.047:
        raise ValueError("Predictive q2 value falls outside the configured panel range")

    order = {key: idx for idx, key in enumerate(EXPECTED_KEYS)}
    units["display_order"] = [order[key] for key in map(tuple, units[key_columns].itertuples(index=False, name=None))]
    units = units.sort_values("display_order", ascending=False).reset_index(drop=True)

    condition_labels = units["cond"].map({"open": "EO", "closed": "EC"})
    if condition_labels.isna().any():
        raise ValueError("Unexpected condition label")
    domain_labels = units["aqvar"].map(
        {
            "detail": "attention to detail",
            "switching": "attention switching",
            "communication": "communication",
            "imagination": "imagination",
        }
    )
    if domain_labels.isna().any():
        raise ValueError("Unexpected AQ domain label")
    units["label"] = condition_labels + " | " + domain_labels
    units["color"] = units["cond"].map({"open": OPEN_COLOR, "closed": CLOSED_COLOR})
    return units


def build_figure(units: pd.DataFrame) -> plt.Figure:
    """Create in-sample and nested-LOO panels from the same seven units."""

    y = np.arange(len(units))
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.2, 3.35),
        sharey=True,
        gridspec_kw={"width_ratios": [1.0, 1.15], "wspace": 0.12},
    )

    # Panel a: selection-conditioned in-sample added variance.
    ax = axes[0]
    xerr = np.vstack(
        [
            units["dR2_main"] - units["dR2_boot_lo"],
            units["dR2_boot_hi"] - units["dR2_main"],
        ]
    )
    ax.barh(
        y,
        units["dR2_main"],
        xerr=xerr,
        color=units["color"],
        alpha=0.9,
        capsize=2.5,
        error_kw={"elinewidth": 0.8, "capthick": 0.8},
    )
    ax.set_yticks(y, labels=units["label"])
    ax.set_xlim(0, 0.58)
    ax.set_xlabel(r"Added variance, $\Delta R^2$ (bootstrap 95% CI)")
    ax.set_title("Selection-conditioned in-sample fit")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.text(-0.17, 1.04, "a", transform=ax.transAxes, fontweight="bold", fontsize=9)

    # Panel b: signed nested-LOO predictive q-squared for both models.
    ax = axes[1]
    for idx, row in units.iterrows():
        ax.plot(
            [row["q2_base"], row["q2_unit"]],
            [idx, idx],
            color="#B8B8B8",
            linewidth=1.0,
            zorder=1,
        )
        ax.scatter(
            row["q2_base"],
            idx,
            s=23,
            facecolor="white",
            edgecolor=BASELINE_COLOR,
            linewidth=0.9,
            marker="o",
            zorder=2,
        )
        ax.scatter(
            row["q2_unit"],
            idx,
            s=24,
            color=row["color"],
            edgecolor="white",
            linewidth=0.5,
            marker="s",
            zorder=3,
        )
        ax.text(
            0.047,
            idx,
            rf"$\Delta q^2$={row['dq2']:+.3f}",
            va="center",
            ha="left",
            fontsize=6.1,
        )

    ax.axvline(0, color="black", linewidth=0.75)
    ax.set_xlim(-0.30, 0.16)
    ax.set_xlabel(r"Restricted-scope LOO $q^2$ (vs sample mean)")
    ax.set_title("Conditional internal validation (○ age + sex; ■ age + sex + brain)")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.text(-0.08, 1.04, "b", transform=ax.transAxes, fontweight="bold", fontsize=9)

    fig.suptitle("Post-selection follow-up: in-sample fit and restricted-scope validation", y=0.995, fontsize=9)
    fig.subplots_adjust(left=0.19, right=0.985, bottom=0.17, top=0.84)
    return fig


def save_figure(fig: plt.Figure) -> None:
    """Export editable vector files and a high-resolution manuscript preview."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_name = "fig4_post_selection_performance_v2"
    with tempfile.TemporaryDirectory(prefix="figure4_v2_", dir=OUTPUT_DIR) as temp_dir:
        temp_stem = Path(temp_dir) / output_name
        fig.savefig(
            temp_stem.with_suffix(".svg"),
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
            metadata={"Date": None},
        )
        fig.savefig(
            temp_stem.with_suffix(".pdf"),
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
            metadata={"CreationDate": None, "ModDate": None, "Creator": "Matplotlib"},
        )
        fig.savefig(
            temp_stem.with_suffix(".png"),
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
        )

        tiff_buffer = io.BytesIO()
        fig.savefig(
            tiff_buffer,
            format="tiff",
            dpi=600,
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
        )
        tiff_buffer.seek(0)
        with Image.open(tiff_buffer) as source_image:
            rgb_image = source_image.convert("RGB")
            rgb_image.load()
        tiff_buffer.close()
        try:
            rgb_image.save(
                temp_stem.with_suffix(".tiff"),
                compression="tiff_lzw",
                dpi=(600, 600),
            )
        finally:
            rgb_image.close()

        for extension in (".svg", ".pdf", ".png", ".tiff"):
            temporary_output = temp_stem.with_suffix(extension)
            if not temporary_output.exists() or temporary_output.stat().st_size == 0:
                raise RuntimeError(f"Missing or empty figure export: {temporary_output}")
        with Image.open(temp_stem.with_suffix(".tiff")) as check_tiff:
            if check_tiff.mode != "RGB" or check_tiff.info.get("dpi") != (600.0, 600.0):
                raise RuntimeError("TIFF export failed RGB/600-dpi validation")

        for extension in (".svg", ".pdf", ".png", ".tiff"):
            temp_stem.with_suffix(extension).replace(OUTPUT_DIR / f"{output_name}{extension}")


def main() -> None:
    configure_style()
    units = load_confirmatory_units()
    fig = build_figure(units)
    try:
        save_figure(fig)
    finally:
        plt.close(fig)
    source_sha256 = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    print(f"Wrote revised Figure 4 to {OUTPUT_DIR}")
    print(
        "Runtime: "
        f"Python {sys.version.split()[0]}, NumPy {np.__version__}, "
        f"pandas {pd.__version__}, Matplotlib {mpl.__version__}, Pillow {PIL.__version__}"
    )
    print(f"Source SHA-256: {source_sha256}")


if __name__ == "__main__":
    main()
