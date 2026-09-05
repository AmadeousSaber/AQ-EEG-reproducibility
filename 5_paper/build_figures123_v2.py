"""Build revised manuscript Figures 1–3 from immutable statistical outputs.

The full entry point needs private individual residuals for Figure 3.
For the public Figure 1/2/4 workflow, install requirements-figures.txt and
run scripts/render_public_figures.py from the repository root.
"""

from __future__ import annotations

import hashlib
import io
from pathlib import Path
import tempfile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent / "reference_outputs" / "4_corr"
OUTPUT_DIR = HERE / "figures_v2"

AQ_ORDER = ["communication", "switching", "detail", "imagination", "social", "total"]
STAR = "*"
OPEN_COLOR = "#4C78A8"
CLOSED_COLOR = "#E28E5B"
EXPECTED_SOURCE_MANIFEST_SHA256 = "833378d520c199f13acf3ee29b843e8233284ec00f73fa0933995a32c05ce5b2"


def configure_style() -> None:
    """Set publication-safe fonts and deterministic vector identifiers."""

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
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


def source_manifest_hash() -> str:
    """Hash the exact source CSV bytes used by Figures 1–3."""

    paths = [
        *(DATA_DIR / f"results_{level}_{cond}_{conn}.csv"
          for cond in ("open", "closed")
          for conn in ("dwpli", "aec")
          for level in ("global", "network", "node")),
        DATA_DIR / "plot_data_all.csv",
    ]
    digest = hashlib.sha256()
    for path in sorted(paths):
        digest.update(path.name.encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def load_results() -> pd.DataFrame:
    """Load all global, network-mean, and nodal tables and validate hit count."""

    frames: list[pd.DataFrame] = []
    for condition in ("open", "closed"):
        for estimator in ("dwpli", "aec"):
            for level in ("global", "network", "node"):
                table = pd.read_csv(DATA_DIR / f"results_{level}_{condition}_{estimator}.csv")
                table["cond"] = condition
                table["conn"] = estimator
                table["level"] = level
                frames.append(table)
    results = pd.concat(frames, ignore_index=True)
    if int((results["q4"] < 0.05).sum()) != 14:
        raise ValueError("Expected exactly 14 locally corrected tests")
    return results


def save_bundle(fig: plt.Figure, name: str) -> None:
    """Write a validated, atomic SVG/PDF/PNG/TIFF export bundle."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"{name}_", dir=OUTPUT_DIR) as temp_dir:
        stem = Path(temp_dir) / name
        fig.savefig(
            stem.with_suffix(".svg"),
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
            metadata={"Date": None},
        )
        fig.savefig(
            stem.with_suffix(".pdf"),
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
            metadata={"CreationDate": None, "ModDate": None, "Creator": "Matplotlib"},
        )
        fig.savefig(
            stem.with_suffix(".png"),
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
        )

        buffer = io.BytesIO()
        fig.savefig(
            buffer,
            format="tiff",
            dpi=600,
            bbox_inches="tight",
            facecolor="white",
            transparent=False,
        )
        buffer.seek(0)
        with Image.open(buffer) as source_image:
            rgb_image = source_image.convert("RGB")
            rgb_image.load()
        buffer.close()
        try:
            rgb_image.save(stem.with_suffix(".tiff"), compression="tiff_lzw", dpi=(600, 600))
        finally:
            rgb_image.close()

        for extension in (".svg", ".pdf", ".png", ".tiff"):
            path = stem.with_suffix(extension)
            if not path.exists() or path.stat().st_size == 0:
                raise RuntimeError(f"Missing or empty figure export: {path}")
        with Image.open(stem.with_suffix(".tiff")) as check_tiff:
            if check_tiff.mode != "RGB" or check_tiff.info.get("dpi") != (600.0, 600.0):
                raise RuntimeError("TIFF export failed RGB/600-dpi validation")

        for extension in (".svg", ".pdf", ".png", ".tiff"):
            stem.with_suffix(extension).replace(OUTPUT_DIR / f"{name}{extension}")


def select_row(
    results: pd.DataFrame,
    condition: str,
    estimator: str,
    level: str,
    selection: dict[str, object],
    aq_score: str,
) -> pd.Series:
    """Select one uniquely identified row from the frozen results."""

    rows = results[
        results["cond"].eq(condition)
        & results["conn"].eq(estimator)
        & results["level"].eq(level)
    ]
    for column, value in selection.items():
        rows = rows[rows[column].eq(value)]
    rows = rows[rows["aqvar"].eq(aq_score)]
    if len(rows) != 1:
        raise ValueError(
            f"Expected one row for {condition}/{estimator}/{level}/{selection}/{aq_score}, "
            f"found {len(rows)}"
        )
    return rows.iloc[0]


def build_figure1(results: pd.DataFrame) -> plt.Figure:
    """Show selected measures across every AQ score with local-family stars."""

    columns = [
        ("Default label\nα clustering", "dwpli", "node", {"node": 42, "band": "alpha", "quantity": "auc_cc"}),
        ("Control mean\nβ degree", "aec", "network", {"network": "Cont", "band": "beta", "quantity": "auc_degree"}),
        ("Control label\nβ degree", "aec", "node", {"node": 82, "band": "beta", "quantity": "auc_degree"}),
        ("Control label\nα degree", "dwpli", "node", {"node": 89, "band": "alpha", "quantity": "auc_degree"}),
        ("Attention label\nα clustering", "dwpli", "node", {"node": 19, "band": "alpha", "quantity": "auc_cc"}),
        ("Visual mean\nβ degree", "aec", "network", {"network": "Vis", "band": "beta", "quantity": "auc_degree"}),
        ("Visual label\nβ degree", "aec", "node", {"node": 55, "band": "beta", "quantity": "auc_degree"}),
        ("Somatomotor label\nβ degree", "aec", "node", {"node": 65, "band": "beta", "quantity": "auc_degree"}),
        ("Global γ\nclustering", "aec", "global", {"metric": "aec_gamma_auc_gcc"}),
        ("Global γ\nsmall-worldness", "aec", "global", {"metric": "aec_gamma_auc_sw"}),
        ("Attention label\nγ clustering", "dwpli", "node", {"node": 23, "band": "gamma", "quantity": "auc_cc"}),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15), sharey=True)
    for axis, condition in zip(axes, ("open", "closed"), strict=True):
        values = np.full((len(AQ_ORDER), len(columns)), np.nan)
        stars = np.full((len(AQ_ORDER), len(columns)), "", dtype=object)
        for column_index, (_, estimator, level, selection) in enumerate(columns):
            for row_index, aq_score in enumerate(AQ_ORDER):
                row = select_row(results, condition, estimator, level, selection, aq_score)
                values[row_index, column_index] = row["r"]
                if row["q4"] < 0.05:
                    stars[row_index, column_index] = STAR
        axis.imshow(values, cmap="RdBu_r", vmin=-0.65, vmax=0.65, aspect="auto")
        for row_index in range(len(AQ_ORDER)):
            for column_index in range(len(columns)):
                axis.text(
                    column_index,
                    row_index,
                    f"{values[row_index, column_index]:+.2f}{stars[row_index, column_index]}",
                    ha="center",
                    va="center",
                    fontsize=5.3,
                    color="black",
                )
        axis.set_xticks(range(len(columns)))
        axis.set_xticklabels([column[0] for column in columns], rotation=47, ha="right", fontsize=5.7)
        axis.set_title(f"eyes-{condition}")
        axis.set_yticks(range(len(AQ_ORDER)))
        axis.set_yticklabels(AQ_ORDER, fontsize=6.5)

    fig.suptitle(
        "Selected locally corrected tests across AQ scores (partial r; * local q < .05)",
        fontsize=8.5,
    )
    fig.tight_layout(rect=[0, 0, 0.93, 0.95])
    color_axis = fig.add_axes([0.945, 0.20, 0.012, 0.62])
    scalar = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-0.65, 0.65))
    fig.colorbar(scalar, cax=color_axis, label="partial r")
    return fig


def build_figure2(results: pd.DataFrame) -> plt.Figure:
    """Plot the overlapping attention-to-detail tests without replication language."""

    rows = [
        ("open", "Visual-label mean β degree", "aec", "network", {"network": "Vis", "band": "beta", "quantity": "auc_degree"}),
        ("open", "Visual label 1 β degree", "aec", "node", {"node": 55, "band": "beta", "quantity": "auc_degree"}),
        ("open", "Visual label 2 β degree", "aec", "node", {"node": 58, "band": "beta", "quantity": "auc_degree"}),
        ("open", "Visual label 3 β degree", "aec", "node", {"node": 54, "band": "beta", "quantity": "auc_degree"}),
        ("open", "Visual label 4 β degree", "aec", "node", {"node": 5, "band": "beta", "quantity": "auc_degree"}),
        ("closed", "Global γ clustering", "aec", "global", {"metric": "aec_gamma_auc_gcc"}),
        ("closed", "Global γ small-worldness", "aec", "global", {"metric": "aec_gamma_auc_sw"}),
        ("closed", "Somatomotor label β degree", "aec", "node", {"node": 65, "band": "beta", "quantity": "auc_degree"}),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    for axis, condition in zip(axes, ("open", "closed"), strict=True):
        labels: list[str] = []
        values: list[float] = []
        lower_errors: list[float] = []
        upper_errors: list[float] = []
        for row_condition, label, estimator, level, selection in rows:
            if row_condition != condition:
                continue
            row = select_row(results, condition, estimator, level, selection, "detail")
            labels.append(label)
            values.append(float(row["r"]))
            lower_errors.append(float(row["r"] - row["ci_lo"]))
            upper_errors.append(float(row["ci_hi"] - row["r"]))
        y_positions = np.arange(len(labels))[::-1]
        axis.barh(
            y_positions,
            values,
            xerr=[lower_errors, upper_errors],
            color=OPEN_COLOR if condition == "open" else CLOSED_COLOR,
            alpha=0.9,
            capsize=2.5,
            error_kw={"elinewidth": 0.8, "capthick": 0.8},
        )
        axis.set_yticks(y_positions, labels=labels)
        axis.axvline(0, color="black", linewidth=0.7)
        axis.set_xlabel("partial r with attention to detail")
        axis.set_title(f"eyes-{condition}")
        for y_position, value in zip(y_positions, values, strict=True):
            axis.text(
                value + np.sign(value) * 0.025,
                y_position,
                STAR,
                va="center",
                ha="left" if value > 0 else "right",
                fontsize=9,
            )
    fig.suptitle("Overlapping attention-to-detail tests (* local q < .05)", fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


def build_figure3() -> plt.Figure:
    """Plot six representative age- and sex-residualized associations."""

    data = pd.read_csv(DATA_DIR / "plot_data_all.csv")
    picks = [
        ("open", "dwpli", "node", "communication", "default label α clustering × communication"),
        ("open", "aec", "network", "detail", "visual-label mean β degree × attention to detail"),
        ("open", "aec", "node", "switching", "control label β degree × attention switching"),
        ("closed", "dwpli", "node", "switching", "control label α degree × attention switching"),
        ("closed", "aec", "global", "detail", "global γ clustering × attention to detail"),
        ("closed", "dwpli", "node", "imagination", "attention label γ clustering × imagination"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.2))
    used: set[object] = set()
    for axis, (condition, estimator, level, aq_score, display_label) in zip(
        axes.flat, picks, strict=True
    ):
        candidates = data[
            data["cond"].eq(condition)
            & data["conn"].eq(estimator)
            & data["level"].eq(level)
            & data["aqvar"].eq(aq_score)
            & ~data["hit_id"].isin(used)
        ]
        if candidates.empty:
            raise ValueError(f"No residual data for {condition}/{estimator}/{level}/{aq_score}")
        if level == "node" and candidates["hit_id"].nunique() > 1:
            selected_hit = candidates.sort_values("q4").iloc[0]["hit_id"]
            candidates = candidates[candidates["hit_id"].eq(selected_hit)]
        hit_id = candidates.iloc[0]["hit_id"]
        used.add(hit_id)
        group = candidates[candidates["hit_id"].eq(hit_id)]
        for gender, marker, color, label in (
            (1, "o", OPEN_COLOR, "M"),
            (2, "s", CLOSED_COLOR, "F"),
        ):
            subset = group[group["gender"].eq(gender)]
            axis.scatter(
                subset["resid_brain"],
                subset["resid_aq"],
                marker=marker,
                c=color,
                s=15,
                alpha=0.82,
                label=label,
            )
        slope, intercept = np.polyfit(group["resid_brain"], group["resid_aq"], 1)
        x_values = np.linspace(group["resid_brain"].min(), group["resid_brain"].max(), 50)
        axis.plot(x_values, slope * x_values + intercept, color="black", linewidth=0.8)
        first = group.iloc[0]
        axis.set_title(
            f"{condition} | {display_label}\nr={first['r']:+.2f}, local q={first['q4']:.3f}",
            fontsize=5.6,
        )
        axis.set_xlabel("brain residual", fontsize=6)
        axis.set_ylabel("AQ residual", fontsize=6)
        axis.tick_params(labelsize=5.6)
    axes.flat[0].legend(fontsize=6, loc="upper left")
    fig.suptitle("Representative partial-residual associations (selected local tests)", fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


def main() -> None:
    configure_style()
    observed_manifest_hash = source_manifest_hash()
    if observed_manifest_hash != EXPECTED_SOURCE_MANIFEST_SHA256:
        raise ValueError(
            "Frozen Figure 1–3 source manifest changed: "
            f"expected {EXPECTED_SOURCE_MANIFEST_SHA256}, observed {observed_manifest_hash}"
        )
    results = load_results()
    figures = [
        (build_figure1(results), "fig1_selected_local_tests_v2"),
        (build_figure2(results), "fig2_attention_detail_local_tests_v2"),
        (build_figure3(), "fig3_representative_residuals_v2"),
    ]
    try:
        for figure, name in figures:
            save_bundle(figure, name)
    finally:
        for figure, _ in figures:
            plt.close(figure)
    print(f"Wrote revised Figures 1–3 to {OUTPUT_DIR}")
    print(f"Source manifest SHA-256: {observed_manifest_hash}")


if __name__ == "__main__":
    main()
