"""Render only aggregate Figures 1, 2 and 4; never read participant residuals."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from verify_public_results import audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-inputs", action="store_true")
    args = parser.parse_args()
    audit()
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "5_paper"))
    import build_figures123_v2 as first
    import build_figure4_v2 as fourth
    import matplotlib.pyplot as plt

    results = first.load_results()
    units = fourth.load_confirmatory_units()
    if args.check_inputs:
        print(f"PASS: {len(results)} aggregate tests and {len(units)} reported follow-up units")
        return
    output = root / "outputs" / "figures"
    output.mkdir(parents=True, exist_ok=True)
    first.configure_style()
    for name, builder in (("figure1", first.build_figure1), ("figure2", first.build_figure2)):
        figure = builder(results)
        try:
            figure.savefig(output / f"{name}.png", dpi=300, facecolor="white")
        finally:
            plt.close(figure)
    fourth.configure_style()
    figure = fourth.build_figure(units)
    try:
        figure.savefig(output / "figure4.png", dpi=300, facecolor="white")
    finally:
        plt.close(figure)
    print(f"Wrote aggregate Figure 1, 2 and 4 previews to {output}")


if __name__ == "__main__":
    main()
