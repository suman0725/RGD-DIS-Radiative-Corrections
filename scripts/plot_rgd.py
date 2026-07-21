#!/usr/bin/env python3
"""Plot RG-D DIS radiative-correction checks from the reference CSV."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TARGET_ORDER = ["LD2", "C1", "C2", "Cu", "Sn"]
Q2_EDGES = [1.0, 1.5, 2.0, 3.0, 4.5, 6.5, 10.0]
Q2_LABELS = [f"{Q2_EDGES[i]:.1f}-{Q2_EDGES[i + 1]:.1f}" for i in range(len(Q2_EDGES) - 1)]
COLORS = ["black", "red", "blue", "green", "purple", "orange"]
MARKERS = ["o", "s", "^", "D", "v", "P"]


def read_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            parsed: dict[str, object] = {
                "target": row["target"],
                "radiation_mode": row["radiation_mode"],
            }
            for key, value in row.items():
                if key in parsed:
                    continue
                parsed[key] = float(value)
            parsed["rad_to_born"] = float(parsed["sigma_rad"]) / float(parsed["sigma_born"])
            parsed["q2_bin"] = q2_bin(float(parsed["Q2"]))
            rows.append(parsed)
    return rows


def q2_bin(q2: float) -> str | None:
    for lo, hi, label in zip(Q2_EDGES[:-1], Q2_EDGES[1:], Q2_LABELS):
        if lo < q2 <= hi:
            return label
    return None


def by_target_mode(rows: list[dict[str, object]]) -> dict[tuple[str, str], list[dict[str, object]]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["target"]), str(row["radiation_mode"]))].append(row)
    for values in grouped.values():
        values.sort(key=lambda row: (float(row["Q2"]), float(row["xB"])))
    return grouped


def row_key(row: dict[str, object]) -> tuple[float, float, float]:
    return (round(float(row["Eprime"]), 4), round(float(row["theta"]), 4), round(float(row["Q2"]), 4))


def paired_modes(rows: list[dict[str, object]], target: str) -> list[tuple[dict[str, object], dict[str, object]]]:
    grouped = by_target_mode(rows)
    internal = {row_key(row): row for row in grouped[(target, "internal")]}
    external = {row_key(row): row for row in grouped[(target, "internal-external")]}
    pairs = []
    for key in sorted(internal, key=lambda item: (item[2], item[0], item[1])):
        if key in external:
            pairs.append((internal[key], external[key]))
    return pairs


def rows_for_bin(rows: list[dict[str, object]], label: str) -> list[dict[str, object]]:
    return sorted(
        [row for row in rows if row["q2_bin"] == label],
        key=lambda row: float(row["xB"]),
    )


def setup_axes(nrows: int, ncols: int, figsize: tuple[int, int]):
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True)
    return fig, list(axes.flatten())


def plot_mode_comparison(rows: list[dict[str, object]], out_dir: Path) -> Path:
    grouped = by_target_mode(rows)
    fig, axes = setup_axes(2, 3, (18, 10))
    for ax, target in zip(axes, TARGET_ORDER):
        for idx, label in enumerate(Q2_LABELS):
            color = COLORS[idx % len(COLORS)]
            marker = MARKERS[idx % len(MARKERS)]
            for mode, style, linewidth in [
                ("internal", "--", 1.2),
                ("internal-external", "-", 1.8),
            ]:
                subset = rows_for_bin(grouped[(target, mode)], label)
                if not subset:
                    continue
                ax.plot(
                    [row["xB"] for row in subset],
                    [row["rad_to_born"] for row in subset],
                    color=color,
                    marker=marker,
                    linestyle=style,
                    linewidth=linewidth,
                    markersize=5,
                    label=f"{mode} {label}",
                )
        ax.axhline(1.0, color="gray", linestyle=":", linewidth=1)
        ax.set_title(target)
        ax.set_xlabel(r"$x_B$")
        ax.set_ylabel(r"$\sigma_{rad}/\sigma_{born}$")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    axes[-1].axis("off")
    fig.suptitle("RG-D RC factor: internal-only vs internal+external")
    fig.tight_layout()
    path = out_dir / "rgd_rc_int_vs_int_ext_by_target.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_external_effect(rows: list[dict[str, object]], out_dir: Path) -> Path:
    fig, axes = setup_axes(2, 3, (18, 10))
    for ax, target in zip(axes, TARGET_ORDER):
        pairs = paired_modes(rows, target)
        plot_rows = []
        for internal, external in pairs:
            ratio = float(external["rad_to_born"]) / float(internal["rad_to_born"])
            item = dict(external)
            item["external_effect"] = ratio
            plot_rows.append(item)
        for idx, label in enumerate(Q2_LABELS):
            subset = rows_for_bin(plot_rows, label)
            if not subset:
                continue
            ax.plot(
                [row["xB"] for row in subset],
                [row["external_effect"] for row in subset],
                color=COLORS[idx % len(COLORS)],
                marker=MARKERS[idx % len(MARKERS)],
                linestyle="-",
                linewidth=1.6,
                markersize=5,
                label=label,
            )
        ax.axhline(1.0, color="gray", linestyle=":", linewidth=1)
        ax.set_title(target)
        ax.set_xlabel(r"$x_B$")
        ax.set_ylabel("external effect")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, title=r"$Q^2$")
    axes[-1].axis("off")
    fig.suptitle("External radiation effect: RC(int+ext) / RC(int)")
    fig.tight_layout()
    path = out_dir / "rgd_external_effect_by_target.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_nuclear_ratio(rows: list[dict[str, object]], out_dir: Path) -> Path:
    grouped = by_target_mode(rows)
    fig, axes = setup_axes(2, 2, (15, 10))
    targets = ["C1", "C2", "Cu", "Sn"]
    for ax, target in zip(axes, targets):
        for mode, style, linewidth in [("internal", "--", 1.2), ("internal-external", "-", 1.8)]:
            ld2_by_key = {row_key(row): row for row in grouped[("LD2", mode)]}
            ratios = []
            for row in grouped[(target, mode)]:
                ld2 = ld2_by_key.get(row_key(row))
                if ld2 is None:
                    continue
                item = dict(row)
                item["ratio"] = float(ld2["rad_to_born"]) / float(row["rad_to_born"])
                ratios.append(item)
            for idx, label in enumerate(Q2_LABELS):
                subset = rows_for_bin(ratios, label)
                if not subset:
                    continue
                ax.plot(
                    [row["xB"] for row in subset],
                    [row["ratio"] for row in subset],
                    color=COLORS[idx % len(COLORS)],
                    marker=MARKERS[idx % len(MARKERS)],
                    linestyle=style,
                    linewidth=linewidth,
                    markersize=5,
                    label=f"{mode} {label}",
                )
        ax.axhline(1.0, color="gray", linestyle=":", linewidth=1)
        ax.set_title(f"LD2 / {target}")
        ax.set_xlabel(r"$x_B$")
        ax.set_ylabel(r"$RC_{LD2}/RC_A$")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    fig.suptitle("RC factor entering nuclear ratios")
    fig.tight_layout()
    path = out_dir / "rgd_nuclear_ratio_rc_factor.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_net_external_on_ratio(rows: list[dict[str, object]], out_dir: Path) -> Path:
    grouped = by_target_mode(rows)
    fig, axes = setup_axes(2, 2, (15, 10))
    targets = ["C1", "C2", "Cu", "Sn"]
    ld2_int = {row_key(row): row for row in grouped[("LD2", "internal")]}
    ld2_ext = {row_key(row): row for row in grouped[("LD2", "internal-external")]}
    for ax, target in zip(axes, targets):
        target_ext = {row_key(row): row for row in grouped[(target, "internal-external")]}
        values = []
        for target_int in grouped[(target, "internal")]:
            key = row_key(target_int)
            if key not in ld2_int or key not in ld2_ext or key not in target_ext:
                continue
            ratio_int = float(ld2_int[key]["rad_to_born"]) / float(target_int["rad_to_born"])
            ratio_ext = float(ld2_ext[key]["rad_to_born"]) / float(target_ext[key]["rad_to_born"])
            item = dict(target_int)
            item["net_external_effect"] = ratio_ext / ratio_int
            values.append(item)
        for idx, label in enumerate(Q2_LABELS):
            subset = rows_for_bin(values, label)
            if not subset:
                continue
            ax.plot(
                [row["xB"] for row in subset],
                [row["net_external_effect"] for row in subset],
                color=COLORS[idx % len(COLORS)],
                marker=MARKERS[idx % len(MARKERS)],
                linestyle="-",
                linewidth=1.8,
                markersize=5,
                label=label,
            )
        ax.axhline(1.0, color="gray", linestyle=":", linewidth=1)
        ax.set_title(f"LD2 / {target}")
        ax.set_xlabel(r"$x_B$")
        ax.set_ylabel("net external effect")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, title=r"$Q^2$")
    fig.suptitle("Net external effect on nuclear-ratio RC")
    fig.tight_layout()
    path = out_dir / "rgd_net_external_effect_on_nuclear_ratio.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="results/tables/dis_rc_reference_grid.csv")
    parser.add_argument("--output-dir", default="results/plots")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    input_path = root / args.input
    out_dir = root / args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = read_rows(input_path)
    modes = sorted({str(row["radiation_mode"]) for row in rows})
    targets = sorted({str(row["target"]) for row in rows})
    if modes != ["internal", "internal-external"]:
        raise RuntimeError(f"expected both radiation modes, found {modes}")
    if "LD2" not in targets:
        raise RuntimeError("LD2 rows are required for nuclear-ratio plots")
    if any(not math.isfinite(float(row["rad_to_born"])) for row in rows):
        raise RuntimeError("non-finite correction factor in input table")

    outputs = [
        plot_mode_comparison(rows, out_dir),
        plot_external_effect(rows, out_dir),
        plot_nuclear_ratio(rows, out_dir),
        plot_net_external_on_ratio(rows, out_dir),
    ]
    print(f"Read {len(rows)} rows from {input_path}")
    for path in outputs:
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
