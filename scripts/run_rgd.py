#!/usr/bin/env python3
"""Run RG-D EXTERNALS targets and extract validated cross-section tables."""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
from pathlib import Path


TARGETS = {
    "LD2": "rgdD2_rgd",
    "C1": "rgdC1_rgd",
    "C2": "rgdC2_rgd",
    "Cu": "rgdCu_rgd",
    "Sn": "rgdSn_rgd",
}

RADIATION_MODES = {
    "internal": "0",
    "internal-external": "1",
}

COLUMNS = [
    "Ebeam", "Eprime", "theta", "xB", "Q2",
    "sigma_born", "sigma_born_in", "sigma_born_qe",
    "sigma_rad", "sigma_rad_el", "sigma_rad_qe", "sigma_rad_dis",
    "coulomb_factor",
]


def read_bin_map(path: Path | None, expected: int) -> list[dict[str, object]] | None:
    if path is None:
        return None
    rows: list[dict[str, object]] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"binning_version", "bin_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for source in reader:
            rows.append({key: value for key, value in source.items() if key not in {"Ebeam", "Eprime", "theta"}})
    if len(rows) != expected:
        raise RuntimeError(f"{path}: found {len(rows)} bin rows; expected {expected}")
    return rows


def physics_rows(stdout: str) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in stdout.splitlines():
        fields = line.split()
        if len(fields) != len(COLUMNS):
            continue
        try:
            values = [float(value) for value in fields]
        except ValueError:
            continue
        if values[0] <= 0 or values[1] <= 0 or values[1] >= values[0]:
            continue
        rows.append(values)
    return rows


def runtime_input(root: Path, stem: str, run_plan: Path, radiation_mode: str) -> Path:
    source = root / "INP" / f"{stem}.inp"
    lines = source.read_text().splitlines()
    if len(lines) < 4:
        raise RuntimeError(f"unexpected target input format: {source}")
    lines[1] = run_plan.relative_to(root).as_posix()

    path = root / "INP" / f"{stem}_{radiation_mode.replace('-', '_')}.inp"
    path.write_text("\n".join(lines) + "\n")
    return path.relative_to(root)


def run_target(
    root: Path,
    binary: Path,
    target: str,
    stem: str,
    expected: int,
    radiation_mode: str,
    bin_map: list[dict[str, object]] | None,
    run_plan: Path,
) -> list[dict[str, object]]:
    input_path = root / "INP" / f"{stem}.inp"
    if not input_path.exists():
        raise RuntimeError(f"missing target input: {input_path}")
    run_input = runtime_input(root, stem, run_plan, radiation_mode)

    env = os.environ.copy()
    env["RGD_DOEXT"] = RADIATION_MODES[radiation_mode]
    try:
        completed = subprocess.run(
            [str(binary)],
            cwd=root,
            input=f"{run_input}\n",
            text=True,
            capture_output=True,
            check=False,
            env=env,
        )
    finally:
        (root / run_input).unlink(missing_ok=True)
    log_dir = root / "results" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_stem = f"{stem}_{radiation_mode.replace('-', '_')}"
    (log_dir / f"{log_stem}.stdout.log").write_text(completed.stdout)
    (log_dir / f"{log_stem}.stderr.log").write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError(
            f"{target}: externals_all exited with {completed.returncode}; "
            f"see {log_dir / f'{log_stem}.stdout.log'} and {log_dir / f'{log_stem}.stderr.log'}"
        )

    parsed = physics_rows(completed.stdout)
    if len(parsed) != expected:
        raise RuntimeError(f"{target}: extracted {len(parsed)} physics rows; expected {expected}")

    rows: list[dict[str, object]] = []
    for index, values in enumerate(parsed):
        row: dict[str, object] = {
            "target": target,
            "radiation_mode": radiation_mode,
            **dict(zip(COLUMNS, values)),
        }
        if bin_map is not None:
            row = {**bin_map[index], **row}
        born = float(row["sigma_born"])
        radiated = float(row["sigma_rad"])
        if not (math.isfinite(born) and math.isfinite(radiated) and born > 0 and radiated > 0):
            raise RuntimeError(f"{target}: invalid cross section at xB={row['xB']}, Q2={row['Q2']}")
        row["dis_rc_factor"] = born / radiated
        row["dis_rc_factor_err"] = 0.0
        rows.append(row)
    return rows


def run_plan_points(path: Path) -> int:
    count = 0
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            [float(value) for value in fields[:3]]
        except ValueError:
            continue
        count += 1
    if count == 0:
        raise RuntimeError(f"no kinematic points found in {path}")
    return count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-plan", default="RUNPLAN/rgd_kin.inp")
    parser.add_argument("--binary", default="./externals_all")
    parser.add_argument("--output", default="results/tables/dis_rc_reference_grid.csv")
    parser.add_argument(
        "--bin-map",
        default=None,
        help="optional CSV with binning_version and bin_id rows matching the run plan",
    )
    parser.add_argument(
        "--radiation",
        choices=["internal", "internal-external", "both"],
        default="both",
        help="radiation mode to run; default runs internal and internal-external",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    (root / "OUT").mkdir(exist_ok=True)
    run_plan = root / args.run_plan
    expected = run_plan_points(run_plan)
    bin_map = read_bin_map(root / args.bin_map if args.bin_map else None, expected)
    binary = (root / args.binary).resolve()
    if not binary.exists():
        raise SystemExit(f"missing executable: {binary}; run make first")

    rows: list[dict[str, object]] = []
    modes = list(RADIATION_MODES) if args.radiation == "both" else [args.radiation]
    for mode in modes:
        for target, stem in TARGETS.items():
            print(f"Running {target} ({stem}, {mode})...", flush=True)
            rows.extend(run_target(root, binary, target, stem, expected, mode, bin_map, run_plan))

    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["target", "radiation_mode", *COLUMNS, "dis_rc_factor", "dis_rc_factor_err"]
    if bin_map is not None:
        fields = [*bin_map[0].keys(), *fields]
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Validated targets: {', '.join(TARGETS)}")
    print(f"Radiation modes: {', '.join(modes)}")
    print(f"Kinematic points per target: {expected}")
    print(f"Total factor rows: {len(rows)}")
    print(f"Output: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
