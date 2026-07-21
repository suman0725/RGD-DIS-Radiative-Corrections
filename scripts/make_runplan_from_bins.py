#!/usr/bin/env python3
"""Convert an xB-Q2 bin definition CSV into an EXTERNALS run plan."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


NUCLEON_MASS_GEV = 0.9382720813


def center(row: dict[str, str], name: str) -> float:
    direct = row.get(f"{name}_center") or row.get(f"{name}_mean")
    if direct not in (None, ""):
        return float(direct)
    return 0.5 * (float(row[f"{name}_min"]) + float(row[f"{name}_max"]))


def kinematics(ebeam: float, xb: float, q2: float) -> tuple[float, float]:
    nu = q2 / (2.0 * NUCLEON_MASS_GEV * xb)
    eprime = ebeam - nu
    if eprime <= 0.0:
        raise ValueError(f"invalid eprime={eprime:.6g} for xB={xb}, Q2={q2}")
    sin2_half_theta = q2 / (4.0 * ebeam * eprime)
    if not 0.0 < sin2_half_theta < 1.0:
        raise ValueError(
            f"invalid sin^2(theta/2)={sin2_half_theta:.6g} for xB={xb}, Q2={q2}"
        )
    theta = math.degrees(2.0 * math.asin(math.sqrt(sin2_half_theta)))
    return eprime, theta


def read_bins(path: Path, ebeam: float, version: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"bin_id", "xB_min", "xB_max", "Q2_min", "Q2_max"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(f"{path} is missing columns: {', '.join(sorted(missing))}")

        for source in reader:
            xb = center(source, "xB")
            q2 = center(source, "Q2")
            eprime, theta = kinematics(ebeam, xb, q2)
            rows.append(
                {
                    "binning_version": version,
                    "bin_id": source["bin_id"],
                    "xB_min": float(source["xB_min"]),
                    "xB_max": float(source["xB_max"]),
                    "Q2_min": float(source["Q2_min"]),
                    "Q2_max": float(source["Q2_max"]),
                    "xB_center": xb,
                    "Q2_center": q2,
                    "Ebeam": ebeam,
                    "Eprime": eprime,
                    "theta": theta,
                }
            )
    if not rows:
        raise RuntimeError(f"no bins found in {path}")
    return rows


def write_run_plan(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write("RUNPLAN for RG-D xB-Q2 bin centers\n")
        stream.write(" inclusive:\n\n\n\n")
        stream.write("   E     Ep    theta     W     y     x      Q2\n")
        for row in rows:
            stream.write(
                f"{float(row['Ebeam']):5.2f}  "
                f"{float(row['Eprime']):6.4f}  "
                f"{float(row['theta']):7.4f}\n"
            )


def write_bin_map(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "binning_version",
        "bin_id",
        "xB_min",
        "xB_max",
        "Q2_min",
        "Q2_max",
        "xB_center",
        "Q2_center",
        "Ebeam",
        "Eprime",
        "theta",
    ]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bins_csv", type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--ebeam", type=float, default=10.5)
    parser.add_argument("--run-plan", required=True, type=Path)
    parser.add_argument("--bin-map", required=True, type=Path)
    args = parser.parse_args()

    rows = read_bins(args.bins_csv, args.ebeam, args.version)
    write_run_plan(args.run_plan, rows)
    write_bin_map(args.bin_map, rows)
    print(f"Read {len(rows)} bins from {args.bins_csv}")
    print(f"Wrote run plan: {args.run_plan}")
    print(f"Wrote bin map: {args.bin_map}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
