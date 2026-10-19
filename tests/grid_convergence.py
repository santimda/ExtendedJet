#!/usr/bin/env python3
"""Run full-model grid convergence experiments.

This script varies the runtime ``Ns`` and ``numbins`` settings in isolated
copies of the M87 reference namelist. It compares every populated aperture
and total-flux column of both jet files over the full SED, omitting only rows
below 1e-3 of that column's reference maximum.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = ROOT / "Code"
CONFIG = ROOT / "configs" / "m87_reference.nml"
EXECUTABLE = CODE_DIR / "run.x"
FILES = ("Flux_e.dat", "FluxCJ_e.dat")


def make_config(directory: Path, parameter: str, value: int, code_dir: Path) -> Path:
    output = directory / "output"
    output.mkdir(parents=True)
    lines = []
    replaced = False
    for line in CONFIG.read_text().splitlines():
        if line.strip().startswith(f"{parameter} ="):
            lines.append(f"  {parameter} = {value}")
            replaced = True
        else:
            line = line.replace("outdir = '../output'", f"outdir = '{output}'")
            line = line.replace("scriptdir = '../Code'", f"scriptdir = '{code_dir}'")
            lines.append(line)
    if not replaced:
        raise RuntimeError(f"Could not find {parameter} in {CONFIG}")
    config = directory / f"{parameter}_{value}.nml"
    config.write_text("\n".join(lines) + "\n")
    return config


def load_flux(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text().splitlines():
        try:
            rows.append([float(value) for value in line.split()])
        except ValueError:
            continue
    return np.asarray(rows)


def compare(reference: dict[str, np.ndarray], directory: Path) -> tuple[float, float, dict[str, float]]:
    maximum_aperture = 0.0
    maximum_total = 0.0
    details = {}
    for filename in FILES:
        ref = reference[filename]
        new = load_flux(directory / "output" / filename)
        if ref.shape != new.shape or not np.allclose(ref[:, 0], new[:, 0]):
            raise RuntimeError(f"incompatible output grid in {filename}")
        for column in range(1, ref.shape[1]):
            column_ref = ref[:, column]
            threshold = 1e-3 * np.max(column_ref)
            selected = column_ref > threshold
            if np.any(selected):
                difference = np.max(np.abs(new[selected, column] / column_ref[selected] - 1.0))
                if column == ref.shape[1] - 1:
                    maximum_total = max(maximum_total, float(difference))
                else:
                    maximum_aperture = max(maximum_aperture, float(difference))
                details[f"{filename}:column_{column}"] = float(difference)
    return maximum_aperture, maximum_total, details


def run(parameter: str, values: list[int], reference_value: int, numdt: int | None) -> None:
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    with tempfile.TemporaryDirectory(prefix=f"extendedjet_{parameter}_") as temporary:
        root = Path(temporary)
        executable = EXECUTABLE
        code_dir = CODE_DIR
        if numdt is not None:
            code_dir = root / "Code"
            shutil.copytree(
                CODE_DIR,
                code_dir,
                ignore=shutil.ignore_patterns("build", "*.o", "*.mod", "run.x", "*.so", "JetCode*"),
            )
            globals_file = code_dir / "globals.f90"
            source = globals_file.read_text()
            source = source.replace("integer, parameter :: numdt = 100000", f"integer, parameter :: numdt = {numdt}")
            source = source.replace(
                "jlast = jz   ! last computed step (arrays are only filled up to here)",
                "jlast = jz   ! last computed step (arrays are only filled up to here)\n"
                "         if (jlast >= numdt) print *, 'WARNING: transport history reached numdt'",
            )
            globals_file.write_text(source)
            subprocess.run(["make"], cwd=code_dir, env=env, check=True, stdout=subprocess.DEVNULL)
            executable = code_dir / "run.x"

        reference_config = make_config(root / "reference", parameter, reference_value, code_dir)
        subprocess.run([str(executable), str(reference_config)], cwd=code_dir, env=env, check=True)
        reference = {name: load_flux(root / "reference" / "output" / name) for name in FILES}
        print(f"{parameter} reference={reference_value}")
        for value in values:
            case = root / str(value)
            case.mkdir()
            config = make_config(case, parameter, value, code_dir)
            start = time.perf_counter()
            subprocess.run([str(executable), str(config)], cwd=code_dir, env=env, check=True)
            elapsed = time.perf_counter() - start
            aperture_difference, total_difference, details = compare(reference, case)
            print(
                f"{parameter}={value:4d}  max_aperture={aperture_difference:.6e}"
                f"  max_total={total_difference:.6e}  seconds={elapsed:.2f}"
            )
            if parameter == "Ns":
                print("  " + "  ".join(f"{key}={value:.3e}" for key, value in details.items()))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", choices=("Ns", "numbins"), required=True)
    parser.add_argument("--values", nargs="+", type=int, required=True)
    parser.add_argument("--reference", type=int, required=True)
    parser.add_argument(
        "--numdt",
        type=int,
        help="Build an isolated executable with this transport-history limit.",
    )
    args = parser.parse_args()
    if not EXECUTABLE.is_file():
        raise SystemExit(f"Missing executable: {EXECUTABLE}; run 'make -C Code' first")
    run(args.grid, args.values, args.reference, args.numdt)


if __name__ == "__main__":
    main()
