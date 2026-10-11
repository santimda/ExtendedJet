#!/usr/bin/env python3
"""Run the reduced Fortran model and compare its fluxes with references."""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CODE_DIR = ROOT / "Code"
CONFIG = ROOT / "configs" / "test.nml"
EXECUTABLE = CODE_DIR / "run.x"
REGRESS = HERE / "regress.py"
REFERENCE_DIR = HERE / "reference" / "test"


def make_scratch_config(scratch_dir: Path) -> Path:
    """Copy the reduced namelist and redirect all generated output."""
    output_dir = scratch_dir / "output"
    output_dir.mkdir()
    config = scratch_dir / "test.nml"
    replaced = False
    lines = []
    for line in CONFIG.read_text().splitlines():
        if "outdir" in line and "=" in line and not replaced:
            lines.append(f"  outdir = '{output_dir}'")
            replaced = True
        else:
            lines.append(line)
    if not replaced:
        raise RuntimeError(f"Could not find outdir in {CONFIG}")
    config.write_text("\n".join(lines) + "\n")
    return config


def compare(reference: Path, result: Path, tolerance: float) -> None:
    subprocess.run(
        [
            "python3",
            str(REGRESS),
            str(reference),
            str(result),
            "--tol",
            str(tolerance),
        ],
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tol", type=float, default=0.05)
    args = parser.parse_args()

    if not EXECUTABLE.is_file():
        raise SystemExit(f"Missing executable: {EXECUTABLE}; run 'make -C Code' first")

    with tempfile.TemporaryDirectory(prefix="extendedjet_regression_") as tmp:
        scratch_dir = Path(tmp)
        config = make_scratch_config(scratch_dir)
        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = "1"
        subprocess.run([str(EXECUTABLE), str(config)], cwd=CODE_DIR, env=env, check=True)

        for name in ("Flux_e.dat", "FluxCJ_e.dat"):
            result = scratch_dir / "output" / name
            reference = REFERENCE_DIR / name
            if not result.is_file():
                raise SystemExit(f"Fortran run did not produce {result}")
            if not reference.is_file():
                raise SystemExit(f"Missing regression reference: {reference}")
            print(f"Comparing {name}")
            compare(reference, result, args.tol)


if __name__ == "__main__":
    main()
