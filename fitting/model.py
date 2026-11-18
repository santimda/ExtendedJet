"""Exact-frequency ExtendedJet model interface."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "Code"


class ExtendedJetModel:
    """Evaluate circular-aperture synchrotron flux at requested frequencies."""

    log_parameters = {"Lj"}
    parameter_names = {
        "sigma_0", "zmax", "Lj", "a_inj", "zmax_inj", "gammaIndex",
        "geomIndex_out", "z0", "fnT", "facc", "p",
    }

    def __init__(self, namelist, output_dir, ns=None, numbins=None):
        sys.path.insert(0, str(CODE))
        import JetCode

        self.jetcode = JetCode
        self.system_parameters = JetCode.system_parameters
        self.system_parameters.read_namelist(str(Path(namelist).resolve()))
        if ns is not None:
            self.system_parameters.ns = int(ns)
        if numbins is not None:
            self.system_parameters.numbins = int(numbins)
        self.system_parameters.derive_parameters()
        # Fitting apertures are defined as projected circular apertures.
        JetCode.controls.aperture_area = True
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.source_name = Path(namelist).stem
        self.h = float(JetCode.constants.h)
        self.parameters = {
            "sigma_0": float(self.system_parameters.sigma_0),
            "zmax": float(self.system_parameters.zmax),
            "Lj": float(self.system_parameters.lj),
            "a_inj": float(self.system_parameters.a_inj),
            "zmax_inj": float(self.system_parameters.zmax_inj),
            "gammaIndex": float(self.system_parameters.gammaindex),
            "geomIndex_out": float(self.system_parameters.geomindex_out),
            "z0": float(self.system_parameters.z0),
            "fnT": float(self.system_parameters.fnt),
            "facc": float(self.system_parameters.facc),
            "p": float(self.system_parameters.p),
        }

    def __call__(self, frequencies_ghz, apertures_rg, parameters=None, output_dir=None):
        parameters = dict(parameters or {})
        log_parameter_names = {f"log_{name}" for name in self.log_parameters}
        unknown = set(parameters) - self.parameter_names - log_parameter_names
        if unknown:
            raise ValueError(f"Unknown model parameters: {sorted(unknown)}")
        for name in self.log_parameters:
            if name in parameters and f"log_{name}" in parameters:
                raise ValueError(f"provide either {name} or log_{name}, not both")
        values = self.parameters.copy()
        values.update(parameters)
        for name in self.log_parameters:
            log_name = f"log_{name}"
            if log_name in values:
                values[name] = 10.0 ** float(values.pop(log_name))
        frequencies_ghz = np.asarray(frequencies_ghz, dtype=float)
        apertures_rg = np.asarray(apertures_rg, dtype=float)
        if frequencies_ghz.ndim != 1 or np.any(frequencies_ghz <= 0):
            raise ValueError("frequencies_ghz must be one-dimensional and positive")
        if apertures_rg.ndim != 1 or np.any(apertures_rg <= 0):
            raise ValueError("apertures_rg must be one-dimensional and positive")
        if np.any(np.diff(frequencies_ghz) <= 0) or np.any(np.diff(apertures_rg) <= 0):
            raise ValueError("frequencies_ghz and apertures_rg must be strictly sorted")
        if len(frequencies_ghz) > 300 or len(apertures_rg) > 10:
            raise ValueError("the exact evaluator supports at most 300 frequencies and 10 apertures")
        target_dir = Path(output_dir) if output_dir is not None else self.output_dir
        target_dir.mkdir(parents=True, exist_ok=True)
        self.system_parameters.set_outdir(str(target_dir.resolve()))
        raw = self.jetcode.syn_model_apertures(
            values["sigma_0"], values["zmax"], values["Lj"], values["a_inj"],
            values["zmax_inj"], values["gammaIndex"], values["geomIndex_out"],
            values["z0"], values["fnT"], values["facc"], values["p"],
            self.h * frequencies_ghz * 1e9, apertures_rg,
        )
        result = np.asarray(raw)[:len(frequencies_ghz), :len(apertures_rg)]
        if result.shape != (len(frequencies_ghz), len(apertures_rg)) or not np.all(np.isfinite(result)):
            raise ValueError("exact evaluator returned an invalid flux array")
        return result
