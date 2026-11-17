"""Likelihood for aperture SED measurements."""

from __future__ import annotations

import time

import bilby
import numpy as np

from .config import FluxData

ERROR_FRACTION = 0.15


class SEDLikelihood(bilby.Likelihood):
    def __init__(self, model, data: FluxData, fit_parameters, fixed_parameters=None):
        super().__init__()
        self.model = model
        self.data = data
        self.frequencies = data.frequencies_ghz
        self.apertures = data.apertures_rg
        self.fit_parameters = tuple(fit_parameters)
        self.fixed_parameters = dict(fixed_parameters or {})
        self.observed = data.observed
        if data.errors_cgs is None:
            floor = 1e-6 * np.nanmax(data.fluxes_cgs)
            self.sigma = np.where(
                self.observed,
                np.maximum(ERROR_FRACTION * data.fluxes_cgs, floor),
                np.nan,
            )
        else:
            self.sigma = data.errors_cgs
        self._evaluation = 0
        self.evaluation_seconds = 0.0

    def evaluate_model(self, parameters, frequencies=None, apertures=None, output_dir=None):
        model_parameters = self.fixed_parameters.copy()
        model_parameters.update(parameters)
        return self.model(
            self.frequencies if frequencies is None else frequencies,
            self.apertures if apertures is None else apertures,
            model_parameters,
            output_dir=output_dir,
        )

    def log_likelihood(self, parameters):
        self._evaluation += 1
        started = time.perf_counter()
        try:
            prediction = self.evaluate_model(parameters)
        except (OSError, RuntimeError, ValueError):
            return -1e30
        finally:
            self.evaluation_seconds += time.perf_counter() - started
        residual = (self.data.fluxes_cgs - prediction)[self.observed] / self.sigma[self.observed]
        sigma = self.sigma[self.observed]
        return float(-0.5 * np.sum(residual**2 + np.log(2.0 * np.pi * sigma**2)))
