# AGENTS.md: ExtendedJet (branch m87-paper)

1D multizone jet code (Fortran 90 + OpenMP, optional f2py module). Used for the M87 RIAF+jet paper.

## Build and run
- `cd Code && make` builds `run.x` (`make depend` regenerates `Makefile.depend`; `make f2py` builds the Python module `JetCode`).
- Set `OMP_NUM_THREADS` explicitly when benchmarking or calling through f2py. On the current machine,
  `configs/m87_reference.nml` takes roughly 1.4-1.6 s after hoisting the injection normalization out of the
  energy loop; 2-4 threads do not currently improve on one thread. These timings are machine-specific; use
  this full configuration rather than `configs/test.nml` for profiling.
- Namelist groups `&model` (physics; lengths in R_g, energies in eV) and `&run` (grid, apertures, switches, `outdir`,
  `sync_approx='full'|'delta'`, `fix_norm_taper`, `fix_inj_emin`). Defaults live in `Code/system_parameters.f90`.
  Local configs: `configs/m87_reference.nml`, `configs/test.nml`.
- SSC is compiled off (`controls.f90`); do not enable.
- `SEDssc_data` is intentionally retained as an incomplete future SSC routine; its unused-argument warnings are expected.
- The Fortran modules hold mutable global state, so concurrent f2py calls in one process are not assumed safe.

## Outputs (in `outdir`)
- `Flux_e.dat`, `FluxCJ_e.dat`: col 0 = nu [Hz]; then the flux integrated up to each aperture R* (a zone counts when
  z sin i <= R*; apertures listed in the 2nd header line); last col = total. Jet and counter-jet in separate files.
- `Ne_cuts.dat` (E, N(E), segment), `quant_FF*.dat`, `trajectory*.dat`, `U_e.dat`, `SEDe_syn_cuts*.dat`.

## Tests
`tests/regress.py ref.dat new.dat` prints the max relative difference per column in the radio band.
`make -C Code test` runs the numerical unit tests and the reduced `run.x` regression against
the frozen references in `tests/reference/test/`.
For runtime checks and compiler warnings, run `make -C Code clean && make -C Code DEBUG=Y test`.
For profiling symbols with production optimization, run `make -C Code clean && make -C Code PROFILE=Y`.
`tests/f2py_check.py` currently requires an external Jet+RIAF namelist whose location is pending confirmation.
