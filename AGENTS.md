# AGENTS.md: ExtendedJet (branch m87-paper)

1D multizone jet code (Fortran 90 + OpenMP, optional f2py module). Used for the M87 RIAF+jet paper.

## Build and run
- `cd Code && make` builds `run.x` (`make depend` regenerates `Makefile.depend`; `make f2py` builds the Python module `JetCode`).
- `OMP_NUM_THREADS=1 ./run.x params.nml`. One full model takes about 40 s single-threaded (more threads are slower).
- Namelist groups `&model` (physics; lengths in R_g, energies in eV) and `&run` (grid, apertures, switches, `outdir`,
  `sync_approx='full'|'delta'`, `fix_norm_taper`, `fix_inj_emin`). Defaults live in `Code/system_parameters.f90`.
  Local configs: `configs/m87_reference.nml`, `configs/test.nml`.
- SSC is compiled off (`controls.f90`); do not enable.

## Outputs (in `outdir`)
- `Flux_e.dat`, `FluxCJ_e.dat`: col 0 = nu [Hz]; then the flux integrated up to each aperture R* (a zone counts when
  z sin i <= R*; apertures listed in the 2nd header line); last col = total. Jet and counter-jet in separate files.
- `Ne_cuts.dat` (E, N(E), segment), `quant_FF*.dat`, `trajectory*.dat`, `U_e.dat`, `SEDe_syn_cuts*.dat`.

## Tests
`tests/regress.py ref.dat new.dat` prints the max relative difference per column in the radio band.
`make -C Code test` runs the numerical unit tests and the reduced `run.x` regression against
the frozen references in `tests/reference/test/`.
