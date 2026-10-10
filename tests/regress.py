#!/usr/bin/env python3
"""Compare two Flux_e.dat-style files (col 0 nu, then apertures, last = total).

usage: regress.py ref.dat new.dat [--nu-min 1e9 --nu-max 1e12] [--tol 0.05]
Prints the max relative difference per column over the band (rows where the
reference is > 1e-3 of its column maximum, to skip numerically empty tails).
Exit status 1 if any column exceeds --tol.
"""
import argparse, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('ref'); ap.add_argument('new')
ap.add_argument('--nu-min', type=float, default=1e9)
ap.add_argument('--nu-max', type=float, default=1e12)
ap.add_argument('--tol', type=float, default=0.05)
a = ap.parse_args()

def load(f):  # skip '#' lines and the unmarked text header of legacy files
    rows = []
    for l in open(f):
        try: rows.append([float(x) for x in l.split()])
        except ValueError: continue
    return np.array(rows)
r = load(a.ref); n = load(a.new)
if r.shape != n.shape:
    sys.exit(f'shape mismatch {r.shape} vs {n.shape}')
m = (r[:, 0] >= a.nu_min) & (r[:, 0] <= a.nu_max)
bad = False
for k in range(1, r.shape[1]):
    ref = r[m, k]; new = n[m, k]
    sel = ref > 1e-3 * ref.max() if ref.max() > 0 else np.zeros_like(ref, bool)
    if not sel.any():
        print(f'col {k}: reference empty in band'); continue
    d = np.max(np.abs(new[sel] / ref[sel] - 1))
    print(f'col {k}: max rel diff {d:.3e}  ({sel.sum()} pts)')
    bad |= d > a.tol
sys.exit(1 if bad else 0)
