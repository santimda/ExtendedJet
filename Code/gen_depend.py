#!/usr/bin/env python3
"""Write Makefile.depend from the `use` statements of the *.f90 files (run in Code/)."""
import re, glob
mods = {}
for f in glob.glob('*.f90'):
    for l in open(f):
        m = re.match(r'\s*module\s+(\w+)\s*$', l, re.I)
        if m: mods[m.group(1).lower()] = f
out = []
for f in sorted(glob.glob('*.f90')):
    deps = set()
    for l in open(f):
        m = re.match(r'\s*use\s+(\w+)', l, re.I)
        if m and m.group(1).lower() in mods and mods[m.group(1).lower()] != f:
            deps.add(mods[m.group(1).lower()][:-4] + '.o')
    if deps: out.append(f"{f[:-4]}.o: {' '.join(sorted(deps))}")
open('Makefile.depend', 'w').write("\n".join(out) + "\n")
