#!/bin/bash

#Number of photon energy steps for the synchrotron SED
NSED_sy=$1

if [ -f "./SEDe_syn_cuts.dat" ]; then
    # Sort by the third column (segment index), then by the first (E)
    sort -k3,3n -k1,1g SEDe_syn_cuts.dat -o SEDe_syn_cuts_int.dat

    # Split into blocks every NSED_sy lines per segment group
    awk -v n=$NSED_sy 'NR > 1 && (NR-1) % n == 0 {print ""} {print}' SEDe_syn_cuts_int.dat > SEDe_syn_cuts_ord.dat

    echo "   --> SEDe_syn_cuts.dat sorted and grouped in SEDe_syn_cuts_ord.dat"
fi


if [ -f "./SEDe_syn_cutsCJ.dat" ]; then
    # Sort by the third column (segment index), then by the first (E)
	sort -k3,3n -k1,1g 'SEDe_syn_cutsCJ.dat' -o 'SEDe_syn_cuts_int.dat'

    # Split into blocks every NSED_sy lines per segment group
	awk -v n=$NSED_sy 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'SEDe_syn_cuts_int.dat' > 'SEDe_syn_cuts_ordCJ.dat'
	
	echo '   --> SEDe_syn_cutsCJ.dat sorted by segment number in SEDe_syn_cuts_ordCJ.dat'
fi

rm 'SEDe_syn_cuts_int.dat'