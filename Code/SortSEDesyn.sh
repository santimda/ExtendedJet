#!/bin/bash

#Number of photon energy steps for the synchrotron SED
NSED_sy=$1

if [ -f "./SEDe_syn_cuts.dat" ]; then
	sort -k4,4n -k1,1g 'SEDe_syn_cuts.dat' -o 'SEDe_syn_cuts_int.dat'
	awk -v n=$NSED_sy 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'SEDe_syn_cuts_int.dat' > 'SEDe_syn_cuts_ord.dat'
	echo '   --> SEDe_syn_cuts.dat sorted by segment number in SEDe_syn_cuts_ord.dat'	
fi


if [ -f "./SEDe_syn_cutsCJ.dat" ]; then
	sort -k4,4n -k1,1g 'SEDe_syn_cutsCJ.dat' -o 'SEDe_syn_cuts_int.dat'
	awk -v n=$NSED_sy 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'SEDe_syn_cuts_int.dat' > 'SEDe_syn_cuts_ordCJ.dat'
	echo '   --> SEDe_syn_cutsCJ.dat sorted by segment number in SEDe_syn_cuts_ordCJ.dat'
fi

rm 'SEDe_syn_cuts_int.dat'