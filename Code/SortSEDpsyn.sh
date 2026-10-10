#!/bin/bash

#Number of photon energy steps for the synchrotron SED
NSED_sy=$1

if [ -f "./SEDp_syn_cuts.dat" ]; then
	sort -k4,4n -k1,1g 'SEDp_syn_cuts.dat' -o 'SEDp_syn_cuts_int.dat'
	awk -v n=$NSED_sy 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'SEDp_syn_cuts_int.dat' > 'SEDp_syn_cuts_ord.dat'
	echo '   --> SEDp_syn_cuts.dat sorted by segment number in SEDp_syn_cuts_ord.dat'	
fi


if [ -f "./SEDp_syn_cutsCJ.dat" ]; then
	sort -k4,4n -k1,1g 'SEDp_syn_cutsCJ.dat' -o 'SEDp_syn_cuts_int.dat'
	awk -v n=$NSED_sy 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'SEDp_syn_cuts_int.dat' > 'SEDp_syn_cuts_ordCJ.dat'
	echo '   --> SEDp_syn_cutsCJ.dat sorted by segment number in SEDp_syn_cuts_ordCJ.dat'
fi

rm 'SEDp_syn_cuts_int.dat'
