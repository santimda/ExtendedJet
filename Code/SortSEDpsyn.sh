#!/bin/bash

#Number of photon energy steps for the synchrotron SED
NSED_sy=$1

if [ -f "./SEDpsyn_cuts.dat" ]; then
	sort -k4,4n -k1,1g 'SEDpsyn_cuts.dat' -o 'SEDpsyn_cuts_int.dat'
	awk -v n=$NSED_sy 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'SEDpsyn_cuts_int.dat' > 'SEDpsyn_cuts_ord.dat'
	echo '   --> SEDpsyn_cuts.dat sorted by segment number in SEDpsyn_cuts_ord.dat'	
fi


if [ -f "./SEDpsyn_cutsCJ.dat" ]; then
	sort -k4,4n -k1,1g 'SEDpsyn_cutsCJ.dat' -o 'SEDpsyn_cuts_int.dat'
	awk -v n=$NSED_sy 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'SEDpsyn_cuts_int.dat' > 'SEDpsyn_cuts_ordCJ.dat'
	echo '   --> SEDpsyn_cutsCJ.dat sorted by segment number in SEDpsyn_cuts_ordCJ.dat'
fi

rm 'SEDpsyn_cuts_int.dat'
