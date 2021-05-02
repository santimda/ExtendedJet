#!/bin/bash

#Number of photon energy steps for the synchrotron SED
NSED_sy=$1

sort -k4,4n -k1,1g 'SEDsyn_cuts.dat' -o 'SEDsyn_cuts_int.dat'
awk -v n=$NSED_sy 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'SEDsyn_cuts_int.dat' > 'SEDsyn_cuts_ord.dat'
#sort -k4,4n -k1,1g 'SEDsyn_cutsCJ.dat' -o 'SEDsyn_cuts_int.dat'
#awk -v n=$NSED_sy 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'SEDsyn_cuts_int.dat' > 'SEDsyn_cuts_ordCJ.dat'

rm 'SEDsyn_cuts_int.dat'
echo '   --> SEDsyn_cuts.dat sorted by segment number in SEDsyn_cuts_ord.dat'