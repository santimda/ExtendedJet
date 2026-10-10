#!/bin/bash

#Number of photon energy steps for the synchrotron SED
NSED_sy=$1

if [ -f "./SEDe_syn_cuts.dat" ]; then
    # Extract header and sort the rest
    { 
        head -n 1 SEDe_syn_cuts.dat
        tail -n +2 SEDe_syn_cuts.dat | sort -k3,3n -k1,1g
    } > SEDe_syn_cuts_int.dat

    # Split into blocks every NSED_sy lines per segment group (skip header)
    awk -v n=$NSED_sy 'NR == 1 {print; next} (NR > 2 && (NR-2) % n == 0) {print ""} {print}' SEDe_syn_cuts_int.dat > SEDe_syn_cuts_ord.dat

    echo "   --> SEDe_syn_cuts.dat sorted and grouped in SEDe_syn_cuts_ord.dat"
fi



if [ -f "./SEDe_syn_cutsCJ.dat" ]; then
    # Extract header and sort the rest
    { 
        head -n 1 SEDe_syn_cutsCJ.dat
        tail -n +2 SEDe_syn_cutsCJ.dat | sort -k3,3n -k1,1g
    } > SEDe_syn_cutsCJ_int.dat

    # Split into blocks every NSED_sy lines per segment group (skip header)
    awk -v n=$NSED_sy 'NR == 1 {print; next} (NR > 2 && (NR-2) % n == 0) {print ""} {print}' SEDe_syn_cutsCJ_int.dat > SEDe_syn_cutsCJ_ord.dat

    echo "   --> SEDe_syn_cutsCJ.dat sorted and grouped in SEDe_syn_cutsCJ_ord.dat"
fi


rm 'SEDe_syn_cuts_int.dat'