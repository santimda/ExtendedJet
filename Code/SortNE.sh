#!/bin/bash

#Number of energy steps for N(E)
nNE=$1

#Sort by segment number and introduce blank line between each segment
sort -k3,3n -k1,1n 'NE_cuts.dat' -o 'NE_cuts_int.dat'
awk -v n=$nNE 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'NE_cuts_int.dat' > 'NE_cuts_ord.dat'

rm 'NE_cuts_int.dat'
