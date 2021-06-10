#!/bin/bash

#Number of energy steps for N(E)
nNE=$1

#First sort by segment number and then introduce a blank line between each segment
sort -k3,3n -g -k1,1 'NE_cuts.dat' -o 'NE_cuts_int.dat'
awk -v n=$nNE 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'NE_cuts_int.dat' > 'NE_cuts_ord.dat'

rm 'NE_cuts_int.dat'

echo '   --> NE_cuts.dat sorted by segment number in NE_cuts_ord.dat'
