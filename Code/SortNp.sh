#!/bin/bash

#Number of energy steps for N(E)
nNE=$1

#First sort by segment number and then introduce a blank line between each segment
sort -k3,3n -g -k1,1 'Np_cuts.dat' -o 'Np_cuts_int.dat'
awk -v n=$nNE 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'Np_cuts_int.dat' > 'Np_cuts_ord.dat'

rm 'Np_cuts_int.dat'

echo '   --> Np_cuts.dat sorted by segment number in Np_cuts_ord.dat'
