#!/bin/bash

#Number of energy steps for N(E)
nNE=$1

#First sort by segment number and then introduce a blank line between each segment
sort -k3,3n -g -k1,1 'Ne_cuts.dat' -o 'Ne_cuts_int.dat'
awk -v n=$nNE 'FNR != 1 && (FNR-1) % n == 0{print ""}{print}' 'Ne_cuts_int.dat' > 'Ne_cuts_ord.dat'

rm 'Ne_cuts_int.dat'

echo '   --> Ne_cuts.dat sorted by segment number in Ne_cuts_ord.dat'
