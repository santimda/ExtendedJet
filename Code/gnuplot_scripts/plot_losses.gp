reset
set xlabel "E' [eV]"
set ylabel "t' [s]"
set logscale
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1e6:1e14]
#set yrange [1e0:1e5]
set key right

plot '../losses.dat' u 1:2 w l dt 1 lw 2 lc rgb 'red' t 'Synchrotron' ,\
'../losses.dat' u 1:3 w l dt 1 lw 2 lc rgb 'blue' t 'IC' ,\
'../losses.dat' u 1:4 w l dt 1 lw 2 lc rgb 'green' t 'Adiabatic' ,\
'../losses.dat' u 1:5 w l dt 1 lw 3 lc rgb 'black' t 'Total' ,\
'../losses.dat' u 1:6 w l dt 2 lw 2 lc rgb 'cyan' t 'Advection (cell)' ,\
'../losses.dat' u 1:7 w l dt 2 lw 2 lc rgb 'orange' t 'Acceleration' ,\
'../losses.dat' u 1:8 w l dt 2 lw 2 lc rgb 'purple' t 'Diffusion'

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output 'losses.ps'
replot
!./fixbb losses.ps

