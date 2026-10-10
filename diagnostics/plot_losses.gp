reset
set xlabel "E' [eV]"
set ylabel "t' [s]"
set logscale
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1e6:1e14]
#set yrange [1e0:1e5]
set key right

plot '../output/losses_***.dat' u 1:2 w l dt 1 lw 2 lc rgb 'red' t 'Synchrotron' ,\
'../output/losses_***.dat' u 1:3 w l dt 1 lw 2 lc rgb 'blue' t 'IC' ,\
'../output/losses_***.dat' u 1:4 w l dt 1 lw 2 lc rgb 'green' t 'Adiabatic' ,\
'../output/losses_***.dat' u 1:5 w l dt 1 lw 3 lc rgb 'black' t 'Total' ,\
'../output/losses_***.dat' u 1:6 w l dt 2 lw 2 lc rgb 'cyan' t 'Advection (cell)' ,\
'../output/losses_***.dat' u 1:7 w l dt 2 lw 2 lc rgb 'orange' t 'Acceleration' ,\
'../output/losses_***.dat' u 1:8 w l dt 2 lw 2 lc rgb 'purple' t 'Diffusion'

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output 'losses***.ps'
replot
!./fixbb losses***.ps

reset
set xlabel "E' [eV]"
set ylabel "t' [s]"
set logscale
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1e6:1e14]
#set yrange [1e0:1e5]
set key right

plot '../output/losses_1.dat' u 1:2 w l dt 1 lw 2 lc rgb 'red' t 'Synchrotron' ,\
'../output/losses_1.dat' u 1:3 w l dt 1 lw 2 lc rgb 'blue' t 'IC' ,\
'../output/losses_1.dat' u 1:4 w l dt 1 lw 2 lc rgb 'green' t 'Adiabatic' ,\
'../output/losses_1.dat' u 1:5 w l dt 1 lw 3 lc rgb 'black' t 'Total' ,\
'../output/losses_1.dat' u 1:6 w l dt 2 lw 2 lc rgb 'cyan' t 'Advection (cell)' ,\
'../output/losses_1.dat' u 1:7 w l dt 2 lw 2 lc rgb 'orange' t 'Acceleration' ,\
'../output/losses_1.dat' u 1:8 w l dt 2 lw 2 lc rgb 'purple' t 'Diffusion'

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output 'losses1.ps'
replot
!./fixbb losses1.ps


reset
set xlabel "E' [eV]"
set ylabel "t' [s]"
set logscale
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1e6:1e14]
#set yrange [1e0:1e5]
set key right

plot '../output/losses_271.dat' u 1:2 w l dt 1 lw 2 lc rgb 'red' t 'Synchrotron' ,\
'../output/losses_271.dat' u 1:3 w l dt 1 lw 2 lc rgb 'blue' t 'IC' ,\
'../output/losses_271.dat' u 1:4 w l dt 1 lw 2 lc rgb 'green' t 'Adiabatic' ,\
'../output/losses_271.dat' u 1:5 w l dt 1 lw 3 lc rgb 'black' t 'Total' ,\
'../output/losses_271.dat' u 1:6 w l dt 2 lw 2 lc rgb 'cyan' t 'Advection (cell)' ,\
'../output/losses_271.dat' u 1:7 w l dt 2 lw 2 lc rgb 'orange' t 'Acceleration' ,\
'../output/losses_271.dat' u 1:8 w l dt 2 lw 2 lc rgb 'purple' t 'Diffusion'

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output 'losses271.ps'
replot
!./fixbb losses271.ps


reset
set xlabel "E' [eV]"
set ylabel "t' [s]"
set logscale
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1e6:1e14]
#set yrange [1e0:1e5]
set key right

plot '../output/losses_543.dat' u 1:2 w l dt 1 lw 2 lc rgb 'red' t 'Synchrotron' ,\
'../output/losses_543.dat' u 1:3 w l dt 1 lw 2 lc rgb 'blue' t 'IC' ,\
'../output/losses_543.dat' u 1:4 w l dt 1 lw 2 lc rgb 'green' t 'Adiabatic' ,\
'../output/losses_543.dat' u 1:5 w l dt 1 lw 3 lc rgb 'black' t 'Total' ,\
'../output/losses_543.dat' u 1:6 w l dt 2 lw 2 lc rgb 'cyan' t 'Advection (cell)' ,\
'../output/losses_543.dat' u 1:7 w l dt 2 lw 2 lc rgb 'orange' t 'Acceleration' ,\
'../output/losses_543.dat' u 1:8 w l dt 2 lw 2 lc rgb 'purple' t 'Diffusion'

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output 'losses543.ps'
replot
!./fixbb losses543.ps
