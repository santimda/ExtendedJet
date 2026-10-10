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

plot '../output/losses_50.dat' u 1:2 w l dt 1 lw 2 lc rgb 'red' t 'Synchrotron' ,\
'../output/losses_50.dat' u 1:3 w l dt 1 lw 2 lc rgb 'blue' t 'IC' ,\
'../output/losses_50.dat' u 1:4 w l dt 1 lw 2 lc rgb 'green' t 'Adiabatic' ,\
'../output/losses_50.dat' u 1:5 w l dt 1 lw 3 lc rgb 'black' t 'Total' ,\
'../output/losses_50.dat' u 1:6 w l dt 2 lw 2 lc rgb 'cyan' t 'Advection (cell)' ,\
'../output/losses_50.dat' u 1:7 w l dt 2 lw 2 lc rgb 'orange' t 'Acceleration' ,\
'../output/losses_50.dat' u 1:8 w l dt 2 lw 2 lc rgb 'purple' t 'Diffusion'

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output 'losses50.ps'
replot
!./fixbb losses50.ps


reset
set xlabel "E' [eV]"
set ylabel "t' [s]"
set logscale
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1e6:1e14]
#set yrange [1e0:1e5]
set key right

plot '../output/losses_100.dat' u 1:2 w l dt 1 lw 2 lc rgb 'red' t 'Synchrotron' ,\
'../output/losses_100.dat' u 1:3 w l dt 1 lw 2 lc rgb 'blue' t 'IC' ,\
'../output/losses_100.dat' u 1:4 w l dt 1 lw 2 lc rgb 'green' t 'Adiabatic' ,\
'../output/losses_100.dat' u 1:5 w l dt 1 lw 3 lc rgb 'black' t 'Total' ,\
'../output/losses_100.dat' u 1:6 w l dt 2 lw 2 lc rgb 'cyan' t 'Advection (cell)' ,\
'../output/losses_100.dat' u 1:7 w l dt 2 lw 2 lc rgb 'orange' t 'Acceleration' ,\
'../output/losses_100.dat' u 1:8 w l dt 2 lw 2 lc rgb 'purple' t 'Diffusion'

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output 'losses100.ps'
replot
!./fixbb losses100.ps


reset
set xlabel "E' [eV]"
set ylabel "t' [s]"
set logscale
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1e6:1e14]
#set yrange [1e0:1e5]
set key right

plot '../output/losses_200.dat' u 1:2 w l dt 1 lw 2 lc rgb 'red' t 'Synchrotron' ,\
'../output/losses_200.dat' u 1:3 w l dt 1 lw 2 lc rgb 'blue' t 'IC' ,\
'../output/losses_200.dat' u 1:4 w l dt 1 lw 2 lc rgb 'green' t 'Adiabatic' ,\
'../output/losses_200.dat' u 1:5 w l dt 1 lw 3 lc rgb 'black' t 'Total' ,\
'../output/losses_200.dat' u 1:6 w l dt 2 lw 2 lc rgb 'cyan' t 'Advection (cell)' ,\
'../output/losses_200.dat' u 1:7 w l dt 2 lw 2 lc rgb 'orange' t 'Acceleration' ,\
'../output/losses_200.dat' u 1:8 w l dt 2 lw 2 lc rgb 'purple' t 'Diffusion'

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output 'losses200.ps'
replot
!./fixbb losses200.ps
