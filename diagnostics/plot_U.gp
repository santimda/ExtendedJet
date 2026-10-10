reset
#set encoding iso_8859_1
#set size square
#set title 'Dynamic evolution of the cooling dependences'
#set grid
#set tics scale 2
#set key box opaque width 0.95 height 1.0 bottom left
set key samplen 3.5 spacing 1.3 font ',20'
set format x '10^{%L}'
set format y '10^{%L}'
set logscale
#set xrange [1:8]
set yrange [1e-8:]
#set mxtics 10
set mytics 10
set xlabel 'z_j'
set ylabel 'U [erg s^{-1}]'

f1(x) = 1.0/(x/3e15)
f2(x) = 1.0/(x/3e15)**2


plot '../output/U_FF.dat' u 2:4 t 'U_B' w l lw 3 lc 1,\
	'../output/U_FF.dat' u 2:3 t 'U_{sy}' w l lw 3 lc 2,\
	'../output/U_e.dat' u 2:3 t 'U_{e}' w l lw 3 lc 3,\
	f1(x) t 'z^{-1}' w l lw 2 dt 2 lc rgb "grey30",\
	f2(x) t 'z^{-2}' w l lw 2 dt 3 lc rgb "grey40"
	

set terminal postscript color enhanced dl 1.2 lw 2 'Helvetica' 20
set output 'U_FF.ps'
replot
!./fixbb U_FF.ps
