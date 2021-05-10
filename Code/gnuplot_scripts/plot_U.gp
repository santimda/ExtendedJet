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
#set xrange [:8]
#set yrange [-1:1]
#set mxtics 10
set mytics 10
set xlabel 'z_j'
set ylabel 'U [erg s^{-1}]'

f(x) = 0

#x1 = 4.e5
#set arrow from x1, graph 0.007 to x1, graph 0.993 nohead lc rgb 'grey88' lw 140

plot '../U_FF.dat' u 2:4 t 'U_B' w l lw 3 lc 1,\
	'../U_FF.dat' u 2:3 t 'U_{sy}' w l lw 3 lc 2,\
	'../U_e.dat' u 2:3 t 'U_{e}' w l lw 3 lc 3    
	
#plot './Quant_FF_norm.dat' u 1:5 t 'cos(wo) jet' w l lw 3 lc 2,\
#	'./Quant_FF_norm.dat' u 1:6 t 'cos(woFF) jet' w l lw 3 lc 2 dt 2,\
#	'./Quant_FF_normCJ.dat' u 1:5 t 'cos(wo) cjet' w l lw 3 lc 1,\
#	'./Quant_FF_normCJ.dat' u 1:6 t 'cos(woFF) cjet' w l lw 3 lc 1 dt 2,\
#	f(x) lw 2 lc 4 notitle
	

set terminal postscript color enhanced dl 1.2 lw 2 'Helvetica' 20
set output 'U_FF.ps'
replot
!./fixbb U_FF.ps

