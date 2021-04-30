reset
#set encoding iso_8859_1
#set size square
#set title 'Dynamic evolution of the cooling dependences'
#set grid
#set tics scale 2
#set key box opaque width 0.95 height 1.0 bottom left
set key samplen 3.5 spacing 1.3 font ',20'
#set format x '10^{%L}'
#set format y '10^{%L}'
#set logscale
set xrange [:8]
set yrange [-1:1]
set mxtics 10
set mytics 10
set xlabel 'l_j/a'
set ylabel ''

f(x) = 0

#x1 = 4.e5
#set arrow from x1, graph 0.007 to x1, graph 0.993 nohead lc rgb 'grey88' lw 140

#plot './Quant_FF_norm.dat' u 1:2 t 'U_B/U_{B,0}' w l lw 3 lc 2,\
#	'./Quant_FF_norm.dat' u 1:3 t 'U_{ph}/U_{ph,0}' w l lw 3 dt 2 lc 4,\
#	'./Quant_FF_norm.dat' u 1:4 t '(v_{/Symbol \136}/r_{j})/(v_{{/Symbol \136},0}/r_{j,0})' w l lw 3 dt 3
	
plot './Quant_FF_norm.dat' u 1:5 t 'cos(wo) jet' w l lw 3 lc 2,\
	'./Quant_FF_norm.dat' u 1:6 t 'cos(woFF) jet' w l lw 3 lc 2 dt 2,\
	'./Quant_FF_normCJ.dat' u 1:5 t 'cos(wo) cjet' w l lw 3 lc 1,\
	'./Quant_FF_normCJ.dat' u 1:6 t 'cos(woFF) cjet' w l lw 3 lc 1 dt 2,\
	f(x) lw 2 lc 4 notitle
	
#,\	'./Quant_FF_norm.dat' u 1:($8/inye) w l lw 3 dt 4 lc 7 t 'L^{,}_{inj}/L_j'
#	'./Quant_FF_norm.dat' u 1:4 t 'r_{j}/a' w l lw 3 dt 3,\
#	'./Quant_FF_norm.dat' u 1:5 t 'v_{perp}/c' w l lc 7 lw 3 dt 3,\

set terminal postscript color enhanced dl 1.2 lw 2 'Helvetica' 20
set output 'Quant_FF.ps'
replot
!./fixbb Quant_FF.ps

