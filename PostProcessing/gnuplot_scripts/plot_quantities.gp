reset
#set encoding iso_8859_1
#set size square
#set title 'Dynamic evolution of the cooling dependences'
#set grid
#set tics scale 2
set key box opaque bottom left
set key samplen 3.5 spacing 1.4 font ',20'
set format x '10^{%L}'
set format y '10^{%L}'
set logscale
#set xrange [:8]
set yrange [1e-2:2e1]
#set mxtics 10
set mytics 10
set xlabel 'z_j/R_g'
set ylabel ''

f(x) = 0

#x1 = 4.e5
#set arrow from x1, graph 0.007 to x1, graph 0.993 nohead lc rgb 'grey88' lw 140

plot '../../Code/quant_FF_norm.dat' u 1:2 t 'B/B_0' w l lw 3 lc 1,\
	'../../Code/quant_FF_norm.dat' u 1:3 t 'v_{/Symbol \136}/v_{{/Symbol \136},0}' w l lw 3 dt 3 lc 2,\
	'../../Code/quant_FF_norm.dat' u 1:4 t 'r/r_{0}' w l lw 3 dt 4 lc 3,\
	'../../Code/quant_FF_norm.dat' u 1:6 t '{/Symbol G}/{/Symbol G}_{0}' w l lw 3 dt 5 lc 5,\
	'../../Code/quant_FF_norm.dat' u 1:7 t '{/Symbol d}/{/Symbol d}_{0}' w l lw 3 dt 6 lc 6


	#'../../Code/quant_FF_norm.dat' u 1:5 t 'V/V_{0}' w l lw 3 dt 2 lc 4,\
	#'../../Code/quant_FF_norm.dat' u 1:($2**2) t 'U_B/U_{B,0}' w l lw 3 lc 7,\
	#'../../Code/quant_FF_norm.dat' u 1:($3/$4) t '(v_{/Symbol \136}/r_{j})/(v_{{/Symbol \136},0}/r_{j,0})' w l lw 3 dt 9 lc 9,\
	#'../../Code/quant_FF_norm.dat' u 1:($7**2) t '{/Symbol d}^2/{/Symbol d}^2_{0}' w l lw 3 dt 8 lc 8
	

set terminal postscript color enhanced dl 1.2 lw 2 'Helvetica' 20
set output 'quant_FF.ps'
replot
!./fixbb quant_FF.ps

