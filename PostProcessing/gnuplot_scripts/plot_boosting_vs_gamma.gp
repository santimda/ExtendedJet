reset
set encoding iso_8859_1
#set tics scale 2
set key width 0.95 height 1.0 bottom left #box opaque 
set key samplen 3.5 spacing 1.3 font ',20'
#set format x '10^{%L}'
#set format y '10^{%L}'
set logscale
set xrange [0.1:50]
set yrange [0.005:60]
set mxtics 10
set mytics 10
set xlabel 'l/a'
set ylabel '{/Symbol d}^3_{obs}/{/Symbol g}_j'

f(x) = 1

plot f(x) t '' w l lw 1 lc 0,\
	'./Deltas_2.0_30.dat' u 3:4 t '{/Symbol g}_j=2.0' w l lw 3 lc 1,\
	'./Deltas_3.0_30.dat' u 3:4 t '{/Symbol g}_j=3.0' w l lw 3 lc 2,\
	'./Deltas_2.0_60.dat' u 3:4 t '' w l lw 3 lc 1 dt 2,\
	'./Deltas_3.0_60.dat' u 3:4 t '' w l lw 3 lc 2 dt 2
	
	
#	'./Deltas_1.5_30.dat' u 3:4 t '{/Symbol g}_j=1.5' w l lw 3 lc 1,\
#	'./Deltas_1.5_60.dat' u 3:4 t '' w l lw 3 lc 1 dt 2,\

set terminal postscript color enhanced dl 1.2 lw 2 'Helvetica' 20
set output 'boosting_vs_gamma.ps'
replot
!./fixbb boosting_vs_gamma.ps

