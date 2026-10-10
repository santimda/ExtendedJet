reset
set xlabel '{/Symbol e} [eV]'
set ylabel '{/Symbol e} L_{/Symbol e} [erg/s]'
set logscale xy
set format x '10^{%T}'
set format y '10^{%T}'
#set xrange [1e-6:1e12]
#set yrange [1e24:1e34]
set yrange [1e32:]
set key noautotitles top right spacing 1.35 width -2

plot '../output/SEDssc_jet.dat' u 1:3 w l dt 2 lw 3 lc 2 t 'z < 10 z_{0}',\
	'../output/SEDssc_jet.dat' u 1:4 w l dt 3 lw 3 lc 2 t 'z < 100 z_{0}',\
	'../output/SEDssc_jet.dat' u 1:2 w l dt 1 lw 3 lc 2 t 'Total',\
	'../output/SEDssc_cuts.dat' u 1:2:3 w l dt 1 lw 1 lc palette

set term postscript color enhanced dl 1.3 lw 2 'Helvetica' 20
set output 'SED_ssc.ps'
replot
!./fixbb SED_ssc.ps
