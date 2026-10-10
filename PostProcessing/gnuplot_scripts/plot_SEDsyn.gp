reset
set xlabel '{/Symbol e} [eV]'
set ylabel '{/Symbol e} L_{/Symbol e} [erg/s]'
set logscale xy
set format x '10^{%T}'
set format y '10^{%T}'
set xrange [2e-8:1e9]
set yrange [1e28:1e39]
#set yrange [1e32:]
set key noautotitles bottom right spacing 1.35 width -2

plot '../../Code/SEDe_syn_cuts_ord.dat' u 1:2:3 w l dt 1 lw 1 lc palette,\
	'../../Code/SEDe_syn_cutsCJ_ord.dat' u 1:2:3 w l dt 1 lw 1 lc palette

#	'../../Code/SEDe_syn_jet.dat' u 1:2 w l dt 2 lw 3 lc 2 t 'z < 0.01 z_{max}' ,\
#	'../../Code/SEDe_syn_jet.dat' u 1:3 w l dt 3 lw 3 lc 2 t 'z < 0.1 z_{max}' ,\
#	'../../Code/SEDe_syn_jet.dat' u 1:4 w l dt 3 lw 2 lc rgb 'black' t 'Total(unabs)',\
#	'../../Code/SEDe_syn_jet.dat' u 1:($9/$1) w l dt 1 lw 3 lc 2 t 'Total(abs)' ,\
#   '../../Code/SEDe_syn_jet.dat' u 1:9 w l dt 1 lw 3 lc 2 t 'Total(abs2)'

set term postscript color enhanced dl 1.3 lw 2 'Helvetica' 20
set output 'SED_syn.ps'
replot
!./fixbb SED_syn.ps
