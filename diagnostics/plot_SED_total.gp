set xlabel '{/Symbol e} [eV]'
set ylabel '{/Symbol e} L_{/Symbol e} [erg s^{-1}]'
set logscale xy
set format x '10^{%T}'
set format y '10^{%T}'
set xrange [1e3:1e16]
set yrange [1e38:1e44]
set xtics 1e1

#set key noautotitles top left Left spacing 1.35 width -2



plot '../Code/SEDsyne_total.dat' u 1:2 w l dt 3 lw 2 lc 2 t 'Syn-Jet(unabs)',\
	'../Code/SEDe_syn_jet.dat' u 1:9 w l dt 1 lw 3 lc 2 t 'Syn-Jet(abs)' ,\
	'../Code/SEDsyn_jetCJ.dat' u 1:8 w l dt 3 lw 2 lc 4 t 'Syn-CJ(unabs)',\
	'../../Codse/SEDsyn_jetCJ.dat' u 1:9 w l dt 1 lw 3 lc 4 t 'Syn-CJ(abs)',\
	'../Code/SEDssc_jet.dat' u 1:2 w l dt 1 lw 2 lc 6 t 'SSC-Jet',\
	'../Code/SEDssc_jetCJ.dat' u 1:2 w l dt 1 lw 3 lc 7 t 'SSC-CJ'
	
set arrow from 0, 5*10**39 to 1*10**9 ,5*10**39 nohead lc rgb "#000000"

set term postscript color enhanced dl 1.3 lw 1 'Helvetica' 20
set output 'SED_total.eps'
replot
!./fixbb SED_total.eps
