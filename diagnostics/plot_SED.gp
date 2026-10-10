reset
set xlabel '{/Symbol e} [eV]'
set ylabel '{/Symbol e} L_{/Symbol e} [erg/s]'
set logscale xy
set format x '10^{%T}'
set format y '10^{%T}'
set xrange [1e-6:1e14]
set yrange [1e28:1e45]
set xtics 1e2
set key noautotitles top left Left spacing 1.35 width -2

plot '../Code/SEDe_syn_jet.dat' u 1:3 w l dt 3 lw 2 lc 1 ,\
'../Code/SEDe_syn_jet.dat' u 1:4 w l dt 2 lw 1 lc rgb 'black' ,\
'../Code/SEDe_syn_jet.dat' u 1:7 w l dt 1 lw 2 lc 1 t 'Jet' ,\
'SEDIC_final_jet.dat' u 1:3 w l dt 3 lw 2 lc 1 ,\
'SEDIC_final_jet.dat' u 1:4 w l dt 2 lw 1 lc rgb 'black' ,\
'SEDIC_final_jet.dat' u 1:6 w l dt 1 lw 2 lc 1 ,\
'SEDsyn_final_CJ.dat' u 1:3 w l dt 3 lw 2 lc 2 ,\
'SEDsyn_final_CJ.dat' u 1:4 w l dt 2 lw 1 lc rgb 'black' ,\
'SEDsyn_final_CJ.dat' u 1:7 w l dt 1 lw 2 lc 2 t 'Counter-jet' ,\
'SEDIC_final_CJ.dat' u 1:3 w l dt 3 lw 2 lc 2 ,\
'SEDIC_final_CJ.dat' u 1:4 w l dt 2 lw 1 lc rgb 'black' ,\
'SEDIC_final_CJ.dat' u 1:6 w l dt 1 lw 2 lc 2

set term postscript color enhanced dl 1.3 lw 1 'Helvetica' 20
set output 'SED.ps'
replot
!./fixbb SED.ps
