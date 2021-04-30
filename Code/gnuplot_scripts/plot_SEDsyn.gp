reset
set xlabel '{/Symbol e} [eV]'
set ylabel '{/Symbol e} L_{/Symbol e} [erg/s]'
set logscale xy
set format x '10^{%T}'
set format y '10^{%T}'
set xrange [1e-6:1e10]
set yrange [1e24:1e34]
set key noautotitles top left Left spacing 1.35 width -2

plot 'SEDsyn_final_jet.dat' u 1:2 w l dt 2 lw 3 lc 2 t 'l_j < 0.1D_{*}' ,\
'SEDsyn_final_jet.dat' u 1:3 w l dt 3 lw 3 lc 2 t 'l_j < D_{*}' ,\
'SEDsyn_final_jet.dat' u 1:4 w l dt 3 lw 2 lc rgb 'black' ,\
'SEDsyn_final_jet.dat' u 1:5 w l dt 4 lw 2 lc rgb 'black' ,\
'SEDsyn_final_jet.dat' u 1:6 w l dt 5 lw 2 lc rgb 'black' ,\
'SEDsyn_final_jet.dat' u 1:7 w l dt 1 lw 3 lc 2 t 'Total' ,\
'SEDsyn_final_cuts.dat' u 1:5:6 w l dt 1 lw 1 lc palette

set term postscript color enhanced dl 1.3 lw 2 'Helvetica' 20
set output 'SED_syn.ps'
replot
!./fixbb SED_syn.ps
