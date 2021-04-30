reset
set xlabel '{/Symbol e} [eV]'
set ylabel '{/Symbol e} L_{/Symbol e} [erg/s]'
set logscale xy
set format x '10^{%T}'
set format y '10^{%T}'
set xrange [1e0:1e14]
set yrange [1e25:1e34]
set key noautotitles top left Left spacing 1.35 width -2

plot 'SEDIC_final_jet.dat' u 1:2 w l dt 2 lw 3 lc 2 t 'l_j < 0.1D_{*}' ,\
'SEDIC_final_jet.dat' u 1:3 w l dt 3 lw 3 lc 2 t 'l_j < D_{*}' ,\
'SEDIC_final_jet.dat' u 1:4 w l dt 3 lw 2 lc rgb 'black' ,\
'SEDIC_final_jet.dat' u 1:5 w l dt 4 lw 2 lc rgb 'black' ,\
'SEDIC_final_jet.dat' u 1:6 w l dt 1 lw 3 lc 2 t 'Total' ,\
'SEDIC_final_cuts.dat' u 1:4:5 w l dt 1 lw 1 lc palette

set term postscript color enhanced dl 1.3 lw 2 'Helvetica' 20
set output 'SED_IC.ps'
replot
!./fixbb SED_IC.ps
