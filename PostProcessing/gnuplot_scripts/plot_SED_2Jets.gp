reset
set xlabel '{/Symbol e} [eV]'
set ylabel '{/Symbol e} L_{/Symbol e} [erg/s]'
set logscale xy
set format x '10^{%T}'
set format y '10^{%T}'
set xrange [1e-4:1e14]
set yrange [4e28:2e33]
set xtics 1e2
set key noautotitles top left Left spacing 1.35 width -6

# Gamma = 1.5
#plot 'SEDsyn_final_jet_g1.5.dat' u 1:7 w l dt 3 lw 2 lc 1 t 'Jet; {/Symbol g} = 1.5' ,\
#'SEDIC_final_jet_g1.5.dat' u 1:6 w l dt 3 lw 2 lc 1 ,\
#'SEDsyn_final_CJ_g1.5.dat' u 1:7 w l dt 3 lw 2 lc 2 t 'Counter-jet; {/Symbol g} = 1.5' ,\
#'SEDIC_final_CJ_g1.5.dat' u 1:6 w l dt 3 lw 2 lc 2

# Gamma = 2.0
plot 'SEDsyn_final_jet_g2.0.dat' u 1:7 w l dt 1 lw 2 lc 1 t 'Jet; {/Symbol g} = 2' ,\
'SEDIC_final_jet_g2.0.dat' u 1:6 w l dt 1 lw 2 lc 1 ,\
'SEDsyn_final_CJ_g2.0.dat' u 1:7 w l dt 1 lw 2 lc 2 t 'Counter-jet; {/Symbol g} = 2' ,\
'SEDIC_final_CJ_g2.0.dat' u 1:6 w l dt 1 lw 2 lc 2

# Gamma = 3.0
replot 'SEDsyn_final_jet_g3.0.dat' u 1:7 w l dt 2 lw 2 lc 1 t 'Jet; {/Symbol g} = 3' ,\
'SEDIC_final_jet_g3.0.dat' u 1:6 w l dt 2 lw 2 lc 1 ,\
'SEDsyn_final_CJ_g3.0.dat' u 1:7 w l dt 2 lw 2 lc 2 t 'Counter-jet; {/Symbol g} = 3' ,\
'SEDIC_final_CJ_g3.0.dat' u 1:6 w l dt 2 lw 2 lc 2

set term postscript color enhanced dl 1.3 lw 1 'Helvetica' 20
set output 'SED_2Jets.ps'
replot
!./fixbb SED_2Jets.ps
