reset
#set encoding iso_8859_1
#set size square
#set title 'Dynamic evolution of the cooling dependences'
#set grid
#set tics scale 2
set key top left
set key samplen 3.5 spacing 1.4 font ',20'
set format x '10^{%L}'
set format y '10^{%L}'
set logscale
#set xrange [:8]
#set yrange [-1:1]
#set mxtics 10
set mytics 10
set xlabel 'l'
set ylabel 'E_{e,max}(l)'

plot '../output/Ee_max_FF.dat' u 1:3 t 'E_{max}' w l lw 3 lc 1

set terminal postscript color enhanced dl 1.2 lw 2 'Helvetica' 20
set output 'Ee_max_FF.ps'
replot
!./fixbb Ee_max_FF.ps
