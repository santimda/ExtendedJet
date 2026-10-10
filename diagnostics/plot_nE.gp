reset
set xlabel "E' [eV]"
set ylabel "E'^{2}N'(E') [erg/cm^3]"
set logscale xy
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1e6:1e13]
set yrange [1e-20:]
#set yrange [1e30:1e37]
set cbtics 100
set key noautotitles top right Right spacing 1.35 width -2

eV = 1.602176e-12

plot '../Code/Ne_cuts_ord.dat' u ($1/eV):(($1**2)*$2):3 w l dt 1 lw 1 lc palette

set term postscript color enhanced dl 1.3 lw 2 'Helvetica' 20
set output 'nE.ps'
replot
!./fixbb nE.ps
