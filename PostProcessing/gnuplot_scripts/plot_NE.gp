reset
set xlabel "E' [eV]"
set ylabel "E'^{2}N'(E') [erg]"
set logscale xy
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1e6:1e13]
set yrange [1e38:]
#set yrange [1e30:1e37]
set cbtics 300
set key noautotitles top right Right spacing 1.35 width -2

eV = 1.602176e-12

plot '../../Code/NE_jet.dat' u 1:2 w l dt 3 lw 2 lc 2 t 'z < 0.01 z_{max}' ,\
	'../../Code/NE_jet.dat' u 1:3 w l dt 2 lw 2 lc 2 t 'z < 0.1 z_{max}' ,\
	'../../Code/NE_jet.dat' u 1:4 w l dt 1 lw 2 lc 2 t 'Total' ,\
	'../../Code/NE_cuts_ord.dat' u ($1/eV):(($1**2)*$2):3 w l dt 1 lw 1 lc palette
#for [i=1:Ni] 'NE_cuts.dat' every Ni::i u 1:2:(i) w l dt 1 lc palette

set term postscript color enhanced dl 1.3 lw 2 'Helvetica' 20
set output 'NE.ps'
replot
!./fixbb NE.ps
