set encoding utf8  #iso_8859_1
set xlabel "E' [eV]"; set ylabel offset 0.0,-0.1
set ylabel "E'^{2}N'(E') [erg]"; set ylabel offset 0.1,0.0
set logscale xy
set format y '10^{%T}'
set format x '10^{%T}'
set xrange [1.01e6:1e14]
set yrange [1e38:1e49]
set key noautotitles top right Right spacing 1.35 width -2
#set ytics 2 ; set mytics 2
#set xtics 1 ; set mxtics 2


set palette maxcolors 6 # 6 = N_bins
set cblabel '' 
set cbtics ("I1" 1.4, "I2" 2.2,"I3" 3.,"I4" 3.8,"I5" 4.6,"I6" 5.4)
set palette defined ( 0 "#FF8C00", 1 "#8B0000" ) # From dark-orange to dark-red
#set datafile missing

eV = 1.602176e-12

plot '../output/Ne_bin.dat' u 2:((($2*eV)**2)*$3):1 w l lw 3.0 dt 1 lc palette

set terminal postscript color solid enhanced dl 1.5 lw 2.5 'Helvetica' 21
set output 'NE_bins.ps'
replot
!./fixbb NE_bins.ps
