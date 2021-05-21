set encoding iso_8859_1
set title ''
set key box height 1.1 bottom left
set ylabel 'log_{10}(E^2 N(E) [erg])'
set xlabel 'log_{10}(E [eV])'
set xrange [6:14]
set yrange [35:43]
eV = 1.602176e-12

plot '../dist_e.dat' u (log10($1)):(log10((($1)*eV)**2*($2))) w l lw 2 lt 1 lc rgb "red" t 'N(E)'

set terminal postscript color enhanced lw 2.5 "Helvetica" 20
set output 'dist_e.ps'
replot
!./fixbb dist_e.ps