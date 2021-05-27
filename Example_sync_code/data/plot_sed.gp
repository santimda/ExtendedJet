set encoding iso_8859_1
set key top left
set xrange [-7:14]
set yrange [28:]
set xtics 3; set mxtics 3
set ytics 1; set mytics 2
set xlabel 'log_{10}({/Symbol e} [eV])'
set ylabel 'log_{10}({/Symbol e} F({/Symbol e}) [erg/s/cm2/erg])'
eV = 1.602176e-12 

plot '../radsyn.dat' u (log10($1)):(log10($1*eV*$2)) t 'Sync' w l lw 3 lc rgb "#DC143C",\
	'../radssc.dat' u (log10($1)):(log10($1*eV*$2)) t 'SSC' w l lw 3 dt 2 lc rgb "blue",\

set terminal postscript colour enhanced dl 1.5 lw 2. "Helvetica" 21
set output 'sed.eps'
replot
!./fixbb sed.eps
