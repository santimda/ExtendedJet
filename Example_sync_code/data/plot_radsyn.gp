set encoding iso_8859_1
set key top right
set xrange [-7.2:3]
set yrange [-7:]
set xtics 3; set mxtics 3
set ytics 1; set mytics 2
set xlabel 'log_{10}({/Symbol e} [eV])'
set ylabel 'log_{10}({/Symbol e} F({/Symbol e}) [erg/s/cm2/erg])'

plot '../radsyn.dat' u (log10($1)):(log10($2)) t 'Sync' w l lw 3 lc rgb "#DC143C",\
	'data.dat' u (log10($1)):(log10($2)):($3/$2) w yerrorbars pt 7 ps 1.5 t 'data'

set terminal postscript colour enhanced dl 2. lw 2. "Helvetica" 21
set output 'radsyn.ps'
replot
!./fixbb radsyn.ps
