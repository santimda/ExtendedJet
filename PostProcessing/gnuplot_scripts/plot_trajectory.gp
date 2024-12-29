reset
set size ratio 1./4.5
set yrange [-10:10]
set xrange [-45:45]
set ylabel 'x/a'
set xlabel 'z/a ; y/a'
set ytics 4
set xtics 10
set mxtics 2
#set key at 18,4

# Orbital separation and stellar radius
a = 3e12
Rstar = 7.4e11
set object 1 circle at 0,0 size Rstar/a front fc rgb 'blue' fs solid 1.0

plot '../../Code/trajectory.dat' u ($2/a):($1/a) w l dt 2 lw 2 lc 1 t 'x-y' ,\
	 '../../Code/trajectory.dat' u ($3/a):($1/a) w l dt 1 lw 2 lc 2 t 'x-z' ,\
	 'TrajectoryCJ.dat' u ($3/a):($1/a) w l dt 1 lw 2 lc 2 notitle

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output 'Trajectory.ps'
replot
!./fixbb Trajectory.ps
