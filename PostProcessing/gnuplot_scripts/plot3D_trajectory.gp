reset
set view 60, 45, 1, 2  # rot_x, rot_z, scale, scale_z
set xrange [-12.5:12.5]
set yrange [-12.5:12.5]
set zrange [0:50]
set xlabel 'x/a'
set ylabel 'y/a'
set zlabel 'z/a'
set ztics 5
set ytics 4
set xtics 4
set mxtics 2
set ticslevel 0
#set key at 18,4

# Orbital separation and stellar radius
a = 3e12
Rstar = 7.4e11
RCO = 1e11 # exagerated

set object 1 circle at 0,0,0 size Rstar/a front fc rgb 'blue' fs solid 1.0
set object 2 circle at 1,0,0 size RCO/a front fc rgb 'black' fs solid 1.0

splot 'Trajectory.dat' u ($1/a):($2/a):($3/a) w l dt 1 lw 2 lc 1 t ''

set term postscript color enhanced dl 1.5 lw 2.0 'Helvetica' 20
set output '3D_Trajectory.ps'
replot
!./fixbb 3D_Trajectory.ps
