set logscale xy
set yrange [1e-4:10]

plot 'Quantities.dat' u 1:2 w l lw 2, 'Quantities.dat' u 1:3 w l lw 2, 'Quantities.dat' u 1:4 w l lw 2, 'Quantities.dat' u 1:5 w l lw 2, 'Quantities.dat' u 1:6 w l lw 2
