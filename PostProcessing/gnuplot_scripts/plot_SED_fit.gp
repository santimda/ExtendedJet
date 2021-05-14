reset
set ylabel 'log_{10}({/Symbol e} L_{/Symbol e} [erg s^{-1}])' offset 0.25,0
set xlabel 'log_{10}({/Symbol e} [eV])'
set xrange [-6:14]
set yrange [25:34]
set xtics 2; set mxtics 2
set key noautotitles top left Left spacing 1.35 width -2


# Fit de la SED del jet

x1=5.5
x1_cut = 8
g(x)=b0+b1*(x-x1)+b2*(x-x1)**2+b3*(x-x1)**3+b4*(x-x1)**4+b5*(x-x1)**5+b6*(x-x1)**6+b7*(x-x1)**7+b8*(x-x1)**8+b9*(x-x1)**9+b10*(x-x1)**10
fit g(x) "< awk '($1 >2.e-4 && $1<2.e8) { print $1, $7 }' ./SEDsyn_final_jet.dat" u (log10($1)):(log10($2)) via b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10
g1(x)=x>x1_cut? g(x1_cut)*exp(-(x-x1_cut)**2.) : x<-3.8?0 : g(x)

x0=11.5
x0_cut=12.5
x0_low=6
f(x)=a0+a1*(x-x0)+a2*(x-x0)**2+a3*(x-x0)**3+a4*(x-x0)**4+a5*(x-x0)**5+a6*(x-x0)**6+a7*(x-x0)**7+a8*(x-x0)**8+a9*(x-x0)**9+a10*(x-x0)**10+a11*(x-x0)**11
fit f(x) "< awk '($1>1.e6 && $1<7.e12) { print $1, $6 }' ./SEDIC_final_jet.dat" u (log10($1)):(log10($2)) via a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
f1(x)=x>x0_cut? f(x0_cut)*exp(-(x-x0_cut)) : x<x0_low? f(x0_low)*exp(-(x0_low-x)) : f(x)

# Fit de la SED del counter-jet

x4=5.5
x4_cut=8
i(x)=d0+d1*(x-x4)+d2*(x-x4)**2+d3*(x-x4)**3+d4*(x-x4)**4+d5*(x-x4)**5+d6*(x-x4)**6+d7*(x-x4)**7+d8*(x-x4)**8+d9*(x-x4)**9+d10*(x-x4)**10
fit i(x) "< awk '($1 >5.e-4 && $1<2.e8) { print $1, $7 }' ./SEDsyn_final_CJ.dat" u (log10($1)):(log10($2)) via d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10
i1(x)=x>x4_cut? i(x4_cut)*exp(-(x-x4_cut)**2.) : x<-3.5?0 : i(x)

x3=9.5
x3_cut=12.5
x3_low=6
h(x)=c0+c1*(x-x3)+c2*(x-x3)**2+c3*(x-x3)**3+c4*(x-x3)**4+c5*(x-x3)**5+c6*(x-x3)**6+c7*(x-x3)**7+c8*(x-x3)**8+c9*(x-x3)**9+c10*(x-x3)**10
fit h(x) "< awk '($1>1.e6 && $1<7.e12) { print $1, $6 }' ./SEDIC_final_CJ.dat" u (log10($1)):(log10($2)) via c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10
h1(x)=x>x3_cut? h(x3_cut)*exp(-(x-x3_cut)) : x<x3_low? h(x3_low)*exp(-(x3_low-x)) : h(x)

#-----------------------------------------------------

# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 está en GHz y Sv en $2 en mJy
dist = 3.0*3.0856e21
a = 4.0*pi*dist**2
h = 6.6260755e-27
eV = 1.602176e-12

#-----------------------------------------------------
plot log10(10**f1(x)+10**g1(x)) ls 1 dt 1 lw 3 lc 7 t 'Jet',\
	 log10(10**h1(x)+10**i1(x)) ls 5 dt 1 lw 3 lc 4 t 'Counter-jet',\
	'SEDsyn_final_jet.dat' u (log10($1)):(log10($7)) w l dt 2 lw 2 lc 1 notitle,\
	'SEDsyn_final_jet.dat' u (log10($1)):(log10($4)) w l dt 4 lw 1 lc rgb 'black',\
	'SEDIC_final_jet.dat' u (log10($1)):(log10($6)) w l dt 3 lw 2 lc 1 notitle,\
	'SEDIC_final_jet.dat' u (log10($1)):(log10($4)) w l dt 4 lw 1 lc rgb 'black',\
	'SEDsyn_final_CJ.dat' u (log10($1)):(log10($7)) w l dt 2 lw 2 lc 2 notitle,\
	'SEDsyn_final_CJ.dat' u (log10($1)):(log10($4)) w l dt 4 lw 1 lc rgb 'black',\
	'SEDIC_final_CJ.dat' u (log10($1)):(log10($6)) w l dt 3 lw 2 lc 2 notitle,\
 	'SEDIC_final_CJ.dat' u (log10($1)):(log10($4)) w l dt 4 lw 1 lc rgb 'black',\
	'windemiff.dat' u (log10((($1)*h/eV))):(log10(($1)*a*($2)*1e-26)) w l lw 2 t 'wind f-f'

#	'SEDsyn_jet.dat' u (log10($1)):(log10($4)) w l dt 2 lw 1 lc rgb 'black',\
#	'SEDIC_final_jet.dat' u (log10($1)):(log10($4)) w l dt 4 lw 1 lc rgb 'black',\
#	'SEDsyn_CJ.dat' u (log10($1)):(log10($4)) w l dt 2 lw 1 lc rgb 'black',\
#	'SEDIC_final_CJ.dat' u (log10($1)):(log10($4)) w l dt 4 lw 1 lc rgb 'black',\
#	log10(10**f1(x)+10**g1(x)) ls 1 dt 1 lw 3 lc 7 t 'Jet', f1(x) ls 2 lc 7 dt 2 t '', g1(x) ls 2 lc 7 dt 4 t '',\
#	log10(10**h1(x)+10**i1(x)) ls 5 dt 1 lw 3  lc 7 t 'Counter-jet', h1(x) ls 6 lc 7  dt 2 t '', i1(x) ls 6 lc 7  dt 4 t ''

set term postscript color enhanced dl 1.3 lw 2 'Helvetica' 20
set output 'SEDs_fit.ps'
replot
!./fixbb SEDs_fit.ps
