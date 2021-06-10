#!/bin/bash

####################################### SCRIPT FOR THE SED COMPUTATION #######################################

# Exit the script when pressing Ctrl+C
trap "exit" SIGINT

# Set all the needed control options to YES
sed -i "s/ = no/ = yes/" controls.f90
sed -i "s/NE_calc = yes/NE_calc = no/" controls.f90

#### Jet and counter-jet contributions separately ####
# Set the value of the orbital phase, the system inclination and the fraction of magnetic energy density
i=60
etaB=0.01
phase=90
sed -i "s/orbPhase = [0.-9.][0.-9.]*/orbPhase = $phase/" JetStar.f90
sed -i "s/Norb = [0-9]*/Norb = 1/" JetStar.f90
sed -i "s/inc = [0.-9.]*/inc = $i/" globals.f90
sed -i "s/fB = [0.-9.]*/fB = $etaB/" globals.f90

for gamma in 1.5 2.0 3.0; do
	
	# Set the corresonding value of the jet Lorentz factor in the globals file
	sed -i "s/gj = [0.-9.]*/gj = $gamma/" globals.f90	

	echo -e '\n======= Gamma = '$gamma' ======='

	# Compile and run the program
	ifort -w -no-wrap-margin -fast controls.f90 globals.f90 trajectory.f90 evolution.f90 losses.f90	\
								   edistribution.f90 synrad.f90 icrad.f90 opac_ff.f90 opac_gg.f90 	\
								   occultation.f90 fluxes.f90 JetStar.f90 -o JetStar.x
	./JetStar.x
		
	# Only compute the absorption tables in the first program run
	sed -i "s/remake_table_gg = yes/remake_table_gg = no/" controls.f90
	sed -i "s/remake_table_ff = yes/remake_table_ff = no/" controls.f90
	
	# Rename the results
	mv SEDsyn_final_jet.dat SEDsyn_final_jet_g$gamma.dat
	mv SEDIC_final_jet.dat SEDIC_final_jet_g$gamma.dat
	mv SEDsyn_final_CJ.dat SEDsyn_final_CJ_g$gamma.dat
	mv SEDIC_final_CJ.dat SEDIC_final_CJ_g$gamma.dat
		
done

# Generate plots and copy them to the results folder
gnuplot plot_SED_2Jets.gp
cp SED_2Jets.ps Results/SED_2Jets_fB$etaB'_i'$i.ps

#### Total SED ####
# Set the value of the jet Lorentz factor
gamma=2.0
sed -i "s/gj = [0.-9.]*/gj = $gamma/" globals.f90

for i in 30 60; do
	
	# Set the corresonding value of the system inclination in the globals file
	sed -i "s/inc = [0.-9.]*/inc = $i/" globals.f90
	
	for etaB in 0.01 1.0; do

		# Set the corresonding value of the magnetic energy density faction in the globals file
		sed -i "s/fB = [0.-9.]*/fB = $etaB/" globals.f90

		echo -e '\n======= Gamma = '$gamma' ======='

		# Compile and run the program
		ifort -w -no-wrap-margin -fast controls.f90 globals.f90 trajectory.f90 evolution.f90 losses.f90	\
									   edistribution.f90 synrad.f90 icrad.f90 opac_ff.f90 opac_gg.f90 	\
									   occultation.f90 fluxes.f90 JetStar.f90 -o JetStar.x
		./JetStar.x
		
		# Rename the results
		mv SEDsyn_final.dat SEDsyn_final_fB$etaB.dat
		mv SEDIC_final.dat SEDIC_final_fB$etaB.dat
		
	done
	
	# Generate plots and copy them to the results folder
	gnuplot plot_SED_Total.gp
	cp SED_Total.ps Results/SED_Total_g$gamma'_i'$i.ps
		
done

# Erase intermediate files
shopt -s extglob
rm SED!(*_final_*).dat
