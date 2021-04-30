#!/bin/bash

############################### SCRIPT FOR THE ENERGY DISTRIBUTION COMPUTATION ###############################

# Exit the script when pressing Ctrl+C
trap "exit" SIGINT

# Disable computation for different orbital phases
sed -i "s/Norb = [0-9]*/Norb = 0/" JetStar.f90

# Set all the control options to NO, except for the ones related with the energy disribution
sed -i "s/ = yes/ = no/" controls.f90
sed -i "s/curvedJet = no/curvedJet = yes/" controls.f90
sed -i "s/NE_calc = no/NE_calc = yes/" controls.f90

# Set the value of the jet Lorentz factor
gamma=2.0
sed -i "s/gj = [0.-9.]*/gj = $gamma/" globals.f90

for etaB in 0.01 1.0; do

	echo -e '\n======= fB = '$etaB' ======='

	# Set the corresonding value of the fraction of magnetic energy in the globals file
	sed -i "s/fB = [0.-9.]*/fB = $etaB/" globals.f90
	
	# Compile and run the program
	ifort -w -no-wrap-margin -fast controls.f90 globals.f90 trajectory.f90 evolution.f90 losses.f90	\
								   edistribution.f90 synrad.f90 icrad.f90 opac_ff.f90 opac_gg.f90 	\
								   occultation.f90 fluxes.f90 JetStar.f90 -o JetStar.x
	./JetStar.x
	
	# Generate plots and copy them to the results folder
	gnuplot plot_NE.gp
	cp NE.ps Results/NE_g$gamma'_fB'$etaB.ps
	
done

# Erase intermediate files
rm NE_cuts.dat
	
	

