#!/bin/bash

#################################### SCRIPT FOR THE TRAJECTORY COMPUTATION ###################################

# Exit the script when pressing Ctrl+C
trap "exit" SIGINT

# Set all the needed control options to YES
sed -i "s/ = yes/ = no/" controls.f90
sed -i "s/curvedJet = no/curvedJet = yes/" controls.f90

# Set the value of the orbital phase
phase=0
sed -i "s/orbPhase = [0.-9.][0.-9.]*/orbPhase = $phase/" JetStar.f90
sed -i "s/Norb = [0-9]*/Norb = 0/" JetStar.f90

for gamma in 1.5 2.0 3.0; do
	
	# Set the corresponding value of the jet Lorentz factor in the globals file
	sed -i "s/gj = [0.-9.]*/gj = $gamma/" globals.f90	

	echo -e '\n======= Gamma = '$gamma' ======='

	# Compile and run the program
	ifort -w -no-wrap-margin -fast controls.f90 globals.f90 trajectory.f90 evolution.f90 losses.f90	\
								   edistribution.f90 synrad.f90 icrad.f90 opac_ff.f90 opac_gg.f90 	\
								   occultation.f90 fluxes.f90 JetStar.f90 -o JetStar.x
	./JetStar.x
		
	# Rename the results
	mv Trajectory.dat Trajectory_g$gamma.dat
	mv TrajectoryCJ.dat TrajectoryCJ_g$gamma.dat
		
done

# Generate plot and copy them to the results folder
gnuplot plot_trajectory.gp
cp Trajectory.ps Results/Trajectory.ps
