#!/bin/bash

############################### SCRIPT FOR THE ORBITAL VARIABILITY COMPUTATION ###############################

# Exit the script when pressing Ctrl+C
trap "exit" SIGINT

#### Helical jets ####
# Number of orbital phases considered between 0 and 1
Norb=21
sed -i "s/orbPhase = [0.-9.][0.-9.]*/orbPhase = 0/" JetStar.f90
sed -i "s/Norb = [0-9]*/Norb = $Norb/" JetStar.f90

# Set all the needed control options to YES
sed -i "s/ = no/ = yes/" controls.f90
sed -i "s/NE_calc = yes/NE_calc = no/" controls.f90

# Set the value of the fraction of magnetic energy density
etaB=0.01
sed -i "s/fB = [0.-9.]*/fB = $etaB/" globals.f90

for gamma in 1.5 2.0 3.0; do
	
	# Set the corresonding value of the jet Lorentz factor in the globals file and compute N(E) if necesary
	sed -i "s/gj = [0.-9.]*/gj = $gamma/" globals.f90	
	sed -i "s/NEint_calc = no/NEint_calc = yes/" controls.f90

	for i in 30 60; do

		echo -e '\n==============================================='
		echo ' Curved jets: Gamma = '$gamma' ; i = '$i'º'
		echo '==============================================='

		# Set the corresonding value of the system inclination in the globals file
		sed -i "s/inc = [0.-9.]*/inc = $i/" globals.f90
		
		# Compile and run the program
		ifort -w -no-wrap-margin -fast controls.f90 globals.f90 trajectory.f90 evolution.f90 losses.f90	\
									   edistribution.f90 synrad.f90 icrad.f90 opac_ff.f90 opac_gg.f90 	\
									   occultation.f90 fluxes.f90 JetStar.f90 -o JetStar.x
		./JetStar.x
		
		# Only compute the absorption tables in the first program run, and N(E) when gamma changes
		sed -i "s/remake_table_gg = yes/remake_table_gg = no/" controls.f90
		sed -i "s/remake_table_ff = yes/remake_table_ff = no/" controls.f90
		sed -i "s/NEint_calc = yes/NEint_calc = no/" controls.f90
		
	done
done

#### Vertical jets ####
sed -i "s/curvedJet = yes/curvedJet = no/" controls.f90

# Number of orbital phases considered between 0 and 0.5
Norb=11
sed -i "s/Norb = [0-9]*/Norb = $Norb/" JetStar.f90

for gamma in 1.5 2.0 3.0; do
	
	# Set the corresonding value of the jet Lorentz factor in the globals file and compute N(E) if necesary
	sed -i "s/gj = [0.-9.]*/gj = $gamma/" globals.f90	
	sed -i "s/NEint_calc = no/NEint_calc = yes/" controls.f90

	for i in 30 60; do

		echo -e '\n==============================================='
		echo ' Vertical jets: Gamma = '$gamma' ; i = '$i'º'
		echo '==============================================='

		# Set the corresonding value of the system inclination in the globals file
		sed -i "s/inc = [0.-9.]*/inc = $i/" globals.f90
		
		# Compile and run the program
		ifort -w -no-wrap-margin -fast controls.f90 globals.f90 trajectory.f90 evolution.f90 losses.f90	\
									   edistribution.f90 synrad.f90 icrad.f90 opac_ff.f90 opac_gg.f90 	\
									   occultation.f90 fluxes.f90 JetStar.f90 -o JetStar.x
		./JetStar.x
		
	done
	
	# Generate plots and copy them to the results folder
	cd Gamma$gamma
	gnuplot plot_SED_orb.gp
	gnuplot multiplot_fluxes.gp
	cp SED_orb.ps ../Results/SED_orb_g$gamma'_i60.ps'
	cp Fluxes.ps ../Results/Fluxes_g$gamma.ps
	cd ..
	
done

# Erase intermediate files
shopt -s extglob
rm SED!(*_final_*).dat
