#!/bin/bash

# Exit the script when pressing Ctrl+C
trap "exit" SIGINT

# Run the program
ifort -w -no-wrap-margin -fast controls.f90 globals.f90 trajectory.f90 evolution.f90 losses.f90	\
							   edistribution.f90 synrad.f90 icrad.f90 opac_ff.f90 opac_gg.f90 	\
							   occultation.f90 fluxes.f90 JetStar.f90 -o JetStar.x
time ./JetStar.x

# Erase intermediate files
shopt -s extglob
rm NE_cuts.dat SED!(*_final_*).dat
