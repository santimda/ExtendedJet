module controles

implicit none
public

!===================================================
! Programs to run 
logical, parameter :: no = .FALSE.
logical, parameter :: yes = .TRUE.

logical, parameter :: write_files = no!yes 

logical, parameter :: dist_calc = yes
logical, parameter :: radsyn_calc = yes
logical, parameter :: radssc_calc = yes
logical, parameter :: radsyn_data = yes

!===================================================

END module controles
