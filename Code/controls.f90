MODULE controls
    implicit none
    
    logical, parameter :: no = .FALSE.
    logical, parameter :: yes = .TRUE.

    !===================================================
    ! Programs to run

    logical, parameter :: run_jet_calc = no
    logical, parameter :: run_sed_calc = yes

    logical, parameter :: jet_calc = yes
    logical, parameter :: NE_calc = yes
    logical, parameter :: radsyn_calc = yes
    logical, parameter :: radic_calc = no!yes
    logical, parameter :: map_calc = no!yes
    
    !===================================================

END MODULE controls
