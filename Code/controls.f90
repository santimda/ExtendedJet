MODULE controls
    implicit none
    
    logical, parameter :: no = .FALSE.
    logical, parameter :: yes = .TRUE.

    !===================================================
    ! Programs to run

    logical, parameter :: jet_calc = yes
    
    logical, parameter :: NE_calc = no
    logical, parameter :: NEint_calc = no!yes
    
    logical, parameter :: radsyn_calc = no!yes
    logical, parameter :: radic_calc = no
    
    !===================================================

END MODULE controls
