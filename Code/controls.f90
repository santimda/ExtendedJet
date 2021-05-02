MODULE controls
    implicit none
    
    logical, parameter :: no = .FALSE.
    logical, parameter :: yes = .TRUE.

    !===================================================
    ! Programs to run

    logical, parameter :: jet_calc = yes
    
    logical, parameter :: NE_calc = yes
    logical, parameter :: NEint_calc = yes
    
    logical, parameter :: radsyn_calc = yes
    logical, parameter :: radic_calc = yes
    
    !===================================================

END MODULE controls
