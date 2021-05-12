MODULE controls
    implicit none
    
    logical, parameter :: no = .FALSE.
    logical, parameter :: yes = .TRUE.

    !===================================================
    ! Programs to run

    logical, parameter :: jet_calc = yes
    
    logical, parameter :: NE_calc = yes    ! For plotting only
    logical, parameter :: NEint_calc = yes ! For precision in integrals
    
    logical, parameter :: radsyn_calc = yes
    logical, parameter :: radic_calc = no!yes

    logical, parameter :: map_calc = no!yes
    
    !===================================================

END MODULE controls
