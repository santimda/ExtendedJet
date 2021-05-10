MODULE controls
    implicit none
    
    logical, parameter :: no = .FALSE.
    logical, parameter :: yes = .TRUE.

    !===================================================
    ! Programs to run

    logical, parameter :: jet_calc = no!yes
    
    logical, parameter :: NE_calc = no!yes    ! For plotting only
    logical, parameter :: NEint_calc = no!yes ! For precision in integrals
    
    logical, parameter :: radsyn_calc = no!yes
    logical, parameter :: radic_calc = no!yes

    logical, parameter :: map_calc = yes
    
    !===================================================

END MODULE controls
