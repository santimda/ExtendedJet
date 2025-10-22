MODULE controls
    implicit none
    
    logical, parameter :: no = .FALSE.
    logical, parameter :: yes = .TRUE.

    !===================================================
    ! Programs to run

    logical, parameter :: run_jet_calc = yes
    logical, parameter :: run_sed_calc = no!yes

    logical, parameter :: JET_ = yes
    logical, parameter :: CJ_ = yes

    logical, parameter :: jet_calc = yes
    logical, parameter :: inner_radius_calc = yes
    logical, parameter :: NE_calc = yes
    logical, parameter :: Np_calc = yes
    logical, parameter :: radsyn_calc = yes
    logical, parameter :: radpsyn_calc = yes
    logical, parameter :: radssc_calc = no
    logical, parameter :: radic_calc = no!yes
    
    !===================================================

END MODULE controls
