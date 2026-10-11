MODULE controls
    implicit none

    logical, parameter :: off = .FALSE.
    logical, parameter :: on = .TRUE.

    !===================================================
    ! Programs to run (compile-time; not part of the namelist)

    logical, parameter :: run_jet_calc = on

    logical, parameter :: jet_calc = on
    logical, parameter :: inner_radius_calc = on

    ! SSC is not supported in this version: kept as a compile-time parameter
    ! (not part of the &run namelist) so it can never be switched on at runtime.
    logical, parameter :: radssc_calc = off

    !===================================================
    ! Runtime switches (defaults here; overridden by the &run namelist group
    ! read by system_parameters::read_namelist). All default to "do nothing"
    ! so that running without a namelist file is a safe no-op.

    logical :: JET_ = off
    logical :: CJ_ = off
    logical :: NE_calc = off
    logical :: Np_calc = off
    logical :: radsyn_calc = off
    logical :: radpsyn_calc = off
    logical :: radic_calc = off

    logical :: print_tlosses = on
    logical :: sort_outputs = off

    ! Physics-fix flags (awaiting author's decision). Reading them does not
    ! change any behaviour yet; see the TODO comments at their call sites.
    logical :: fix_norm_taper = off
    logical :: fix_inj_emin = off

    !===================================================

END MODULE controls
