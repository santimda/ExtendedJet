MODULE controls
    implicit none

    logical, parameter :: no = .FALSE.
    logical, parameter :: yes = .TRUE.

    !===================================================
    ! Programs to run (compile-time; not part of the namelist)

    logical, parameter :: run_jet_calc = yes

    logical, parameter :: jet_calc = yes
    logical, parameter :: inner_radius_calc = yes

    ! SSC is not supported in this version: kept as a compile-time parameter
    ! (not part of the &run namelist) so it can never be switched on at runtime.
    logical, parameter :: radssc_calc = no

    !===================================================
    ! Runtime switches (defaults here; overridden by the &run namelist group
    ! read by system_parameters::read_namelist). All default to "do nothing"
    ! so that running without a namelist file is a safe no-op.

    logical :: JET_ = .false.
    logical :: CJ_ = .false.
    logical :: NE_calc = .false.
    logical :: Np_calc = .false.
    logical :: radsyn_calc = .false.
    logical :: radpsyn_calc = .false.
    logical :: radic_calc = .false.

    logical :: print_tlosses = .true.
    logical :: sort_outputs = .false.

    ! Physics-fix flags (awaiting author's decision). Reading them does not
    ! change any behaviour yet; see the TODO comments at their call sites.
    logical :: fix_norm_taper = .false.
    logical :: fix_inj_emin = .false.

    !===================================================

END MODULE controls
