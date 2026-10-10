PROGRAM main
    use constants
    use globals
    use controls
    use system_parameters
    implicit none
    real(dp) :: sigma0_value, Lj_value, z0_value, a_inj_value
    real(dp), allocatable :: E_data(:)
    real(dp) :: S_data(m_data_max)
    logical, parameter :: verbose = .True.
    character(len=512) :: nml_path

    interface
        subroutine run_jet(verbose) 
            logical, intent(in) :: verbose
        end subroutine run_jet
    end interface

    interface
        subroutine sed_model(sigma0_value, Lj_value, a_inj_value, E_vector, S_vector) 
            use globals
            use constants
            real(dp), intent(in) :: sigma0_value, Lj_value, a_inj_value
            real(dp), intent(in) :: E_vector(:)
            real(dp), intent(out) :: S_vector(m_data_max)
        end subroutine sed_model
    end interface 

    interface
        subroutine syn_model(sigma0_value, zmax_value, Lj_value, a_inj_value, &
                             zmax_inj_value, gammaIndex_value, geomIndex_out_value, &
                             z0_value, fnT_value, facc_value, p_value, E_vector, S_vector) 
            use globals
            use constants
            real(dp), intent(in) :: sigma0_value, zmax_value, Lj_value, a_inj_value
            real(dp), intent(in) :: zmax_inj_value, gammaIndex_value, geomIndex_out_value
            real(dp), intent(in) :: z0_value, fnT_value, facc_value, p_value
            real(dp), intent(in) :: E_vector(:)
            real(dp), intent(out) :: S_vector(m_data_max)
        end subroutine syn_model
    end interface 

    interface
        subroutine power_model(sigma0_value, zmax_value, Lj_value, a_inj_value, z_data, L_B_data, L_kin_data, L_e_data, L_e_2_data, L_sy_data, L_sy2_data)
            use globals
            use constants
            real(dp), intent(in) :: sigma0_value, zmax_value, Lj_value, a_inj_value, z_data(:)
            real(dp), intent(out) :: L_B_data(m_data_max), L_kin_data(m_data_max), L_e_data(m_data_max), L_e_2_data(m_data_max), L_sy_data(m_data_max), L_sy2_data(m_data_max)
        end subroutine power_model
    end interface

    !----------------------------------------------------------------
    ! Original Fortran program

    if (command_argument_count() < 1) then
        print *, 'Usage: ./run.x params.nml'
        stop 1
    endif
    call get_command_argument(1, nml_path)
    call read_namelist(nml_path)

    if (run_jet_calc) call run_jet(verbose)

    !----------------------------------------------------------------
    ! Short version to calculate specific fluxes integrated up to 
    ! certain jet heights. Used for fitting SEDs.
    
    ! if (run_sed_calc) then
    !     allocate( E_data(4) )
    !     E_data = (/1d-6*eV, 1d-2*eV, 1.d2*eV, 1.d5*eV/)
    !     sigma0_value = 1.d1
    !     Lj_value = L_Edd*1d-4
    !     a_inj_value = 0.5d0
    !     zmax_value = 1d6*R_g
    !     call syn_model(sigma0_value, zmax_value, Lj_value, a_inj_value, E_data, S_data) ! Energies in cgs
    !     call sed_model(sigma0_value, Lj_value, a_inj_value, E_data, S_data) ! Energies in cgs
    !     deallocate( E_data )
    ! endif


END PROGRAM main


!==================================================================

subroutine run_jet(verbose)
    use controls
    use constants
    use globals
    use initialize
    use trajectory
    use evolution
    use losses
    use edistribution
    use synrad
    use icrad
    use sscrad
    use system_parameters
    implicit none
    logical, intent(in) :: verbose
    
    !real(dp) :: EphSy(NSED_sy), Lsy(NSED_sy), Lsy_CJ(NSED_sy)
    !real(dp) :: EphIC(NSED_IC), LIC(NSED_IC), LIC_CJ(NSED_IC)
    !real(dp) :: Flux(N_bands)
    character(3) :: numNE, numSED

    write(numNE,'(I3)') nIn         !Number of points for the N(E) distribution
    write(numSED,'(I3)') NSED_sy    !Number of points for the synchrotron SED

    ! Calculate everything for the jet first
    if (JET_) then
        if (verbose) print *, '\n Initialize vectors'    
        call initialize_allocate(sigma_0, zmax, Lj, a_inj) !, z0 
        !Jet trajectory
        if (verbose .and. jet_calc) print *, '\n Calculating jet trajectory'
        if (jet_calc) call trajectory_calc('') !'' = jet 
        
        !Evolution of different quantities
        if (verbose .and. jet_calc) print *, '\n Calculating the evolution of quantities along the jet'
        if (jet_calc) call Quantities('') !'' = jet 

        !Particle distribution
        if (NE_calc) then 
            if (verbose) write(*,'(1X,A,I3,A)') '\n Computing the particle energy distribution (',nIn,' points)'
            call Edist('e')        ! Calculate the particle energy distribution
            if (sort_outputs) call execute_command_line('cd '//trim(outdir)//' && '//trim(scriptdir)//'/SortNe.sh '//trim(adjustl(numNE)))   !Sort N(E) file

            call Calc_Ue()      ! Calculate energy densities
            if (print_tlosses) then
                if (verbose) print *, '\n Creating the tloss files for plots'
                call Print_tloss(1, 'e')
                call Print_tloss(Ns/4, 'e')
                call Print_tloss(Ns/2, 'e')
                call Print_tloss(Ns, 'e')
            endif
        endif

        if (Np_calc) then
            call Edist('p')
            if (sort_outputs) call execute_command_line('cd '//trim(outdir)//' && '//trim(scriptdir)//'/SortNp.sh '//trim(adjustl(numNE)))   !Sort N(E) file
        endif

        !Synchrotron SED
        if (radsyn_calc) then
            if (verbose) write(*,'(1X,A,I3,A)') '\n Calculating the synchrotron SED (',NSED_sy,' points) in the jet'
            call SEDsyn('', 'e')      ! Calculate synchrotron SED
            if (sort_outputs) call execute_command_line('cd '//trim(outdir)//' && '//trim(scriptdir)//'/SortSEDesyn.sh '//trim(adjustl(numSED))) !Sort SED file
            !call Calc_Usyn('')   ! Calculate the energy density of synchrotron photons
        endif

        if (radpsyn_calc) then
            call SEDsyn('', 'p')
            if (sort_outputs) call execute_command_line('cd '//trim(outdir)//' && '//trim(scriptdir)//'/SortSEDpsyn.sh '//trim(adjustl(numSED))) !Sort SED file
        endif
        
        !IC SED     
        if (radic_calc) call SEDIC(iso)     ! Calculate the IC-CMB SED

        !SSC SED     
        if (radssc_calc) call SEDssc('')    ! Calculate the SSC SED
    endif

    ! Calculate everything for the counter-jet next
    if (CJ_) then
        if (.not. JET_) then
            if (verbose) print *, '\n Initialize vectors'    
            call initialize_allocate(sigma_0,zmax,Lj,a_inj) !z0,
        endif

        !Jet trajectory
        if (verbose .and. jet_calc) print *, '\n Calculating counter-jet trajectory'
        if (jet_calc) call trajectory_calc('CJ')
        
        !Evolution of different quantities
        if (verbose .and. jet_calc) print *, '\n Calculating the evolution of quantities along the counter-jet'
        if (jet_calc) call Quantities('CJ')

        !Particle distribution
        if (NE_calc .and. (.not. JET_)) then
            print*, 'ERROR: current version cannot calculate the electron distribution for the CJ, run with jet!' 
            stop
            if (verbose) write(*,'(1X,A,I3,A)') '\n Computing the particle energy distribution (',nIn,' points)'
            call Edist('e')        ! Calculate the particle energy distribution
            !call Edist('p')
            if (sort_outputs) call execute_command_line('cd '//trim(outdir)//' && '//trim(scriptdir)//'/SortNe.sh '//trim(adjustl(numNE)))   !Sort N(E) file
            call Calc_Ue()      ! Calculate energy densities
            if (print_tlosses) then
                if (verbose) print *, '\n Creating the tloss files for plots'
                call Print_tloss(1, 'e')
                call Print_tloss(Ns/3, 'e')
                call Print_tloss(Ns, 'e')
            endif
        endif

        !Synchrotron SED
        if (radsyn_calc) then
            if (verbose) write(*,'(1X,A,I3,A)') '\n Calculating the synchrotron SED (',NSED_sy,' points) in the counter-jet'
            call SEDsyn('CJ', 'e')      ! Calculate synchrotron SED
            if (sort_outputs) call execute_command_line('cd '//trim(outdir)//' && '//trim(scriptdir)//'/SortSEDesyn.sh '//trim(adjustl(numSED))) !Sort SED file
            !call Calc_Usyn('CJ')   ! Calculate the energy density of synchrotron photons
        endif

        if (radpsyn_calc) then
            call SEDsyn('CJ', 'p')
            if (sort_outputs) call execute_command_line('cd '//trim(outdir)//' && '//trim(scriptdir)//'/SortSEDpsyn.sh '//trim(adjustl(numSED))) !Sort SED file
        endif
        
        !IC SED     
        if (radic_calc) call SEDIC(iso)     ! Calculate the IC-CMB SED

        !SSC SED     
        if (radssc_calc) call SEDssc('CJ')  ! Calculate the SSC SED
    endif

    if (verbose) print *, '\n Deallocate vectors'    
    call initialize_deallocate
    
    if (verbose) print *, '\n Program finished'    

end subroutine run_jet



! !============================================================================
! ! sed_model(START)
! !============================================================================
! subroutine sed_model(sigma0_value, Lj_value, a_inj_value, E_data, S_data)
!   !===================================================================
!   ! Subroutine to calculate the synchrotron (+ SSC) SED at specific 
!   ! energies for the given input parameters. Calculates both J and CJ
!   ! TO DO: INCORPORATE SSC IN THE COMPUTATION!!!
!   !===================================================================
!     use controls
!     use constants
!     use globals
!     use initialize
!     use trajectory
!     use evolution
!     use losses
!     use edistribution
!     use synrad
!     use sscrad

!     implicit none
!     real(dp), intent(in) :: sigma0_value, Lj_value, a_inj_value
!     real(dp), intent(in) :: E_data(:)
!     real(dp), intent(out) :: S_data(m_data_max)
!     real(dp) :: S_data_J(m_data_max), S_data_CJ(m_data_max)
!     logical, parameter :: verbose = .False.
!     logical, parameter :: calc_SSC = .False.
        
!     S_data = 0.d0 ! Initialize
!     ! Calculate everything for the jet first
!     if (JET_) then  
!         call initialize_allocate(sigma0_value, Lj_value, a_inj_value) !z0_value, 

!         !Jet trajectory
!         call trajectory_calc('')  ! jet 
        
!         !Evolution of different quantities
!         call Quantities('')       ! jet 

!         !Particle distribution
!         call Edist('e')   

!         !Synchrotron SED at the given values  
!         call SEDsyn_data(E_data, S_data_J, '', 'e')  
!         S_data = S_data_J

!         if (calc_SSC) then
!         !SSC SED at the given values  
!             call SEDssc_data(E_data, S_data_J, '') 
!             S_data = S_data + S_data_J
!         endif   
!         if (verbose) print *, 'S_data[cgs] =', S_data

!         call initialize_deallocate 
!     endif 


!     ! Calculate everything for the counter-jet next
!     if (CJ_) then
!         call initialize_allocate(sigma0_value, Lj_value, a_inj_value) !z0_value, 

!         !Jet trajectory
!         call trajectory_calc('CJ')  ! counter-jet 
        
!         !Evolution of different quantities
!         call Quantities('CJ')       ! counter-jet 

!         !Particle distribution
!         if (.not. JET_) call Edist('e')  ! Calculate the particle energy distribution

!         !Synchrotron SED at the given values
!         call SEDsyn_data(E_data, S_data_CJ, 'CJ', 'e')
!         S_data = S_data + S_data_CJ

!         !SSC SED at the given values  
!         if (calc_SSC) then
!             call SEDssc_data(E_data, S_data_CJ, '') 
!             S_data = S_data + S_data_CJ  
!         endif 
!         if (verbose) print *, 'S_data[cgs] =', S_data
        
!         call initialize_deallocate 
!     endif

!     if (verbose) then
!         print *, 'E_data[eV] =', E_data/eV
!         print *, 'S_data[cgs] =', S_data
!     endif

! end subroutine sed_model

!===================================================================



!============================================================================
! syn_model(START)
!============================================================================
subroutine syn_model(sigma0_value, zmax_value, Lj_value, a_inj_value, &
                     zmax_inj_value, gammaIndex_value, geomIndex_out_value, &
                     z0_value, fnT_value, facc_value, p_value, E_data, S_data)
  !===================================================================
  ! Subroutine to calculate the synchrotron SED at specific energies
  ! for the given input parameters. 
  ! Calculates both J and CJ
  !===================================================================
    use controls
    use constants
    use globals
    use initialize
    use trajectory
    use evolution
    use losses
    use edistribution
    use synrad
    use sscrad
    implicit none
    real(dp), intent(in) :: sigma0_value, zmax_value, Lj_value, a_inj_value
    real(dp), intent(in) :: zmax_inj_value, gammaIndex_value, geomIndex_out_value
    real(dp), intent(in) :: z0_value, fnT_value, facc_value, p_value
    real(dp), intent(in) :: E_data(:)
    real(dp), intent(out) :: S_data(m_data_max)
    real(dp) :: S_data_J(m_data_max), S_data_CJ(m_data_max)
    logical, parameter :: verbose = .False.

    ! Set system parameters from input values (lengths go through their *_rg
    ! mirrors because derive_parameters() recomputes the cm values from them)
    call derive_parameters()   ! makes sure R_g is set
    zmax_inj_rg = zmax_inj_value/R_g
    gammaIndex = gammaIndex_value
    geomIndex_out = geomIndex_out_value
    z0_rg = z0_value/R_g
    fnT = fnT_value
    facc = facc_value
    p = p_value
        
    ! Calculate everything for the jet first
    call initialize_allocate(sigma0_value, zmax_value, Lj_value, a_inj_value) !z0_value, 
    call trajectory_calc('')  
    call Quantities('')    
    call Edist('e')   
    !Synchrotron SED at the given values
    call SEDsyn_data(E_data, S_data_J, '', 'e')   
    !call SEDsyn_data(E_data, S_data_J, '', 'e')  

    ! Calculate everything for the counter-jet next
    call trajectory_calc('CJ')  
    call Quantities('CJ')       
    call SEDsyn_data(E_data, S_data_CJ, 'CJ', 'e')
    call initialize_deallocate 

    !write(*,*) Ns, dlogz
    S_data = S_data_J + S_data_CJ

    if (verbose) then
        print *, 'E_data[eV] =', E_data/eV
        print *, 'S_data[cgs] =', S_data
    endif

end subroutine syn_model

!===================================================================

!============================================================================
! power_model(START)
!============================================================================
subroutine power_model(sigma0_value, zmax_value, Lj_value, a_inj_value, z_data, L_B_data, L_kin_data, L_e_data, L_e_2_data, L_sy_data, L_sy2_data)
  !===================================================================
  ! Subroutine to calculate the power components (magnetic, kinetic, electron)
  ! as a function of jet height for the given input parameters.
  ! Returns arrays of size Ns (defined in globals.f90, default 200)
  !===================================================================
    use controls
    use constants
    use globals
    use system_parameters
    use initialize
    use trajectory
    use evolution
    use edistribution
    use synrad

    implicit none
    real(dp), intent(in) :: sigma0_value, zmax_value, Lj_value, a_inj_value
    real(dp), intent(in) :: z_data(:)
    real(dp), intent(out) :: L_B_data(m_data_max), L_kin_data(m_data_max), L_e_data(m_data_max), L_e_2_data(m_data_max), L_sy_data(m_data_max), L_sy2_data(m_data_max)
    integer :: j
    logical, parameter :: verbose = .False.
        
    ! Calculate everything for the jet
    call initialize_allocate(sigma0_value, zmax_value, Lj_value, a_inj_value)
    call trajectory_calc('')  
    call Quantities('')
    call Edist('e')   

    do j = 1, size(z_data)
        L_B_data(j) = L_B(z_data(j))
        L_kin_data(j) = Lp(z_data(j))
        L_e_data(j) = Le_z(z_data(j))
        L_e_2_data(j) = Le_z2(z_data(j))
        L_sy_data(j) = L_syn(z_data(j))
        L_sy2_data(j) = L_syn2(z_data(j))
    end do

    call initialize_deallocate 

    if (verbose) then
        print *, 'z_data[cm] =', z_data
        print *, 'L_B_data[erg/s] =', L_B_data
        print *, 'L_kin_data[erg/s] =', L_kin_data
        print *, 'L_e_data[erg/s] =', L_e_data
    endif

end subroutine power_model

!===================================================================