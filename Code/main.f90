PROGRAM main
    use constants
    use globals
    use controls
    implicit none
    real(dp) :: sigma0_value, Lj_value, z0_value, a_inj_value
    real(dp), allocatable :: E_data(:), z_data(:)
    real(dp) :: S_data(m_data_max)
    logical, parameter :: verbose = .False.

    interface
        subroutine run_jet(verbose) 
            logical, intent(in) :: verbose
        end subroutine run_jet
    end interface

    interface
        subroutine sed_model(sigma0_value, Lj_value, z0_value, a_inj_value, E_vector, S_vector) 
            use globals
            use constants
            real(dp), intent(in) :: sigma0_value, Lj_value, z0_value, a_inj_value
            real(dp), intent(in) :: E_vector(:)
            real(dp), intent(out) :: S_vector(m_data_max)
        end subroutine sed_model
    end interface 

    interface
        subroutine syn_model(sigma0_value, Lj_value, z0_value, a_inj_value, E_vector, z_vector, S_vector) 
            use globals
            use constants
            real(dp), intent(in) :: sigma0_value, Lj_value, z0_value, a_inj_value
            real(dp), intent(in) :: E_vector(:)
            real(dp), intent(in) :: z_vector(:)
            real(dp), intent(out) :: S_vector(m_data_max)
        end subroutine syn_model
    end interface 

    !----------------------------------------------------------------
    ! Original Fortran program

    if (run_jet_calc) call run_jet(verbose)   

    !----------------------------------------------------------------
    ! Short version to calculate specific fluxes integrated up to 
    ! certain jet heights. Used for fitting SEDs.
    
    if (run_sed_calc) then
        allocate( E_data(4), z_data(4) )
        E_data = (/1d-6*eV, 1d-2*eV, 1.d2*eV, 1.d5*eV/)
        z_data = (/5d2*R_g, 1d2*R_g, 3d1*R_g, 1.5d1*R_g/)
        sigma0_value = 1.d1
        Lj_value = L_Edd*1d-4
        z0_value = 10.d0*R_g
        a_inj_value = 1.0001d0        
        call syn_model(sigma0_value, Lj_value, z0_value, a_inj_value, E_data, z_data, S_data) ! Energies in cgs
        call sed_model(sigma0_value, Lj_value, z0_value, a_inj_value, E_data, S_data) ! Energies in cgs
        deallocate( E_data, z_data )
    endif


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
!    use fluxes
    
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
        call initialize_allocate(sigma_0, Lj, z0, a_inj)
        !Jet trajectory
        if (verbose .and. jet_calc) print *, '\n Calculating jet trajectory'
        if (jet_calc) call trajectory_calc('') !'' = jet 
        
        !Evolution of different quantities
        if (verbose .and. jet_calc) print *, '\n Calculating the evolution of quantities along the jet'
        if (jet_calc) call Quantities('') !'' = jet 
        !if (verbose) write(*,'(2X,A12,F8.4)') 'B0_FF[G] = ', B0/gj

        !Particle distribution
        if (NE_calc) then 
            if (verbose) write(*,'(1X,A,I3,A)') '\n Computing the particle energy distribution (',nIn,' points)'
            call Edist('e')        ! Calculate the particle energy distribution
            call execute_command_line('./SortNe.sh '//trim(adjustl(numNE)))   !Sort N(E) file
            call Calc_Ue()      ! Calculate energy densities
            if (verbose) print *, '\n Creating the tloss files for plots'
            call Print_tloss(1, 'e')
            call Print_tloss(Ns/4, 'e')
            call Print_tloss(Ns/2, 'e')
            call Print_tloss(Ns, 'e')
        endif

        if (Np_calc) then
            call Edist('p')
            call execute_command_line('./SortNp.sh '//trim(adjustl(numNE)))   !Sort N(E) file
        endif

        !Synchrotron SED
        if (radsyn_calc) then
            if (verbose) write(*,'(1X,A,I3,A)') '\n Calculating the synchrotron SED (',NSED_sy,' points) in the jet'
            call SEDsyn('', 'e')      ! Calculate synchrotron SED
            call execute_command_line('./SortSEDesyn.sh '//trim(adjustl(numSED))) !Sort SED file
            call Calc_Usyn('')   ! Calculate the energy density of synchrotron photons
        endif

        if (radpsyn_calc) then
            call SEDsyn('', 'p')
            call execute_command_line('./SortSEDpsyn.sh '//trim(adjustl(numSED))) !Sort SED file
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
            call initialize_allocate(sigma_0,Lj,z0,a_inj)
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
            call execute_command_line('./SortNe.sh '//trim(adjustl(numNE)))   !Sort N(E) file
            call Calc_Ue()      ! Calculate energy densities
            if (verbose) print *, '\n Creating the tloss files for plots'
            call Print_tloss(1, 'e')
            call Print_tloss(Ns/3, 'e')
            call Print_tloss(Ns, 'e')
        endif

        !Synchrotron SED
        if (radsyn_calc) then
            if (verbose) write(*,'(1X,A,I3,A)') '\n Calculating the synchrotron SED (',NSED_sy,' points) in the counter-jet'
            call SEDsyn('CJ', 'e')      ! Calculate synchrotron SED
            call execute_command_line('./SortSEDesyn.sh '//trim(adjustl(numSED))) !Sort SED file
            call Calc_Usyn('CJ')   ! Calculate the energy density of synchrotron photons
        endif

        if (radpsyn_calc) then
            call SEDsyn('', 'p')
            call execute_command_line('./SortSEDpsyn.sh '//trim(adjustl(numSED))) !Sort SED file
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




  !============================================================================
  ! sed_model(START)
  !============================================================================
  subroutine sed_model(sigma0_value, Lj_value, z0_value, a_inj_value, E_data, S_data)
  !===================================================================
  ! Subroutine to calculate the synchrotron + SSC SED at specific 
  ! energies for the given input parameters. Calculates both J and CJ
  ! TO DO: INCORPORATE SSC IN THE COMPUTATION!!!
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
    real(dp), intent(in) :: sigma0_value, Lj_value, z0_value, a_inj_value
    real(dp), intent(in) :: E_data(:)
    real(dp), intent(out) :: S_data(m_data_max)
    real(dp) :: S_data_J(m_data_max), S_data_CJ(m_data_max)
    logical, parameter :: verbose = .False.
        
    S_data = 0.d0 ! Initialize
    ! Calculate everything for the jet first
    if (JET_) then  
        call initialize_allocate(sigma0_value, Lj_value, z0_value, a_inj_value)

        !Jet trajectory
        call trajectory_calc('')  ! jet 
        
        !Evolution of different quantities
        call Quantities('')       ! jet 

        !Particle distribution
        call Edist('e')   

        !Synchrotron SED at the given values  
        call SEDsyn_data(E_data, S_data_J, '', 'e')  
        S_data = S_data_J

        !SSC SED at the given values  
        !call SEDssc_data(E_data, S_data_J, '', 'e') 
        !S_data = S_data + S_data_J 
        if (verbose) print *, 'S_data[cgs] =', S_data

        call initialize_deallocate 
    endif 


    ! Calculate everything for the counter-jet next
    if (CJ_) then
        call initialize_allocate(sigma0_value, Lj_value, z0_value, a_inj_value)

        !Jet trajectory
        call trajectory_calc('CJ')  ! counter-jet 
        
        !Evolution of different quantities
        call Quantities('CJ')       ! counter-jet 

        !Particle distribution
        if (.not. JET_) call Edist('e')  ! Calculate the particle energy distribution

        !Synchrotron SED at the given values
        call SEDsyn_data(E_data, S_data_CJ, 'CJ', 'e')
        S_data = S_data + S_data_CJ

        !SSC SED at the given values  
        !call SEDssc_data(E_data, S_data_CJ, '', 'e') 
        !S_data = S_data + S_data_CJ  
        if (verbose) print *, 'S_data[cgs] =', S_data
        
        call initialize_deallocate 
    endif

    if (verbose) then
        print *, 'E_data[eV] =', E_data/eV
        print *, 'S_data[cgs] =', S_data
    endif

  end subroutine sed_model

  !===================================================================



  !============================================================================
  ! syn_model(START)
  !============================================================================
  subroutine syn_model(sigma0_value, Lj_value, z0_value, a_inj_value, E_data, z_data, S_data)
  !===================================================================
  ! Subroutine to calculate the synchrotron SED at specific energies
  ! integrated up to certain heights for the given input parameters. 
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
    real(dp), intent(in) :: sigma0_value, Lj_value, z0_value, a_inj_value
    real(dp), intent(in) :: E_data(:), z_data(:)
    real(dp), intent(out) :: S_data(m_data_max)
    real(dp) :: S_data_J(m_data_max), S_data_CJ(m_data_max)
    logical, parameter :: verbose = .False.
        
    ! Calculate everything for the jet first
    call initialize_allocate(sigma0_value, Lj_value, z0_value, a_inj_value)
    call trajectory_calc('')  
    call Quantities('')    
    call Edist('e')   
    !Synchrotron SED at the given values
    call SEDsyn_z_data(E_data, z_data, S_data_J, '', 'e')   
    call SEDsyn_data(E_data, S_data_J, '', 'e')  

    ! Calculate everything for the counter-jet next
    call trajectory_calc('CJ')  
    call Quantities('CJ')       
    call SEDsyn_z_data(E_data, z_data, S_data_CJ, 'CJ', 'e')
    call initialize_deallocate 

    S_data = S_data_J + S_data_CJ

    if (verbose) then
        print *, 'E_data[eV] =', E_data/eV
        print *, 'z_data[Rg] =', z_data/R_g
        print *, 'S_data[cgs] =', S_data
    endif

  end subroutine syn_model

  !===================================================================
