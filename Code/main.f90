PROGRAM main
    use globals
    use controls
    implicit none
    real(kind(1.d0)), allocatable :: E_data(:)
    real(kind(1.d0)) :: S_data(m_data_max)
    logical, parameter :: verbose = .True.

    interface
        subroutine run_jet(verbose) 
            logical, intent(in) :: verbose
        end subroutine run_jet
    end interface

    interface
        subroutine sed_model(E_vector, S_vector) 
            use globals
            real(kind(1.d0)), intent(in) :: E_vector(:)
            real(kind(1.d0)), intent(out) :: S_vector(m_data_max)
        end subroutine sed_model
    end interface 

    !----------------------------------------------------------------
    ! Original Fortran program

    if (run_jet_calc) call run_jet(verbose)   

    !----------------------------------------------------------------
    ! Short version to calculate specific fluxes, used for fitting SEDs
    
    if (run_sed_calc) then
        allocate( E_data(4) )
        E_data = (/1d-19, 1d-14, 1.d-10, 2.d-8/)
        call sed_model(E_data, S_data) ! Energies in cgs
        deallocate( E_data )
    endif


END PROGRAM main


  !==================================================================

  subroutine run_jet(verbose)
    use controls
    use globals
    use initialize
    use trajectory
    use evolution
    use losses
    use edistribution
    use synrad
    use icrad
!    use fluxes
    
    implicit none
    logical, intent(in) :: verbose
    !real(kind(1.d0)) :: EphSy(NSED_sy), Lsy(NSED_sy), Lsy_CJ(NSED_sy)
    !real(kind(1.d0)) :: EphIC(NSED_IC), LIC(NSED_IC), LIC_CJ(NSED_IC)
    !real(kind(1.d0)) :: Flux(N_bands)
    character(3) :: numNE, numSED

    write(numNE,'(I3)') nIn         !Number of points for the N(E) distribution
    write(numSED,'(I3)') NSED_sy    !Number of points for the synchrotron SED

    ! Calculate everything for the jet first
    if (JET_) then 
        if (verbose) write(*,*) '\n Initialize vectors'    
        call initialize_allocate

        !Jet trajectory
        if (verbose .and. jet_calc) write(*,*) '\n Calculating jet trajectory'
        if (jet_calc) call trajectory_calc('') !'' = jet 
        
        !Evolution of different quantities
        if (verbose .and. jet_calc) write(*,*) '\n Calculating the evolution of quantities along the jet'
        if (jet_calc) call Quantities('') !'' = jet 
        if (verbose) write(*,'(2X,A12,F8.4)') 'B0_FF[G] = ', B0/gj

        !Particle distribution
        if (NE_calc) then 
            if (verbose) write(*,'(1X,A,I3,A)') '\n Computing the particle energy distribution (',nIn,' points)'
            call Edist()        ! Calculate the particle energy distribution
            call execute_command_line('./SortNE.sh '//trim(adjustl(numNE)))   !Sort N(E) file
            call Calc_Ue()      ! Calculate energy densities
            if (verbose) write(*,*) '\n Creating the tloss files for plots'
            call Print_tloss(1)
            call Print_tloss(Ns/3)
            call Print_tloss(Ns)
        endif

        !Synchrotron SED
        if (radsyn_calc) then
            if (verbose) write(*,'(1X,A,I3,A)') '\n Calculating the synchrotron SED (',NSED_sy,' points) in the jet'
            call SEDsyn('')      ! Calculate synchrotron SED
            call execute_command_line('./SortSEDsyn.sh '//trim(adjustl(numSED))) !Sort SED file
            call Calc_Usyn('')   ! Calculate the energy density of synchrotron photons
        endif
        
        !IC SED     
        if (radic_calc) call SEDIC(iso)     ! Calculate the IC-CMB SED
    endif

    ! Calculate everything for the counter-jet next
    if (CJ_) then
        if (.not. JET_) then
            if (verbose) write(*,*) '\n Initialize vectors'    
            call initialize_allocate
        endif

        !Jet trajectory
        if (verbose .and. jet_calc) write(*,*) '\n Calculating counter-jet trajectory'
        if (jet_calc) call trajectory_calc('CJ')
        
        !Evolution of different quantities
        if (verbose .and. jet_calc) write(*,*) '\n Calculating the evolution of quantities along the counter-jet'
        if (jet_calc) call Quantities('CJ')

        !Particle distribution
        if (NE_calc .and. (.not. JET_)) then
            if (verbose) write(*,'(1X,A,I3,A)') '\n Computing the particle energy distribution (',nIn,' points)'
            call Edist()        ! Calculate the particle energy distribution
            call execute_command_line('./SortNE.sh '//trim(adjustl(numNE)))   !Sort N(E) file
            call Calc_Ue()      ! Calculate energy densities
            if (verbose) write(*,*) '\n Creating the tloss files for plots'
            call Print_tloss(1)
            call Print_tloss(Ns/3)
            call Print_tloss(Ns)
        endif

        !Synchrotron SED
        if (radsyn_calc) then
            if (verbose) write(*,'(1X,A,I3,A)') '\n Calculating the synchrotron SED (',NSED_sy,' points) in the counter-jet'
            call SEDsyn('CJ')      ! Calculate synchrotron SED
            call execute_command_line('./SortSEDsyn.sh '//trim(adjustl(numSED))) !Sort SED file
            call Calc_Usyn('CJ')   ! Calculate the energy density of synchrotron photons
        endif
        
        !IC SED     
        if (radic_calc) call SEDIC(iso)     ! Calculate the IC-CMB SED
    endif

    if (verbose) write(*,*) '\n Deallocate vectors'    
    call initialize_deallocate
    
    if (verbose) write(*,*) '\n Program finished'    

  end subroutine run_jet




  !============================================================================
  ! sed_model(START)
  !============================================================================
  subroutine sed_model(E_data, S_data)
  !===================================================================
  ! Subroutine to calculate the synchrotron SED at specific energies
  ! for the given input parameters (zero now).
  !===================================================================
    use controls
    use globals
    use initialize
    use trajectory
    use evolution
    use losses
    use edistribution
    use synrad

    implicit none
    real(kind(1.d0)), intent(in) :: E_data(:)
    real(kind(1.d0)), intent(out) :: S_data(m_data_max)
    real(kind(1.d0)) :: S_data_J(m_data_max), S_data_CJ(m_data_max)
    logical, parameter :: verbose = .True.

    S_data = 0.d0 ! Initialize
    ! Calculate everything for the jet first
    if (JET_) then  
        call initialize_allocate

        !Jet trajectory
        call trajectory_calc('')  ! jet 
        
        !Evolution of different quantities
        call Quantities('')       ! jet 

        !Particle distribution
        call Edist()   

        !Synchrotron SED at the given values
        call SEDsyn_data(E_data, S_data_J, '')   
        
        call initialize_deallocate 
        S_data = S_data_J
    endif 


    ! Calculate everything for the counter-jet next
    if (CJ_) then
        call initialize_allocate

        !Jet trajectory
        call trajectory_calc('CJ')  ! counter-jet 
        
        !Evolution of different quantities
        call Quantities('CJ')       ! counter-jet 

        !Particle distribution
        if (.not. JET_) then
            call Edist()        ! Calculate the particle energy distribution
        endif

        !Synchrotron SED at the given values
        call SEDsyn_data(E_data, S_data_CJ, 'CJ')
        S_data = S_data + S_data_CJ
        
        call initialize_deallocate 
    endif

    if (verbose) then
        write(*,*) 'E_data[eV] =', E_data/eV
        write(*,*) 'S_data[cgs] =', S_data
    endif

  end subroutine sed_model

  !===================================================================
