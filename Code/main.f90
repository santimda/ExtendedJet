PROGRAM main
    
    implicit none

    interface
    subroutine jet_calc(verbose) 
      logical, intent(in) :: verbose
    end subroutine jet_calc
    end interface

    call run_jet(.False.)    

END PROGRAM main


  !==================================================================

  subroutine run_jet(verbose)
    use controls
    use globals
    use trajectory
    use evolution
    use losses
    use edistribution
    use synrad
    use icrad
    use map
!    use fluxes
    
    implicit none
    logical, intent(in) :: verbose
    real(dp) :: EphSy(NSED_sy), Lsy(NSED_sy), Lsy_CJ(NSED_sy)
    real(dp) :: EphIC(NSED_IC), LIC(NSED_IC), LIC_CJ(NSED_IC)
    real(dp) :: Flux(N_bands)
    character(3) :: gammaj,incli,numNE,numSED

    write(numNE,'(I3)') nIn         !Number of points for the N(E) distribution
    write(numSED,'(I3)') NSED_sy    !Number of points for the synchrotron SED
    
    !Jet trajectories
    if (verbose .and. jet_calc) write(*,*) '\n Calculating jet trajectory'
    if (jet_calc) call trajectory_calc('') !'' = jet 
    !if (verbose .and. jet_calc) write(*,*) '\n Calculating counter-jet trajectory'
    !if (jet_calc) call trajectory_calc('cj')
    
    !Evolution of different quantities
    if (verbose .and. jet_calc) write(*,*) '\n Calculating the evolution of quantities along the jet'
    if (jet_calc) call Quantities('') !'' = jet 
    if (verbose) write(*,'(2X,A12,F8.4)') 'B0_FF[G] = ', B0/gj

    !if (jet_calc) call Quantities('cj')


    !Particle distribution
    if (NE_calc) then 
        if (verbose) write(*,'(1X,A,I3,A)') '\n Computing the particle energy distribution (',nIn,' points)'
        call Edist()        ! Calculate the particle energy distribution
        call execute_command_line('./SortNE.sh '//trim(adjustl(numNE)))   !Sort N(E) file
        call Calc_Ue()      ! Calculate energy densities
        call Print_tloss(1)
    endif

    !Synchrotron SED
    if (radsyn_calc) then
        if (verbose) write(*,'(1X,A,I3,A)') '\n Calculating the synchrotron SED (',NSED_sy,' points)'
        call SEDsyn()      ! Calculate synchrotron SED
        call execute_command_line('./SortSEDsyn.sh '//trim(adjustl(numSED))) !Sort SED file
        call Calc_Usyn()   ! Calculate the energy density of synchrotron photons
    endif
    
    !IC SED     
    if (radic_calc) call SEDIC(iso)     ! Calculate the IC-CMB SED

    !Emission map
    if (map_calc) call map_sy(1.4d0) ! Frequency in GHz at which to calculate the map
    !if (map_calc) call map_sy_conv() ! Frequency in GHz at which to calculate the map
    
  end subroutine run_jet


  !==================================================================

