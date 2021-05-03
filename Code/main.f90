PROGRAM main
    
    use controls
    use globals
    use trajectory
    use evolution
    use losses
    use edistribution
    use synrad
    use icrad
!    use fluxes
    
    implicit none
    real(dp) :: EphSy(NSED_sy),Lsy(NSED_sy),Lsy_CJ(NSED_sy)
    real(dp) :: EphIC(NSED_IC),LIC(NSED_IC),LIC_CJ(NSED_IC)
    real(dp) :: Flux(N_bands)
    character(3) :: gammaj,incli,numNE,numSED

    write(numNE,'(I3)') nNE         !Number of points for the N(E) plot
    write(numSED,'(I3)') NSED_sy    !Number of points for the synchrotron SED
    
    !Jet trajectories
    if (jet_calc) call trajectory_calc('') !'' = jet 
    !if (jet_calc) call trajectory_calc('cj')
    
    !Evolution of different quantities
    if (jet_calc) call Quantities('') !'' = jet 
    !if (jet_calc) call Quantities('cj')


    !Particle distribution
    if (NE_calc) call Edist()
    if (NE_calc) call execute_command_line('./SortNE.sh '//trim(adjustl(numNE)))        !Sort N(E) file
    if ((radsyn_calc .or. radic_calc) .and. NEint_calc) call Edist_int()
    call Print_tloss(1)

    !Synchrotron SED
    if (radsyn_calc) call SEDsyn()
    if (radsyn_calc) call execute_command_line('./SortSEDsyn.sh '//trim(adjustl(numSED))) !Sort SED file
    
    !IC SED     
    if (radic_calc) call SEDIC(iso)
    

END PROGRAM main
