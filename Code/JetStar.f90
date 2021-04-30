PROGRAM JetStar
	use controls
	use globals
	use trajectory
	use evolution
	use losses
	use edistribution
	use synrad
	use icrad
	use opac_ff
	use opac_gg
	use occultation
	use fluxes
	
	implicit none
	real*8 :: EphSy(NSED_sy),Lsy(NSED_sy),Lsy_CJ(NSED_sy)
	real*8 :: EphIC(NSED_IC),LIC(NSED_IC),LIC_CJ(NSED_IC)
	real*8 :: Flux(N_bands)
	integer :: Norb,iorb
	character(3) :: gammaj,incli,numNE,numSED
	
	!Transform numbers to strings
	write(gammaj,'(F3.1)') gj		!Jet Lorentz factor
	write(incli,'(I3)') int(inc)	!Inclination
	write(numNE,'(I3)') nNE			!Number of points for the N(E) plot
	write(numSED,'(I3)') NSED_sy	!Number of points for the synchrotron SED
	
	orbPhase = 0d0		!Initial orbital phase (deg)
	Norb = 1			!Number of orbital phases (0 for only trajectory and N(E) computation)

	!Jet trajectories
	if (curvedJet) call JetPath()
	if (curvedJet) call JetPathCJ()
	if (.not. curvedJet) call JetPath_nocurv()		!Vertical jets
	if (.not. curvedJet) call JetPathCJ_nocurv()	!
	
	!Evolution of different quantities
	call Quantities()
	call QuantitiesCJ()
	
	!Particle distribution
	if (NE_calc) call Edist()
	if (NE_calc) call execute_command_line('./SortNE.sh '//trim(adjustl(numNE)))		!Sort N(E) file
	if ((radsyn_calc .or. radic_calc) .and. NEint_calc) call Edist_int()
	call Print_tloss(1)
	
	if (curvedJet) then
		open(100,file='Gamma'//gammaj//'/SEDsyn_orb_i'//trim(adjustl(incli))//'.dat')
		open(110,file='Gamma'//gammaj//'/SEDIC_orb_i'//trim(adjustl(incli))//'.dat')
		open(120,file='Gamma'//gammaj//'/Fluxes_i'//trim(adjustl(incli))//'.dat')
	else
		open(100,file='Gamma'//gammaj//'/SEDsyn_orb_i'//trim(adjustl(incli))//'_nocurv.dat')
		open(110,file='Gamma'//gammaj//'/SEDIC_orb_i'//trim(adjustl(incli))//'_nocurv.dat')
		open(120,file='Gamma'//gammaj//'/Fluxes_i'//trim(adjustl(incli))//'_nocurv.dat')
	endif
	
	do iorb = 1,Norb
	
		write(*,'(X,A,F4.2,A)'), '-------------------- Phase = ',orbPhase/360d0,' --------------------'
	
		!Jet trajectories
		if (curvedJet) call JetPath()
		if (curvedJet) call JetPathCJ()
		if (.not. curvedJet) call JetPath_nocurv()
		if (.not. curvedJet) call JetPathCJ_nocurv()
		call Quantities()
		call QuantitiesCJ()
	
		!Synchrotron SED		
		if (radsyn_calc) then
			call SEDsyn()
			call SEDsynCJ()
			call execute_command_line('./SortSEDsyn.sh '//trim(adjustl(numSED)))	!Sort SED file
			call opac_ff_calc()
			call opac_ff_calcCJ()
			call Occult_Sy(EphSy,Lsy)
			call Occult_SyCJ(Lsy_CJ)
			do i = 1,NSED_sy
				write(100,*) orbPhase,EphSy(i)/eV,Lsy(i),Lsy_CJ(i),Lsy(i)+Lsy_CJ(i)
			enddo 
			write(100,*)	!Blank lines
		endif
	
		!IC SED		
		if (radic_calc) then
			call SEDIC()
			call SEDIC_CJ()
			call opac_gg_calc()
			call opac_gg_calcCJ()
			call Occult_IC(EphIC,LIC)
			call Occult_ICCJ(LIC_CJ)
			do i = 1,NSED_IC
				write(110,*) orbPhase,EphIC(i)/eV,LIC(i),LIC_CJ(i),LIC(i)+LIC_CJ(i)
			enddo 
			write(110,*)	!Blank lines
		endif
		
		if (radsyn_calc .and. radic_calc) call TotalEmission()	!Jet + counter-jet added SEDs		
		if (Norb == 1) exit		!Exit the loop when one only wants to compute the SEDs for a given phase
		
		!Fluxes
		if (radsyn_calc .and. radic_calc) call fluxes_calc(EphSy,Lsy,Lsy_CJ,EphIC,LIC,LIC_CJ,Flux)		
		if (radsyn_calc .and. radic_calc) write(120,*) orbPhase,(flux(j), j=1,N_bands)
		
		!Increase in the orbital phase
		if (curvedJet) then
			orbPhase = orbPhase + 360d0/(Norb-1)	
		else
			orbPhase = orbPhase + 180d0/(Norb-1)	!Only compute half orbit for vertical jets (symmetric)
		endif
		
	enddo
	
	close(100)
	close(110)
	close(120)
	
END PROGRAM JetStar
