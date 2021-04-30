MODULE synrad
	use globals
	implicit none

contains

!***************************************************************************************************
!	Computation of the synchrotron SED (jet)
!***************************************************************************************************
	subroutine SEDsyn ()
	real*8 :: E(nIn),NE(nIn,Ns),NEcut(nIn),Ein,Eph,dEph,EphFF(Ns),dEphFF(Ns),Lumcut(Ns),LumcutFF(Ns),	&
	&		  L1,L2,Lum,LumFF,LBolFF,LumcutA(Ns),LumcutAFF(Ns),LumA,LumAFF,LBolAFF
	real*8 :: dummy
	allocate (lenj(Ns),BFF(Ns),rj(Ns),sjoFF(Ns),VFF(Ns),delta_obs(Ns))
	
	!Quantities at each segment
	open(10,file='Quant_FF.dat')
	open(20,file='Deltas.dat')
	do i = 1,Ns
		read(10,*) lenj(i),BFF(i),dummy,rj(i),dummy,dummy,dummy,dummy,dummy,sjoFF(i),VFF(i)
		read(20,*) dummy,delta_obs(i)
	enddo
	close(10)
	close(20)
	
	!N(E) of each cut
	open(30,file='NE_intcuts.dat')
   	do i = 1,nIn
		do j = 1,Ns
			read(30,*) E(i),NE(i,j)
		enddo
	enddo
	close(30)

	open(40,file='SEDsyn_cuts.dat')
	open(50,file='SEDsyn_jet.dat')

	write(*,'(1X,A,I3,A)') 'Synchrotron SED (',NSED_sy,')'
	
	Ein = (Emax_sy/Emin_sy)**(1d0/(NSED_sy-1))	!Energy increase factor
	Eph = Emin_sy								!Initialization

	!Initialization of the bolometric luminosities
	LBolFF = 0d0
	LBolAFF = 0d0

	do l = 1,NSED_sy

		!Initialization of the total luminosity for each photon energy
	 	Lum = 0d0
	 	LumFF = 0d0
	 	LumA = 0d0
	 	LumAFF = 0d0
	 	
	 	dEph = Eph*(Ein-1d0)	!Energy interval

		!Loop for the segments
		do j = 1,Ns
		
			EphFF(j) = Eph/delta_obs(j)		!Photon energy in the FF
			dEphFF(j) = dEph/delta_obs(j)	!Energy interval in the FF
		
			do i = 1,nIn
				!N(E) for the segment j
	   			NEcut(i) = NE(i,j)
	  		enddo

			!Unabsorbed synchrotron luminosity for each segment
			LumcutFF(j) = SyncLum(BFF(j),EphFF(j),E,NEcut)
			Lumcut(j) = delta_obs(j)**2/gj*LumcutFF(j)		!(de)Boosted

			!Unabsorbed (de)boosted luminosity for Eph
			Lum = Lum + Lumcut(j)
			
			LumFF = LumFF + LumcutFF(j)*dEphFF(j)

	 	enddo
	 	
	 	!Unabsorbed bolometric luminosity
	 	LBolFF = LBolFF + LumFF
	 		 	
	 	!Absorbed radiation for each cut
	 	call SSA(lenj,BFF,rj,sjoFF,VFF,EphFF,E,NEcut,LumcutFF,LumcutAFF)
	 	
	 	do j = 1,Ns
	 	
	 		LumcutA(j) = delta_obs(j)**2/gj*LumcutAFF(j)	!(de)Boosted
	 		write(40,*) Eph/eV,Eph*Lumcut(j),Eph*LumcutA(j),j
	 		
	 		!Absorbed (de)boosted luminosity for Eph
			LumA = LumA + LumcutA(j)
			
			LumAFF = LumAFF + LumcutAFF(j)*dEphFF(j)
			
			!Jet and c-jet contribution of the segments up to given lengths
	  		if (lenj(j)-lenj(1) <= 0.1d0*a) L1 = LumA
	  		if (lenj(j)-lenj(1) <= a) L2 = LumA
			
	 	enddo

	 	write(50,*) Eph/eV,Eph*L1,Eph*L2,Eph*Lum,Eph*LumA
	 	
	 	!Absorbed bolometric luminosity
	 	LBolAFF = LBolAFF + LumAFF
	 	
	 	Eph = Eph*Ein		!Increase in the photon energy

	enddo
	
	!write(*,'(3X,A,ES13.6E2,1X,A)') 'Jet unabsorbed synchrotron luminosity (FF):',LBolFF,'erg/s'
	!write(*,'(3X,A,ES13.6E2,1X,A)') 'Jet absorbed synchrotron luminosity (FF):',LBolAFF,'erg/s'

	close(40)
	close(50)
	
	deallocate (lenj,BFF,rj,sjoFF,VFF,delta_obs)
	return
	end subroutine SEDsyn
	
!*************************************************************************************************************
!	Computation of the synchrotron SED (counter-jet)
!*************************************************************************************************************
	subroutine SEDsynCJ ()
	real*8 :: E(nIn),NE(nIn,Ns),NEcut(nIn),Ein,Eph,dEph,EphFF(Ns),dEphFF(Ns),Lumcut(Ns),LumcutFF(Ns),	&
	&		  L1,L2,Lum,LumFF,LBolFF,LumcutA(Ns),LumcutAFF(Ns),LumA,LumAFF,LBolAFF
	real*8 :: dummy
	allocate (lenj(Ns),BFF(Ns),rj(Ns),sjoFF(Ns),VFF(Ns),delta_obs(Ns))
	
	!Quantities at each segment
	open(10,file='Quant_FFCJ.dat')
	open(20,file='DeltasCJ.dat')
	do i = 1,Ns
		read(10,*) lenj(i),BFF(i),dummy,rj(i),dummy,dummy,dummy,dummy,dummy,sjoFF(i),VFF(i)
		read(20,*) dummy,delta_obs(i)
	enddo
	close(10)
	close(20)
	
	!N(E) of each cut
	open(30,file='NE_intcuts.dat')
   	do i = 1,nIn
		do j = 1,Ns
			read(30,*) E(i),NE(i,j)
		enddo
	enddo
	close(30)

	open(40,file='SEDsyn_cutsCJ.dat')
	open(50,file='SEDsyn_CJ.dat')
	
	Ein = (Emax_sy/Emin_sy)**(1d0/(NSED_sy-1))	!Energy increase factor
	Eph = Emin_sy								!Initialization

	!Initialization of the bolometric luminosities
	LBolFF = 0d0
	LBolAFF = 0d0

	do l = 1,NSED_sy

		!Initialization of the total luminosity for each photon energy
	 	Lum = 0d0
	 	LumFF = 0d0
	 	LumA = 0d0
	 	LumAFF = 0d0
	 	
	 	dEph = Eph*(Ein-1d0)	!Energy interval

		!Loop for the segments
		do j = 1,Ns		
		
			EphFF(j) = Eph/delta_obs(j)		!Photon energy in the FF
			dEphFF(j) = dEph/delta_obs(j)	!Energy interval in the FF
		
			do i = 1,nIn
				!N(E) for the segment j
	   			NEcut(i) = NE(i,j)
	  		enddo

			!Unabsorbed synchrotron luminosity for each segment
			LumcutFF(j) = SyncLum(BFF(j),EphFF(j),E,NEcut)
			Lumcut(j) = delta_obs(j)**2/gj*LumcutFF(j)		!(de)Boosted

			!Unabsorbed (de)boosted luminosity for Eph
			Lum = Lum + Lumcut(j)
			
			LumFF = LumFF + LumcutFF(j)*dEphFF(j)

	 	enddo
	 	
	 	!Unabsorbed bolometric luminosity
	 	LBolFF = LBolFF + LumFF
	 		 	
	 	!Absorbed radiation for each cut
	 	call SSA(lenj,BFF,rj,sjoFF,VFF,EphFF,E,NEcut,LumcutFF,LumcutAFF)
	 	
	 	do j = 1,Ns
	 	
	 		LumcutA(j) = delta_obs(j)**2/gj*LumcutAFF(j)	!(de)Boosted
	 		write(40,*) Eph/eV,Eph*Lumcut(j),Eph*LumcutA(j),j
	 		
	 		!Absorbed (de)boosted luminosity for Eph
			LumA = LumA + LumcutA(j)
			
			LumAFF = LumAFF + LumcutAFF(j)*dEphFF(j)
			
			!Jet and c-jet contribution of the segments up to given lengths
	  		if (lenj(j)-lenj(1) <= 0.1d0*a) L1 = LumA
	  		if (lenj(j)-lenj(1) <= a) L2 = LumA
			
	 	enddo

	 	write(50,*) Eph/eV,Eph*L1,Eph*L2,Eph*Lum,Eph*LumA
	 	
	 	!Absorbed bolometric luminosity
	 	LBolAFF = LBolAFF + LumAFF
	 	
	 	Eph = Eph*Ein		!Increase in the photon energy

	enddo
	
	!write(*,'(3X,A,ES13.6E2,1X,A)') 'C-jet unabsorbed synchrotron luminosity (FF):',LBolFF,'erg/s'
	!write(*,'(3X,A,ES13.6E2,1X,A)') 'C-jet SSAbsorbed synchrotron luminosity (FF):',LBolAFF,'erg/s'

	close(40)
	close(50)
	
	deallocate (lenj,BFF,rj,sjoFF,VFF,delta_obs)
	return
	end subroutine SEDsynCJ

!*************************************************************************************************************
!	Synchrotron luminosity per photon energy unit
!*************************************************************************************************************
	real*8 function SyncLum (B,Eph,E,NE)
	real*8, intent(in) :: B,Eph,E(nIn),NE(nIn)
	real*8 :: dE(nIn),Ec,x,Psync

	!Initialization
	SyncLum = 0d0

	!Integral by rectangles
	do i = 1,nIn

	 	if (i /= nIn) dE(i) = E(i+1) - E(i)		!Energy interval

		!Critical photon energy and power emitted for each particle energy
	 	Ec = dsqrt(3d0/2d0) * (qe*h*B*E(i)**2)/(2d0*pi*me**3*c**5)
	 	x = Eph/Ec
	 	Psync = 1.85d0*dsqrt(2d0)*qe**3*B/(h*me*c**2) * x**(1d0/3d0)*dexp(-x)

		!The approximation is only valid for 0.1 < x < 10
		if ((x .lt. 1d-1) .or. (x .gt. 1d1)) Psync = 0d0

		SyncLum = SyncLum + Psync*NE(i)*dE(i)

	enddo

	return
	end function SyncLum
	
!*************************************************************************************************************
!	Luminosity corrected for synchrotron self-absorption
!*************************************************************************************************************
	subroutine SSA (lenj,B,rj,sjo,V,Eph,E,NE,Lum,LumA)
	real*8, intent(in) :: lenj(Ns),B(Ns),rj(Ns),sjo(Ns),V(Ns),Eph(Ns),E(nIn),NE(nIn),Lum(Ns)
	real*8, intent(out) :: LumA(Ns)
	real*8 :: K,Ecf,Sf,depth,tau,Ap

	!Loop for the segments
	do j = 1,Ns

		!Absorption and emission coefficients, and source function, for each segment (Pacholczyck, 1970)
	 	K = 2.62d0*qe**3*h**2*B(j)/(4d0*pi*Eph(j)**2*me*V(j)) * EIntSSA(B(j),Eph(j),E,NE)
	 	Ecf = Lum(j)/(4d0*pi*V(j))		!Per stereoradian
	 	Sf = Ecf/K

		!Length that photons in the observer's direction must travel inside the segment (in the FF)
	 	depth = rj(j)/sjo(j)

		!Optical depth for SSA
	 	tau = K*min(depth,lenj(j)/gj)

		!Estimation of the projected area, perpendicular to the observer's direction, obtained by equating
		!the unabsorbed Lum to its value for optically thin emission: 4*pi*Ap*Sf*tau
		Ap = V(j)/min(depth,lenj(j)/gj)

		!Luminosity corrected for synchrotron self-absorption
	 	if (tau .le. 0.1d0) then
	  		LumA(j) = Lum(j)
	 	else
	  		LumA(j) = 4d0*pi*Ap*Sf*(1d0-dexp(-tau))
	 	endif

	enddo

	return
	end subroutine SSA

!*************************************************************************************************************
!	Integral for the particle energy in SSA
!*************************************************************************************************************
	real*8 function EintSSA (B,Eph,E,NE)
	real*8, intent(in) :: B,Eph,E(nIn),NE(nIn)
	real*8 :: dE(nIn),Ec,x,F
	
	EIntSSA = 0d0	!Initialization
	
	do i = 1,nIn

	 	if (i /= nIn) dE(i) = E(i+1) - E(i)		!Energy interval

		!Critical photon energy and function to integrate for each particle energy
	 	Ec = dsqrt(3d0/2d0) * (qe*h*B*E(i)**2)/(2d0*pi*me**3*c**5)
	 	x = Eph/Ec
	 	F = NE(i)/E(i)*dexp(-x)*x**(1d0/3d0)*(x + 2d0/3d0)

		!The approximation is only valid for 0.1 < x < 10
		if (x < 1d-1 .or. x > 1d1) F = 0d0

		EintSSA = EintSSA + F*dE(i)

	enddo

	return
	end function EintSSA
	
END MODULE synrad
