MODULE occultation
	use globals
	implicit none
	
contains

!*************************************************************************************************************
!	Correction of the synchrotron emission by occultation by the star (jet)
!*************************************************************************************************************
	subroutine Occult_Sy (Eph,Lsy_tot_SSA_FF_occ)
	real*8, intent(out) :: Eph(NSED_sy),Lsy_tot_SSA_FF_occ(NSED_sy)
	real*8 :: Lsy_noabs(NSED_sy),Lsy_SSA(NSED_sy),Lsy_SSA_FF(NSED_sy),Lsy_SSA_FF_occ(NSED_sy)
	real*8 :: Lsy1(NSED_sy),Lsy2(NSED_sy),Lsy_tot_noabs(NSED_sy),Lsy_tot_SSA(NSED_sy),	&
	&		  Lsy_tot_SSA_FF(NSED_sy)
	real*8 :: dummy,dE,LBol,LBolA
	allocate (isOccult(Ns),lenj(Ns))
	
	open(10,file='Trajectory.dat')
	open(20,file='Quant_FF.dat')
	do i = 1,Ns
		read(10,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,isOccult(i)
		read(20,*) lenj(i)
	enddo
	close(10)
	close(20)
	
	!Initialization of the total luminosities
	Lsy_tot_noabs = 0d0
    Lsy_tot_SSA = 0d0
    Lsy_tot_SSA_FF = 0d0
    Lsy_tot_SSA_FF_occ = 0d0

    open(40,file='SEDsyn_abs_cuts.dat')
    open(41,file='SEDsyn_final_cuts.dat')
    
    !Segment loop
	do j = 1,Ns

    	do i = 1,NSED_sy
			read(40,*) Eph(i),Lsy_noabs(i),Lsy_SSA(i),Lsy_SSA_FF(i)
		enddo
		Eph = Eph*eV			!Conversion to erg       	       
       	if (j < Ns) read(40,*)	!Blank space in the text file
       	
       	!Luminosity is null if occulted by the star
       	if (isOccult(j)) then
       		Lsy_SSA_FF_occ = 0d0
       	else
       	 	Lsy_SSA_FF_occ = Lsy_SSA_FF
       	endif
       	
       	!Total luminosities
       	do i = 1,NSED_sy      		
       		
       		Lsy_tot_noabs(i) = Lsy_tot_noabs(i) + Lsy_noabs(i)					!No absorption
       		Lsy_tot_SSA(i) = Lsy_tot_SSA(i) + Lsy_SSA(i)						!SSA
       		Lsy_tot_SSA_FF(i) = Lsy_tot_SSA_FF(i) + Lsy_SSA_FF(i)				!SSA + FFA
       		Lsy_tot_SSA_FF_occ(i) = Lsy_tot_SSA_FF_occ(i) + Lsy_SSA_FF_occ(i)	!SSA + FFA + occultation
       		
    		!Jet contribution of the segments up to given lengths
    		if (lenj(j)-lenj(1) <= 0.1d0*a) Lsy1(i) = Lsy_tot_SSA_FF_occ(i)
    		if (lenj(j)-lenj(1) <= a) Lsy2(i) = Lsy_tot_SSA_FF_occ(i)
       	  
			write(41,*) Eph(i)/eV,Lsy_noabs(i),Lsy_SSA(i),Lsy_SSA_FF(i),Lsy_SSA_FF_occ(i),j
       	  
       	enddo
       	write(41,*)	!Blank space      	
       	
	enddo
	
	close(40)
	close(41)
	
	LBol = 0d0	
    LBolA = 0d0
    open(42,file='SEDsyn_final_jet.dat')
    do i = 1,NSED_sy
    	write(42,*) Eph(i)/eV,Lsy1(i),Lsy2(i),Lsy_tot_noabs(i),Lsy_tot_SSA(i),Lsy_tot_SSA_FF(i),	&
    	&			Lsy_tot_SSA_FF_occ(i)
    	if (i /= NSED_sy) dE = Eph(i+1)-Eph(i)
    	LBol = LBol + Lsy_tot_noabs(i)/Eph(i)*dE     	  
    	LBolA = LBolA + Lsy_tot_SSA_FF_occ(i)/Eph(i)*dE		!Bolometric synchrotron luminosity
    enddo
    close(42)
    
    write(*,'(3X,A,ES11.4E2,1X,A)') 'Jet unabsorbed synchrotron luminosity (obs):',LBol,'erg/s'
    write(*,'(3X,A,ES11.4E2,1X,A)') 'Jet absorbed synchrotron luminosity (obs):',LBolA,'erg/s'
	
	deallocate (isOccult,lenj)
	return
	end subroutine Occult_Sy
	
!*************************************************************************************************************
!	Correction of the IC emission by occultation by the star (jet)
!*************************************************************************************************************
	subroutine Occult_IC (Eph,LIC_tot_GG_occ)
	real*8, intent(out) :: Eph(NSED_IC),LIC_tot_GG_occ(NSED_IC)
	real*8 :: LIC_noabs(NSED_IC),LIC_GG(NSED_IC),LIC_GG_occ(NSED_IC)
	real*8 :: LIC1(NSED_IC),LIC2(NSED_IC),LIC_tot_noabs(NSED_IC),LIC_tot_GG(NSED_IC)
	real*8 :: dummy,dE,LBol,LBolA
	allocate (isOccult(Ns),lenj(Ns))
	
	open(10,file='Trajectory.dat')
	open(20,file='Quant_FF.dat')
	do i = 1,Ns
		read(10,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,isOccult(i)
		read(20,*) lenj(i)
	enddo
	close(10)
	close(20)
	
	!Initialization of the total luminosities
    LIC_tot_noabs = 0d0
    LIC_tot_GG = 0d0
    LIC_tot_GG_occ = 0d0

    open(40,file='SEDIC_abs_cuts.dat')
    open(41,file='SEDIC_final_cuts.dat')
    
    !Segment loop
	do j = 1,Ns

        do i = 1,NSED_IC
        	read(40,*) Eph(i),LIC_noabs(i),LIC_GG(i)
       	enddo       	
       	Eph = Eph*eV			!Conversion to erg
       	if (j < Ns) read(40,*)	!Blank space in the text file
       	
       	!Luminosities are null if occulted by the star
       	if (isOccult(j)) then
       		LIC_GG_occ = 0d0
       	else
       		LIC_GG_occ = LIC_GG
       	endif
       	
       	!Total luminosities
       	do i = 1,NSED_IC      		
       		
       		LIC_tot_noabs(i) = LIC_tot_noabs(i) + LIC_noabs(i)			!No absorption
       		LIC_tot_GG(i) = LIC_tot_GG(i) + LIC_GG(i)					!GG
       		LIC_tot_GG_occ(i) = LIC_tot_GG_occ(i) + LIC_GG_occ(i)		!GG + occultation
       		
       	 	!Jet contribution of the segments up to given lengths
			if (lenj(j)-lenj(1) <= 0.1d0*a) LIC1(i) = LIC_tot_GG_occ(i)
			if (lenj(j)-lenj(1) <= a) LIC2(i) = LIC_tot_GG_occ(i)
       	  
			write(41,*) Eph(i)/eV,LIC_noabs(i),LIC_GG(i),LIC_GG_occ(i),j
       	  
       	enddo
       	write(41,*)	!Blank space       	
       	
	enddo
	
	close(40)
	close(41)
    
    LBol = 0d0
    LBolA = 0d0
    open(42,file='SEDIC_final_jet.dat')
    do i = 1,NSED_IC
    	write(42,*) Eph(i)/eV,LIC1(i),LIC2(i),LIC_tot_noabs(i),LIC_tot_GG(i),LIC_tot_GG_occ(i)
		if (i /= NSED_IC) dE = Eph(i+1)-Eph(i)
		LBol = LBol + LIC_tot_noabs(i)/Eph(i)*dE        	  
		LBolA = LBolA + LIC_tot_GG_occ(i)/Eph(i)*dE		!Bolometric IC luminosity
    enddo
    close(42)

    write(*,'(3X,A,ES9.2E2,1X,A)') 'Jet unabsorbed IC luminosity (obs):',LBol,'erg/s'
    write(*,'(3X,A,ES9.2E2,1X,A)') 'Jet absorbed IC luminosity (obs):',LBolA,'erg/s'	
	
	deallocate (isOccult,lenj)
	return
	end subroutine Occult_IC
	
!*************************************************************************************************************
!	Correction of the synchrotron emission by occultation by the star (counter-jet)
!*************************************************************************************************************
	subroutine Occult_SyCJ (Lsy_tot_SSA_FF_occ)
	real*8, intent(out) :: Lsy_tot_SSA_FF_occ(NSED_sy)
	real*8 :: Eph(NSED_sy),Lsy_noabs(NSED_sy),Lsy_SSA(NSED_sy),Lsy_SSA_FF(NSED_sy),Lsy_SSA_FF_occ(NSED_sy)
	real*8 :: Lsy1(NSED_sy),Lsy2(NSED_sy),Lsy_tot_noabs(NSED_sy),Lsy_tot_SSA(NSED_sy),	&
	&		  Lsy_tot_SSA_FF(NSED_sy)
	real*8 :: dummy,dE,LBol,LBolA
	allocate (isOccult(Ns),lenj(Ns))
	
	open(10,file='TrajectoryCJ.dat')
	open(20,file='Quant_FFCJ.dat')
	do i = 1,Ns
		read(10,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,isOccult(i)
		read(20,*) lenj(i)
	enddo
	close(10)
	close(20)
	
	!Initialization of the total luminosities
	Lsy_tot_noabs = 0d0
    Lsy_tot_SSA = 0d0
    Lsy_tot_SSA_FF = 0d0
    Lsy_tot_SSA_FF_occ = 0d0

    open(40,file='SEDsyn_abs_cutsCJ.dat')
    open(41,file='SEDsyn_final_cutsCJ.dat')
    
    !Segment loop
	do j = 1,Ns

    	do i = 1,NSED_sy
			read(40,*) Eph(i),Lsy_noabs(i),Lsy_SSA(i),Lsy_SSA_FF(i)
		enddo
		Eph = Eph*eV			!Conversion to erg       	       
       	if (j < Ns) read(40,*)	!Blank space in the text file
       	
       	!Luminosity is null if occulted by the star
       	if (isOccult(j)) then
       		Lsy_SSA_FF_occ = 0d0
       	else
       	 	Lsy_SSA_FF_occ = Lsy_SSA_FF
       	endif
       	
       	!Total luminosities
       	do i = 1,NSED_sy	
       		
       		Lsy_tot_noabs(i) = Lsy_tot_noabs(i) + Lsy_noabs(i)					!No absorption
       		Lsy_tot_SSA(i) = Lsy_tot_SSA(i) + Lsy_SSA(i)						!SSA
       		Lsy_tot_SSA_FF(i) = Lsy_tot_SSA_FF(i) + Lsy_SSA_FF(i)				!SSA + FFA
       		Lsy_tot_SSA_FF_occ(i) = Lsy_tot_SSA_FF_occ(i) + Lsy_SSA_FF_occ(i)	!SSA + FFA + occultation
       		
    		!Jet contribution of the segments up to given lengths
    		if (lenj(j)-lenj(1) <= 0.1d0*a) Lsy1(i) = Lsy_tot_SSA_FF_occ(i)
    		if (lenj(j)-lenj(1) <= a) Lsy2(i) = Lsy_tot_SSA_FF_occ(i)
    		
			write(41,*) Eph(i)/eV,Lsy_noabs(i),Lsy_SSA(i),Lsy_SSA_FF(i),Lsy_SSA_FF_occ(i),j
       	  
       	enddo
       	write(41,*)	!Blank space      	
       	
	enddo
	
	close(40)
	close(41)
		
    LBol = 0d0
    LBolA = 0d0
    open(42,file='SEDsyn_final_CJ.dat')
    do i = 1,NSED_sy
    	write(42,*) Eph(i)/eV,Lsy1(i),Lsy2(i),Lsy_tot_noabs(i),Lsy_tot_SSA(i),Lsy_tot_SSA_FF(i),	&
    	&			Lsy_tot_SSA_FF_occ(i)
    	if (i /= NSED_sy) dE = Eph(i+1)-Eph(i)
    	LBol = LBol + Lsy_tot_noabs(i)/Eph(i)*dE        	  
    	LBolA = LBolA + Lsy_tot_SSA_FF_occ(i)/Eph(i)*dE	!Bolometric synchrotron luminosity
    enddo
    close(42)
    
    write(*,'(3X,A,ES11.4E2,1X,A)') 'C-Jet unabsorbed synchrotron luminosity (obs):',LBol,'erg/s'
    write(*,'(3X,A,ES11.4E2,1X,A)') 'C-Jet absorbed synchrotron luminosity (obs):',LBolA,'erg/s'
	
	deallocate (isOccult,lenj)
	return
	end subroutine Occult_SyCJ
	
!*************************************************************************************************************
!	Correction of the IC emission by occultation by the star (counter-jet)
!*************************************************************************************************************
	subroutine Occult_ICCJ (LIC_tot_GG_occ)
	real*8, intent(out) :: LIC_tot_GG_occ(NSED_IC)
	real*8 :: Eph(NSED_IC),LIC_noabs(NSED_IC),LIC_GG(NSED_IC),LIC_GG_occ(NSED_IC)
	real*8 :: LIC1(NSED_IC),LIC2(NSED_IC),LIC_tot_noabs(NSED_IC),LIC_tot_GG(NSED_IC)
	real*8 :: dummy,dE,LBol,LBolA
	allocate (isOccult(Ns),lenj(Ns))
	
	open(10,file='TrajectoryCJ.dat')
	open(20,file='Quant_FFCJ.dat')
	do i = 1,Ns
		read(10,*) dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,isOccult(i)
		read(20,*) lenj(i)
	enddo
	close(10)
	close(20)
	
	!Initialization of the total luminosities
    LIC_tot_noabs = 0d0
    LIC_tot_GG = 0d0
    LIC_tot_GG_occ = 0d0

    open(40,file='SEDIC_abs_cutsCJ.dat')
    open(41,file='SEDIC_final_cutsCJ.dat')
    
    !Segment loop
	do j = 1,Ns

        do i = 1,NSED_IC
        	read(40,*) Eph(i),LIC_noabs(i),LIC_GG(i)
       	enddo       	
       	Eph = Eph*eV			!Conversion to erg
       	if (j < Ns) read(40,*)	!Blank space in the text file
       	
       	!Luminosities are null if occulted by the star
       	if (isOccult(j)) then
       		LIC_GG_occ = 0d0
       	else
       		LIC_GG_occ = LIC_GG
       	endif
       	
       	!Total luminosities
       	do i = 1,NSED_IC      		
       		
       		LIC_tot_noabs(i) = LIC_tot_noabs(i) + LIC_noabs(i)			!No absorption
       		LIC_tot_GG(i) = LIC_tot_GG(i) + LIC_GG(i)					!GG
       		LIC_tot_GG_occ(i) = LIC_tot_GG_occ(i) + LIC_GG_occ(i)		!GG + occultation
       		
       	 	!Jet contribution of the segments up to given lengths
			if (lenj(j)-lenj(1) <= 0.1d0*a) LIC1(i) = LIC_tot_GG_occ(i)
			if (lenj(j)-lenj(1) <= a) LIC2(i) = LIC_tot_GG_occ(i)
       	  
			write(41,*) Eph(i)/eV,LIC_noabs(i),LIC_GG(i),LIC_GG_occ(i),j
       	  
       	enddo
       	write(41,*)	!Blank space       	
       	
	enddo
	
	close(40)
	close(41)
    
    LBol = 0d0
    LBolA = 0d0
    open(42,file='SEDIC_final_CJ.dat')
    do i = 1,NSED_IC
    	write(42,*) Eph(i)/eV,LIC1(i),LIC2(i),LIC_tot_noabs(i),LIC_tot_GG(i),LIC_tot_GG_occ(i)
		if (i /= NSED_IC) dE = Eph(i+1)-Eph(i)
		LBol = LBol + LIC_tot_noabs(i)/Eph(i)*dE   	  
		LBolA = LBolA + LIC_tot_GG_occ(i)/Eph(i)*dE		!Bolometric IC luminosity
    enddo
    close(42)

	write(*,'(3X,A,ES9.2E2,1X,A)') 'C-Jet unabsorbed IC luminosity (obs):',LBol,'erg/s'
    write(*,'(3X,A,ES9.2E2,1X,A)') 'C-Jet absorbed IC luminosity (obs):',LBolA,'erg/s'	
	
	deallocate (isOccult,lenj)
	return
	end subroutine Occult_ICCJ
	
!*************************************************************************************************************
!	Sum of the jet and counter-jet contributions
!*************************************************************************************************************
	subroutine TotalEmission ()
	real*8 :: Eph,Lj1,Lj2,Lj3,Lj4,Lj5,Lj6,Lcj1,Lcj2,Lcj3,Lcj4,Lcj5,Lcj6
	
	open(10,file='SEDsyn_final_jet.dat')
	open(11,file='SEDsyn_final_CJ.dat')
	open(12,file='SEDsyn_final.dat')
	
	open(20,file='SEDIC_final_jet.dat')
	open(21,file='SEDIC_final_CJ.dat')
	open(22,file='SEDIC_final.dat')
	
	do i = 1,NSED_sy
		read(10,*) Eph,Lj1,Lj2,Lj3,Lj4,Lj5,Lj6
		read(11,*) Eph,Lcj1,Lcj2,Lcj3,Lcj4,Lcj5,Lcj6
		write(12,*) Eph,Lj1+Lcj1,Lj2+Lcj2,Lj3+Lcj3,Lj4+Lcj4,Lj5+Lcj5,Lj6+Lcj6
	enddo
       	
	do i = 1,NSED_IC
		read(20,*) Eph,Lj1,Lj2,Lj3,Lj4,Lj5
		read(21,*) Eph,Lcj1,Lcj2,Lcj3,Lcj4,Lcj5
		write(22,*) Eph,Lj1+Lcj1,Lj2+Lcj2,Lj3+Lcj3,Lj4+Lcj4,Lj5+Lcj5
   	enddo       	
	
	close(10)
	close(11)
	close(12)
	close(20)
	close(21)
	close(22)
	
	return
	end subroutine TotalEmission
	
END MODULE occultation
