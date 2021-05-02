MODULE edistribution
    use globals
    use losses
    implicit none

contains

!*************************************************************************************************************
!   Creation of a text file E-E²*N(E) for plotting
!*************************************************************************************************************
    subroutine Edist()
    real(dp) :: Kinj,Ecut,Ein,EFF,NE(Ns),NEtot,NE1,NE2
    real(dp) :: dummy
    allocate ( z(nS), BFF(Ns), rj(Ns), vnorm(Ns), EphMFF(Ns), Kappa(Ns))
    
    write(*,*) '\n Calculating the electron distribution along the jet'
        
    open(10,file='quant_FF.dat')
    do i = 1, Ns
        !Quantities at each segment
        read(10,*) z(i), BFF(i), rj(i), vnorm(i), Kappa(i), dummy, dummy, EphMFF(i)
    enddo
    close(10)
    
    Ecut = CutE(BFF(1),rj(1),vnorm(1),EphMFF(1),Kappa(1))   !Cutoff energy for Q(E) in the FF
    Kinj = NormK(Ecut)                                      !Normalization constant for Q(E)
    write(*,'(X,A,X,F4.2,X,A)') 'Ecutoff =',Ecut/eV/1d12,'TeV'
 
    open(20,file='NE_cuts.dat')
    open(30,file='NE_jet.dat')
    
    Ein = (E2/E1)**(1d0/(nNE-1))    !Energy increase factor
    EFF = E1                        !Initialization

    do i2 = 1,nNE
        
        !Computation of N(E) for each segment   
        call CalcNE(Kinj,EFF,Ecut,z,BFF,rj,vnorm,EphMFF,Kappa,NE)
        
        NEtot = 0d0     !Initialization

        do j = 1,Ns

            write(20,*) EFF/eV,EFF**2*NE(j),j   !Write E-E²*N(E) for each jet cut

            NEtot = NEtot + NE(j)   !Sum for all the jet segments

            !Contribution of the segments up to given lengths
            if (z(j) <= 0.01d0*zmax) NE1 = NEtot
            if (z(j) <= 0.1d0*zmax) NE2 = NEtot

        enddo

        !Write E-E²*N(E) for the whole jet
        write(30,*) EFF/eV,EFF**2*NE1,EFF**2*NE2,EFF**2*NEtot
        !write(*,'(4ES13.3E2)') EFF/eV,EFF**2*NE1,EFF**2*NE2,EFF**2*NEtot
        
        EFF = EFF*Ein   !Exponential energy increase

    enddo
    
    close(20)
    close(30)
    
88  format (4X,ES9.3E2,5X,ES9.3E2,5X,ES9.3E2,5X,ES9.3E2)

    deallocate (z,BFF,rj,vnorm,EphMFF,Kappa)

    end subroutine Edist
    
!*************************************************************************************************************
!   Creation of a text file E-N(E) for the integrals
!*************************************************************************************************************
    subroutine Edist_int ()
    real(dp) :: Kinj,Ecut,Ein,EFF,NE(Ns),NEtot
    real(dp) :: dummy
    allocate (z(Ns),BFF(Ns),rj(Ns),vnorm(Ns),EphMFF(Ns),Kappa(Ns))
        
    open(10,file='quant_FF.dat')
    do i = 1,Ns
        !Quantities at each segment
        read(10,*) z(i), BFF(i), rj(i), vnorm(i), Kappa(i), dummy, dummy, EphMFF(i)
    enddo
    close(10)
    
    Ecut = CutE(BFF(1),rj(1),vnorm(1),EphMFF(1),Kappa(1))   !Cutoff energy for Q(E) in the FF
    Kinj = NormK(Ecut)                                      !Normalization constant for Q(E)
 
    open(20,file='NE_intcuts.dat')
    
    Ein = (E2/E1)**(1d0/(nIn-1))    !Energy increase factor
    EFF = E1                        !Initialization
    
    do i2 = 1,nIn
    
        !Computation of N(E) for each segment
        call CalcNE(Kinj,EFF,Ecut,z,BFF,rj,vnorm,EphMFF,Kappa,NE)
        
        do j = 1,Ns
            write(20,*) EFF,NE(j)   !E and N(E) for each jet cut
        enddo
        
        EFF = EFF*Ein   !Exponential energy increase
        
    enddo
    
    close(20)
    
    write(*,'(X,A)') 'Particle distribution computed'
    
    deallocate (z,BFF,rj,vnorm,EphMFF,Kappa)
    
    end subroutine Edist_int
    
!*************************************************************************************************************
!   Numerical particle distribution by evolving each cut as a 1-zone
!*************************************************************************************************************
    subroutine CalcNE (Kinj,EFF,Ecut,z,BFF,rj,vnorm,EphMFF,Kappa,NE)
    real(dp), intent(in) :: Kinj,EFF,Ecut,z(Ns),BFF(Ns),rj(Ns),vnorm(Ns),EphMFF(Ns),Kappa(Ns)
    real(dp), intent(out) :: NE(Ns)
    real(dp) :: tlossFF,E0i,E0f,E0(Ns),backD,dtpast,prod

    !Initialization of the current energy distribution
    NE = 0d0  !Applies to the whole vector

    !Loop for the injection
    do i = 1,Ns

        !Initial value of the energy for each injection
        E0(i) = EFF

        !Initial value of the lower energy limit for the cut-by-cut one-zone computation
        E0i = EFF
        
        !Initialization
        backD = 0d0 !Backwards distance
        prod = 1d0
        j = i       !Reference for the product
        n = i       !Reference segment for B,nf,rj...

        !Loop for the distance that each population has made
        do while (backD <= z(i)-z(1))
        
            !Null contribution to N(E) if E0 > Emax (not physically possible)
            if (E0(i) > EmaxFF) exit

            !Change the reference segment when backD is long enough
            if (backD > (z(i)-z(n))) n = n - 1
            
            !Selection of the adequate time steps for each position (dtpast already in the FF)
            call Setdt(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0(i),dtpast)
            
            !Increase in the energy
            E0(i) = E0(i) - TLosses(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0(i))*dtpast
            
            !Condition to find the point where the segment i was at the current position of segment j
            if (backD > z(i)-z(j)) then

                !Upper energy limit at the segment j position
                E0f = E0(i)

                !Product of the quotient of the losses
                prod = prod * TLosses(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0f)/        &
                &             TLosses(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0i)

                !Lower energy limit at the segment j-1 position
                E0i = E0(i)

                !Decrease in the segment count
                j = j - 1

            endif

            !Increase in the backwards distance (LAB frame)
            backD = backD + vj*dtpast*gj

        enddo
        
        !Cooling time for energy E0 and segment i
        tlossFF = E0(i)/abs(TLosses(BFF(i),rj(i),vnorm(i),EphMFF(i),Kappa(i),E0(i)))
        
        !Energy distribution for each segment, corrected when tloss' < dt' for a given energy (to avoid 
        !overestimation at high energies due to the 1st segment being too long)
        NE(i) = Q_E(Kinj,Ecut,E0(i))*dtFF*prod * min(1d0,tlossFF/dtFF)
        
        !Correction in N(E) for the first cut when E0 is close to Emax
        if (i == 1 .and. EmaxFF < 2d0*E0(i)) NE(i) = NE(i) * (EmaxFF-E0(i)) * tlossFF/E0(i)     
        
        !if (tlossFF < dtFF .and. E0(i) <= EmaxFF) print *, i,EFF/eV,E0(i)/eV,tlossFF
        
    enddo

    return
    end subroutine CalcNE
    
!*************************************************************************************************************
!   Setting of an appropriate time step for the N(E) computation
!*************************************************************************************************************
    subroutine Setdt (B,rj,vnorm,Eph,Kappa,E,dtpast)
    real(dp), intent(in) :: B,rj,vnorm,Eph,Kappa,E
    real(dp), intent(out) :: dtpast
    real(dp) :: tloss

    !Characteristic time for losses
    tloss = dabs(E/TLosses(B,rj,vnorm,Eph,Kappa,E))

    !Set dt according to the characterstic losses and injection times
    dtpast = 0.3d0 * min(tloss,dtFF)

    return
    end subroutine Setdt
    
!*************************************************************************************************************
!   Determination of the total energy in relativistic particles at each position
!*************************************************************************************************************
    subroutine Int_EN_j ()
    real(dp) :: E(nIn),NE(nIn,Ns),dE,E_NT,sumE
    real(dp) :: x_tr,y_tr,z_tr
   
    open(70,file='NE_intcuts.dat')
    open(80,file='Trajectory.dat')
    open(90,file='ENT_trajectory.dat')
    
    do i = 1,nIn
        do j = 1,Ns
            read(70,*) E(i),NE(i,j)
        enddo
    enddo
    
    close(70)
   
    sumE = 0d0
    do j = 1,Ns
    
        read(80,*) x_tr,y_tr,z_tr
        
        E_NT = 0d0      !Total energy in NT particles at each segment
        
        do i = 1,nIn
            if (i /= nIn) dE = E(i+1)-E(i)  
            E_NT = E_NT + E(i)*NE(i,j)*dE
        enddo
                
        write(90,*) j,x_tr,y_tr,z_tr,E_NT
        
        sumE = sumE + E_NT      !Total jet energy in NT particles
        
    enddo
   
    close(80)
    close(90)
   
    end subroutine Int_EN_j
    
!*************************************************************************************************************
!   Determination of the cutoff energy of the particle injection
!*************************************************************************************************************
    real(dp) function CutE (BFF1,rj1,vnorm1,Eph1,Kappa1)
    real(dp), intent(in) :: BFF1,rj1,vnorm1,Eph1,Kappa1
    real(dp) :: Emax1,Emax2,Eacc,lE1,lE2,lEmid,error,check

    !Injected power (erg/s)
    Eacc = facc*qe*BFF1*c

    !Initialization
    lE1 = log10(1d9*eV)  ! log10(Emin [CGS]) ~ 1 GeV
    lE2 = log10(1d14*eV) ! log10(Emax [CGS]) ~ 100 TeV
    error = 10d0**(-4) ! error < 10^(-4)
        
    !Maximum particle energy
    do while (abs(lE2-lE1) > error)
        lEmid = (lE1+lE2)/2d0
        check =  Eacc - abs(TLosses(BFF1,rj1,vnorm1,Eph1,Kappa1,10d0**lEmid))
        if (check == 0d0) then
            lE2 = lEmid
            exit
        else if (check < 0d0) then
            lE2 = lEmid
        else if (check > 0d0) then
            lE1 = lEmid
        endif
    enddo
    Emax1 = 10d0**(lE2)

    !Maximum energy that particles can get before diffusing away of the injection point
    Emax2 = rj1*sqrt(3d0*qe*BFF1*Eacc/(2d0*c))

    !Ecut obtained by the minimum of the values above
    CutE = min(Emax1,Emax2)

    end function CutE
    
!*************************************************************************************************************
!   Injection function
!*************************************************************************************************************
    real(dp) function Q_E (K,Ecut,E)
    real(dp), intent(in) :: K,Ecut,E
    
    !Q(E) = 0 if outside the injection energy limits
    if (E <= EmaxFF .and. E >= EminFF) then
        Q_E = K*E**(-p) * dexp(-E/Ecut)
    else
        Q_E = 0d0
    endif

    return
    end function Q_E
  
!*************************************************************************************************************
!   Injection normalization
!*************************************************************************************************************  
    real(dp) function NormK (Ecut)
    real(dp), intent(in) :: Ecut
    real(dp) :: E,Ein,dE,QInt
  
    N = 1000                            !Number of subintervals for the integral
    Ein = (EmaxFF/EminFF)**(1d0/(N-1))  !Increase factor of the energy
  
    !Initialization
    E = EminFF
    QInt = 0d0

    ! USE TRAPEZOID RULE INSTEAD!!!!!
    !Integration by the Simpson's rule with variable step size
    do i = 1,N  
        dE = E*(Ein-1d0)                    !Energy interval
        QInt = QInt + E*Q_E(1d0,Ecut,E)*dE  !Integral of E*Q(E)
        E = E*Ein                           !Increase in the energy        
    enddo
  
    NormK = LnTFF/QInt    !Normalization constant
  
    end function NormK

END MODULE edistribution
