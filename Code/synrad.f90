MODULE synrad
    use globals
    implicit none

contains

    !***************************************************************************************************
    !   Computation of the synchrotron SED (jet)
    !***************************************************************************************************
    subroutine SEDsyn()
    real(dp) :: E(nIn),NE(nIn,Ns),NEcut(nIn),Ein,Eph,dEph,EphFF(Ns),dEphFF(Ns),Lumcut(Ns),LumcutFF(Ns), &
    &         L1,L2,Lum,LumFF,LBolFF,LumcutA(Ns),LumcutAFF(Ns),LumA,LumAFF,LBolAFF
    real(dp) :: dummy, boost(Ns)
    allocate ( z(nS), BFF(Ns), rj(Ns), sjoFF(Ns), VFF(Ns), delta_obs(Ns) )
    
    write(*,*) '\n Calculating the synchrotron SED along the jet'
        
    !Quantities at each segment
    open(10,file='quant_FF.dat')
    open(20,file='deltas.dat')
    do i = 1,Ns
        read(10,*) z(i), BFF(i), rj(i), dummy, dummy, sjoFF(i), VFF(i), dummy
        read(20,*) delta_obs(i), boost(i)
    enddo
    close(10)
    close(20)
    
    !N(E) of each cut
    open(30,file='NE_intcuts.dat')
    do i = 1,nIn
        do j = 1, Ns
            read(30,*) E(i),NE(i,j)
        enddo
    enddo
    close(30)

    open(40,file='SEDsyn_cuts.dat')
    open(50,file='SEDsyn_jet.dat')

    write(*,'(1X,A,I3,A)') 'Synchrotron SED (',NSED_sy,')'
    
    Ein = (Emax_sy/Emin_sy)**(1d0/(NSED_sy-1))  !Energy increase factor
    Eph = Emin_sy                               !Initialization

    !Initialization of the bolometric luminosities
    LBolFF = 0d0
    LBolAFF = 0d0

    do l = 1, NSED_sy

        !Initialization of the total luminosity for each photon energy
        Lum = 0d0
        LumFF = 0d0
        LumA = 0d0
        LumAFF = 0d0
        
        dEph = Eph*(Ein-1d0)    !Energy interval

        !Loop for the segments
        do j = 1, Ns
        
            EphFF(j) = Eph/delta_obs(j)     !Photon energy in the FF
            dEphFF(j) = dEph/delta_obs(j)   !Energy interval in the FF
        
            do i = 1, nIn
                !N(E) for the segment j
                NEcut(i) = NE(i,j)
            enddo

            !Unabsorbed synchrotron luminosity for each segment
            LumcutFF(j) = SyncLum(BFF(j),EphFF(j),E,NEcut)
            Lumcut(j) = boost(j)*LumcutFF(j)      !(de)Boosted

            !Unabsorbed (de)boosted luminosity for Eph
            Lum = Lum + Lumcut(j)
            
            LumFF = LumFF + LumcutFF(j)*dEphFF(j)

        enddo
        
        !Unabsorbed bolometric luminosity
        LBolFF = LBolFF + LumFF
                
        !Absorbed radiation for each cut
        call SSA(z,BFF,rj,sjoFF,VFF,EphFF,E,NEcut,LumcutFF,LumcutAFF)
        
        do j = 1, Ns
        
            LumcutA(j) = boost(j)*LumcutAFF(j)    !(de)Boosted
            write(40,*) Eph/eV, Eph*Lumcut(j), Eph*LumcutA(j), j
            
            !Absorbed (de)boosted luminosity for Eph
            LumA = LumA + LumcutA(j)
            
            LumAFF = LumAFF + LumcutAFF(j)*dEphFF(j)
            
            !Jet and c-jet contribution of the segments up to given lengths
            if (z(j) <= 0.01d0*zmax) L1 = LumA
            if (z(j) <= 0.1d0*zmax) L2 = LumA
            
        enddo

        write(50,*) Eph/eV, Eph*L1, Eph*L2, Eph*Lum, Eph*LumA
        
        !Absorbed bolometric luminosity
        LBolAFF = LBolAFF + LumAFF
        
        Eph = Eph*Ein       !Increase in the photon energy

    enddo
    
    !write(*,'(3X,A,ES13.6E2,1X,A)') 'Jet unabsorbed synchrotron luminosity (FF):',LBolFF,'erg/s'
    !write(*,'(3X,A,ES13.6E2,1X,A)') 'Jet absorbed synchrotron luminosity (FF):',LBolAFF,'erg/s'

    close(40)
    close(50)
    
    deallocate( z, BFF, rj, sjoFF, VFF, delta_obs )
    
    end subroutine SEDsyn
    

    !*************************************************************************************************************
    !   Synchrotron luminosity per photon energy unit
    !*************************************************************************************************************
    real(dp) function SyncLum (B,Eph,E,NE)
    real(dp), intent(in) :: B,Eph,E(nIn),NE(nIn)
    real(dp) :: dE(nIn),Ec,x,Psync
    real(dp), parameter :: cte_Ec = sqrt(3d0/2d0)*(qe*h)/(2d0*pi*me**3*c**5)
    real(dp), parameter :: cte_sy = 1.85d0*sqrt(2d0)*qe**3/(h*me*c**2) 

    !Initialization
    SyncLum = 0d0

    !Integral by rectangles
    do i = 1, nIn

        if (i /= nIn) dE(i) = E(i+1) - E(i)     !Energy interval

        !Critical photon energy and power emitted for each particle energy
        Ec = cte_Ec*B*E(i)**2
        x = Eph/Ec
        Psync = cte_sy*B * x**(1d0/3d0)*exp(-x)

        !The approximation is only valid for 0.1 < x < 10
        if ((x < 1d-1) .or. (x > 1d1)) Psync = 0d0

        SyncLum = SyncLum + Psync*NE(i)*dE(i)

    enddo

    return
    end function SyncLum
    
    !*************************************************************************************************************
    !   Luminosity corrected for synchrotron self-absorption
    !*************************************************************************************************************
    subroutine SSA (z,B,rj,sjo,V,Eph,E,NE,Lum,LumA)
    real(dp), intent(in) :: z(Ns),B(Ns),rj(Ns),sjo(Ns),V(Ns),Eph(Ns),E(nIn),NE(nIn),Lum(Ns)
    real(dp), intent(out) :: LumA(Ns)
    real(dp) :: K,Ecf,Sf,depth,tau,Ap

    !Loop for the segments
    do j = 1, Ns

        !Absorption and emission coefficients, and source function, for each segment (Pacholczyck, 1970)
        K = 2.62d0*qe**3*h**2*B(j)/(4d0*pi*Eph(j)**2*me*V(j)) * EIntSSA(B(j),Eph(j),E,NE)
        Ecf = Lum(j)/(4d0*pi*V(j))      !Per stereoradian
        Sf = Ecf/K

        !Length that photons in the observer's direction must travel inside the segment (in the FF)
        depth = rj(j)/sjo(j)

        !Optical depth for SSA
        tau = K*min(depth,z(j)/gj)

        !Estimation of the projected area, perpendicular to the observer's direction, obtained by equating
        !the unabsorbed Lum to its value for optically thin emission: 4*pi*Ap*Sf*tau
        Ap = V(j)/min(depth,z(j)/gj)

        !Luminosity corrected for synchrotron self-absorption
        if (tau < 0.1d0) then
            LumA(j) = Lum(j)
        else
            LumA(j) = 4d0*pi*Ap*Sf*(1d0-exp(-tau))
        endif

    enddo

    return
    end subroutine SSA

!*************************************************************************************************************
!   Integral for the particle energy in SSA
!*************************************************************************************************************
    real(dp) function EintSSA(B,Eph,E,NE)
    real(dp), intent(in) :: B, Eph, E(nIn), NE(nIn)
    real(dp) :: dE(nIn), Ec, x, F
    real(dp), parameter :: cte_Ec = sqrt(3d0/2d0)*(qe*h)/(2d0*pi*me**3*c**5)
    
    EIntSSA = 0d0   !Initialization
    
    do i = 1, nIn

        if (i /= nIn) dE(i) = E(i+1) - E(i)     !Energy interval

        !Critical photon energy and function to integrate for each particle energy
        Ec = cte_Ec*B*E(i)**2
        x = Eph/Ec
        F = NE(i)/E(i)*exp(-x)*x**(1d0/3d0)*(x + 2d0/3d0)

        !The approximation is only valid for 0.1 < x < 10
        if (x < 1d-1 .or. x > 1d1) F = 0d0

        EintSSA = EintSSA + F*dE(i)

    enddo

    return
    end function EintSSA



    !***************************************************************************************************
    !  Computation of the synchrotron energy density in each segment
    !***************************************************************************************************
    subroutine Calc_Usyn()
        real(dp) :: L_tot(Ns), Lsy(nIn, Ns), Esy(nIn), Esy_step, dEsyn
        real(dp) :: z(Ns), Bff(Ns), r(Ns), VFF(Ns), boost(Ns), dummy, U_sy_FF

        ! Read relevant quantities at each segment
        open(10, file='quant_FF.dat')
        open(20, file='deltas.dat')
        do l = 1, Ns
            read(10,*) z(l), BFF(l), r(l), dummy, dummy, dummy, VFF(l), dummy
            read(20,*) dummy, boost(l)
        enddo
        close(10)
        close(20)

        ! Read from a file: Esy(i,l), Lsy_unabs(i,l), Lsy_abs(i,l), l
        open(27, file='SEDsyn_cuts_ord.dat')
        do l = 1, Ns 
            do i = 1, nIn
                read(27,*) Esy(i), dummy, Lsy(i,l), dummy
            enddo
        enddo
        close(27)

        ! Integrate the luminosity of each segment
        L_tot = 0.d0
        Esy_step = Esy(2)/Esy(1)
        open(28, file='U_FF.dat')
        do l = 1, Ns
            do i = 1, nIn
                dEsyn = Esy(i)*(Esy_step-1.d0)
                L_tot(l) = L_tot(l) + Lsy(i,l) * dEsyn
            enddo
            L_tot(l) = L_tot(l)/boost(l) ! We need L_sy_FF, so we unboost
            U_sy_FF = ( L_tot(l)/VFF(l) ) * ( r(l)/c )
            write(28,*) l, z(l), U_sy_FF, BFF(l)**2/(8.d0*pi), U_sy_FF/(BFF(l)**2/(8.d0*pi))
        enddo
        close(28)

    end subroutine Calc_Usyn


!==============================================================


END MODULE synrad
