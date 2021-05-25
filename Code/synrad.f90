MODULE synrad
    use globals
    implicit none

contains

    !============================================================================
    ! SEDsyn (START)
    !============================================================================
    subroutine SEDsyn(jet)
    !============================================================================
    ! This subroutine calculates the synchrotron SED 
    ! The input is the '' (jet) or 'CJ' (counter-jet). 
    !============================================================================
    character(LEN=*), intent(in) :: jet ! 'CJ' or 'jet'
    real(kind(1.d0)) :: E(nIn), NE(nIn,Ns), NEcut(nIn), Ein, Eph, dEph, EphFF(Ns), dEphFF(Ns)
    real(kind(1.d0)) :: Lumcut(Ns), LumcutFF(Ns), L1, L2, L3, L4, L5, L6, Lum, LumFF, LBolFF
    real(kind(1.d0)) :: LumcutA(Ns), LumcutAFF(Ns), LumA, LumAFF, LBolAFF
    real(kind(1.d0)) :: boost(Ns)
    allocate ( delta_obs(Ns) )

    !Quantities at each segment
    open(20,file='deltas'//jet//'.dat')
    do i = 1,Ns
        read(20,*) delta_obs(i), boost(i)
    enddo
    close(20)
    
    !N(E) of each cut
    open(30,file='NE_cuts.dat')
    do i = 1, nIn
        do j = 1, Ns
            read(30,*) E(i), NE(i,j), dummy
        enddo
    enddo
    close(30)

    open(40,file='SEDsyn_cuts'//jet//'.dat')
    open(50,file='SEDsyn_jet'//jet//'.dat')
    
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
        call SSA(BFF,rj,sjoFF,VFF,EphFF,E,NEcut,LumcutFF,LumcutAFF)
        
        do j = 1, Ns
        
            LumcutA(j) = boost(j)*LumcutAFF(j)    !(de)Boosted
            write(40,*) Eph/eV, Eph*Lumcut(j), Eph*LumcutA(j), j
            
            !Absorbed (de)boosted luminosity for Eph
            LumA = LumA + LumcutA(j)
            
            LumAFF = LumAFF + LumcutAFF(j)*dEphFF(j)
            
            !Jet and c-jet contribution of the segments up to given lengths
            if (z(j) <= 7.3*R_g) L1 = LumA
            if (z(j) <= 2.6d1*R_g) L2 = LumA
            if (z(j) <= 6.5d1*R_g) L3 = LumA
            if (z(j) <= 1.95d2*R_g) L4 = LumA
            if (z(j) <= 6.5d2*R_g) L5 = LumA
            if (z(j) <= 3d5*R_g) L6 = LumA
            
        enddo

        write(50,*) Eph/eV, Eph*L1, Eph*L2, Eph*L3, Eph*L4, Eph*L5, Eph*L6, Eph*Lum, Eph*LumA
        
        !Absorbed bolometric luminosity
        LBolAFF = LBolAFF + LumAFF
        
        Eph = Eph*Ein       !Increase in the photon energy

    enddo
    
    !write(*,'(3X,A,ES10.3E2,1X,A)') 'Jet unabsorbed synchrotron luminosity (FF):',LBolFF,'erg/s'
    !write(*,'(3X,A,ES10.3E2,1X,A)') 'Jet absorbed synchrotron luminosity (FF):',LBolAFF,'erg/s'

    close(40)
    close(50)
    
    deallocate( delta_obs )
    
    end subroutine SEDsyn
    !============================================================================
    ! SEDsyn (END)
    !============================================================================
    

    !************************************************************************************
    !   Synchrotron luminosity per photon energy unit
    !************************************************************************************
    real(kind(1.d0)) function SyncLum(B,Eph,E,NE)
    real(kind(1.d0)), intent(in) :: B, Eph, E(nIn), NE(nIn)
    real(kind(1.d0)) :: dE, Ec, x, Psync, E_step
    real(kind(1.d0)), parameter :: cte_Ec = sqrt(3d0/2d0)*(qe*h)/(2d0*pi*me**3*c**5)
    real(kind(1.d0)), parameter :: cte_sy = 1.85d0*sqrt(2d0)*qe**3/(h*me*c**2) 
    !Initialization
    SyncLum = 0d0
    !Integral by rectangles
    E_step = E(2)/E(1)
    do i = 1, nIn
        dE = E(i)*(E_step - 1.d0)     !Energy interval
        !Critical photon energy and power emitted for each particle energy
        Ec = cte_Ec*B*E(i)**2
        x = Eph/Ec
        Psync = cte_sy*B * x**(1d0/3d0)*exp(-x)
        !The approximation is only valid for 0.1 < x < 10
        if ((x < 1d-1) .or. (x > 1d1)) Psync = 0d0
        SyncLum = SyncLum + Psync*NE(i)*dE
    enddo
    end function SyncLum
    

    !============================================================================
    ! SSA (start)
    !============================================================================
    subroutine SSA (B,r_j,sjo,V,Eph,E,NE,Lum,LumA)
    !============================================================================
    ! Luminosity corrected for synchrotron self-absorption
    !============================================================================
    real(kind(1.d0)), intent(in) :: B(Ns),r_j(Ns),sjo(Ns),V(Ns),Eph(Ns),E(nIn),NE(nIn),Lum(Ns)
    real(kind(1.d0)), intent(out) :: LumA(Ns)
    real(kind(1.d0)) :: K,Ecf,Sf,depth,tau,Ap

    !Loop for the segments
    do j = 1, Ns

        !Absorption and emission coefficients, and source function, for each segment (Pacholczyck, 1970)
        K = 2.62d0*qe**3*h**2*B(j)/(4d0*pi*Eph(j)**2*me*V(j)) * EintSSA(B(j),Eph(j),E,NE)
        Ecf = Lum(j)/(4d0*pi*V(j))      !Per stereoradian
        Sf = Ecf/max(K,1.d-30) !!!! CHECK !!!!

        !Length that photons in the observer's direction must travel inside the segment (in the FF)
        depth = r_j(j)/sjo(j)

        !Optical depth for SSA
        tau = K*depth / delta_obs(j)
        !tau = K*min(depth,z(j)/gj)!/delta_obs(j)
        !Estimation of the projected area, perpendicular to the observer's direction, obtained by equating
        !the unabsorbed Lum to its value for optically thin emission: 4*pi*Ap*Sf*tau
        !Ap = V(j)/min(depth,z(j)/gj)
        Ap = V(j)*delta_obs(j)/depth

        !Luminosity corrected for synchrotron self-absorption
        if (tau < 0.01d0) then  !!!! CHECK !!!!
            LumA(j) = Lum(j)
        else
            LumA(j) = 4d0*pi*Ap*Sf*(1d0-exp(-tau))
        endif
        !print*, tau, Lum(j), 4d0*pi*Ap*Sf*(1d0-exp(-tau))
    enddo

    return
    end subroutine SSA
    !============================================================================
    ! SSA (END)
    !============================================================================


    !************************************************************************************
    !   Integral for the particle energy in SSA
    !************************************************************************************
    real(kind(1.d0)) function EintSSA(B,Eph,E,NE)
    real(kind(1.d0)), intent(in) :: B, Eph, E(nIn), NE(nIn)
    real(kind(1.d0)) :: dE, Ec, x, F, E_step
    real(kind(1.d0)), parameter :: cte_Ec = sqrt(3d0/2d0)*(qe*h)/(2d0*pi*me**3*c**5)
    EIntSSA = 0d0   !Initialization
    E_step = E(2)/E(1)
    do i = 1, nIn
        dE = E(i)*(E_step - 1.d0)     !Energy interval
        !Critical photon energy and function to integrate for each particle energy
        Ec = cte_Ec*B*E(i)**2
        x = Eph/Ec
        F = NE(i)/E(i)*exp(-x)*x**(1d0/3d0)*(x + 2d0/3d0)
        !The approximation is only valid for 0.1 < x < 10
        if (x < 1d-3 .or. x > 1d1) F = 0d0 !!!!CHECK THIS!!!!, troublesome in first cell?
        EintSSA = EintSSA + F*dE
    enddo
    end function EintSSA


    !============================================================================
    ! SEDsyn (END)
    !============================================================================
    subroutine Calc_Usyn(jet)
    !============================================================================
    ! Computation of the synchrotron energy density in each segment
    ! The input is the '' (jet) or 'CJ' (counter-jet). 
    !============================================================================
        character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet 
        real(kind(1.d0)) :: L_tot(Ns), L_abs(Ns), Lun(NSED_sy, Ns), Lsy(NSED_sy, Ns)
        real(kind(1.d0)) :: Esy(NSED_sy), Esy_step, dEsy
        real(kind(1.d0)) :: boost(Ns), U_sy_FF

        ! Read relevant quantities at each segment
        open(20, file='deltas'//jet//'.dat')
        do l = 1, Ns
            read(20,*) dummy, boost(l)
        enddo
        close(20)

        ! Read from a file: Esy(i,l), E*Lsy_unabs(i,l), E*Lsy_abs(i,l), l
        open(27, file='SEDsyn_cuts_ord'//jet//'.dat')
        do l = 1, Ns 
            do i = 1, NSED_sy
                read(27,*) Esy(i), Lun(i,l), Lsy(i,l), dummy
                ! We need L_sy_FF, so we unboost
                Lsy(i,l) = Lsy(i,l)/(Esy(i)*boost(l))
                Lun(i,l) = Lun(i,l)/(Esy(i)*boost(l))
            enddo
        enddo
        close(27)

        ! Integrate the luminosity of each segment
        ! Also calculate the total absorbed luminosity of each segment
        L_tot = 0.d0
        L_abs = 0.d0
        Esy_step = Esy(2)/Esy(1)
        open(28, file='U_FF'//jet//'.dat')
        open(29, file='L_inj'//jet//'.dat')
        do l = 1, Ns
            do i = 1, NSED_sy
                dEsy = Esy(i)*(Esy_step-1.d0)
                L_tot(l) = L_tot(l) + Lsy(i,l) * dEsy
                L_abs(l) = L_abs(l) + (Lun(i,l) - Lsy(i,l)) * dEsy
            enddo
            U_sy_FF = ( L_tot(l)/VFF(l) ) * ( rj(l)/c ) !!! CHECK IF <t> = r/c, deltaz/c, or what
            write(28,*) l, z(l), U_sy_FF, BFF(l)**2/(8.d0*pi), U_sy_FF/(BFF(l)**2/(8.d0*pi))
            write(29,*) l, z(l)/R_g, L_tot(l), L_tot(l)/LnTFF
        enddo
        close(28)
        close(29)

    end subroutine Calc_Usyn
    !============================================================================
    ! Calc_Usyn (END)
    !============================================================================


    !============================================================================
    ! SEDsyn_data (START)
    !============================================================================
    subroutine SEDsyn_data(E_data, S_data, jet)
    !============================================================================
    ! This routine calculates the intrinsic synchrotron emission in the fluid frame
    ! at the energies E_data [cgs] (vector), and returns the specific flux at each 
    ! of those energies [cgs] (also a vector). The number of points in E_data has
    ! to be smaller than m_data_max defined in globals.f90.
    ! The input is the '' (jet) or 'CJ' (counter-jet). 
    !============================================================================
    real(kind(1.d0)), intent(in) :: E_data(:)
    real(kind(1.d0)), intent(out) :: S_data(m_data_max)
    character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet
    real(kind(1.d0)), allocatable :: Eph(:), LumA(:)
    real(kind(1.d0)) :: E(nIn), NE(nIn,Ns), NEcut(nIn), EphFF(Ns)
    real(kind(1.d0)) :: LumcutFF(Ns), LumcutAFF(Ns), LumcutA(Ns)
    real(kind(1.d0)) :: boost(Ns), dil
    integer :: m_data
    allocate( delta_obs(Ns) )

    !Quantities at each segment
    open(20,file='deltas'//jet//'.dat')
    do i = 1,Ns
        read(20,*) delta_obs(i), boost(i)
    enddo
    close(20)
    
    !N(E) of each cut
    open(30,file='NE_cuts.dat')
    do i = 1, nIn
        do j = 1, Ns
            read(30,*) E(i), NE(i,j), dummy
        enddo
    enddo
    close(30)

    m_data = SIZE(E_data)
    if (m_data > m_data_max) then
      print*, 'Too many points to fit, make m_data_max larger in global.f90!'
      stop
    endif
    allocate( Eph(m_data), LumA(m_data) )
    Eph = E_data ! Here we can change units if needed          

    do l = 1, m_data
        do j = 1, Ns                        !Loop for the jet segments
            EphFF(j) = Eph(l)/delta_obs(j)  !Photon energy in the FF
            do i = 1, nIn
                NEcut(i) = NE(i,j)          !N(E) for the segment j
            enddo
            !Unabsorbed synchrotron luminosity for each segment
            LumcutFF(j) = SyncLum(BFF(j),EphFF(j),E,NEcut)
        enddo
        
        !Absorbed radiation for each cut
        call SSA(BFF,rj,sjoFF,VFF,EphFF,E,NEcut,LumcutFF,LumcutAFF)
        
        LumA(l) = 0d0 ! Total luminosity at each photon energy
        do j = 1, Ns
            LumcutA(j) = boost(j)*LumcutAFF(j)  !(de)Boosted
            LumA(l) = LumA(l) + LumcutA(j)            !Absorbed (de)boosted luminosity for Eph
        enddo
    enddo

    ! Define the output vector for the first m_data positions
    S_data = 0.d0
    dil = 4.d0*pi*dist**2
    do l = 1, m_data
      S_data(l) = LumA(l)/dil
    enddo
    
    deallocate( Eph, LumA )
    deallocate( delta_obs )
    
    end subroutine SEDsyn_data
    !============================================================================
    ! SEDsyn (END)
    !============================================================================
    

!================================================================================


END MODULE synrad
