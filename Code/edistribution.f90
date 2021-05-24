MODULE edistribution
    use globals
    use losses
    implicit none

contains
    
   !============================================================================
   ! Edist (START)
   !============================================================================
    subroutine Edist()
   !============================================================================
   ! This subroutine calculates the electron energy distribution at each segment
   !============================================================================
    real(dp) :: Kinj, Ecut, Ein, EFF, NE(Ns), NEtot, NE1, NE2
    real(dp) :: NEmatrix(nIn, Ns), NEmat(Ns, nIn), Evec(nIn)
    real(dp) :: NEE(Ns)
    integer :: izz
        
    open(11,file='Emax_FF.dat')
    do i = 1,Ns
        !Quantities at each segment
        Ecut = CutE(BFF(i),rj(i),vnorm(i),EphMFF(i),Kappa(i))   !Cutoff energy for Q(E) in the FF
        if (Ecut > EmaxFF/3.d0) then
            write(*,*) 'EmaxFF in globals.f90 is too low!'
            stop
        else
            write(11,*) i, z(i), Ecut/eV
        endif
    enddo
    close(11)
    
    open(20,file='NE_cuts.dat')
    open(30,file='NE_jet.dat')
    
    ! Calculate logarithmically-spaced energy vector
    call vector_log( EminFF, EmaxFF, Ein, Evec ) 

    do izz = 1,Ns
        NEE(izz) = 0d0
    enddo

    !do i2 = 1, nIn
    !    call CalcNE_analytical(EFF, z, BFF, rj, vnorm, EphMFF, Kappa, NEmatrix, i2)
    !    EFF = EFF * Ein
    !enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ESTA ENERGIA DE CUTOFF DEBERIA CALCULARSE PARA CADA CELDA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Ecut = CutE(BFF(1),rj(1),vnorm(1),EphMFF(1),Kappa(1))   !Cutoff energy for Q(E) in the FF
    !Kinj = NormK(Ecut)                                      !Normalization constant for Q(E)
    ! do i2 = 1, nIn
    !     !Computation of N(E) for each segment
    !     call CalcNE(Kinj,EFF,Ecut,NE)
    !     !call CalcNE_analytical(EFF, z, BFF, rj, vnorm, EphMFF, Kappa, NE)
    !     NEtot = 0.d0    ! Initialization
    !     do j = 1, Ns

    !         write(20,*) EFF, NE(j), j   !E and N(E) for each jet cut
    !         NEtot = NEtot + NE(j)   !Sum for all the jet segments
    !         NEE(j) = NEE(j) + NE(j) * EFF * (Ein-1d0)

    !         !Contribution of the segments up to given lengths
    !         if (z(j) <= 0.01d0*zmax) NE1 = NEtot
    !         if (z(j) <= 0.1d0*zmax) NE2 = NEtot
    !     enddo
    !     write(30,*) EFF/eV, EFF**2*NE1, EFF**2*NE2, EFF**2*NEtot
    !     EFF = EFF*Ein   !Exponential energy increase
    ! enddo

    ! dummy = 0d0
    ! do j=1,Ns
    !     do i2=1,nIn
    !         dummy = dummy + dQ_Edz(Evec(i2), z(j), BFF(j), rj(j), vnorm(j))*deltaZ(j)*Evec(i2)**2*(Ein-1d0)
    !     enddo
    ! enddo
    !write(*,*) dummy
    if (.true.) then
        NEmatrix = 0d0

        ! !$OMP PARALLEL SHARED(NEmatrix, z) DEFAULT(Private)
        ! !$OMP DO
        do i2 = 1, nIn
            write(*,*) i2
            call CalcNE_analytical(Evec(i2), i2, NEmatrix)
            EFF = EFF*Ein   !Exponential energy increase
        enddo
        ! !$OMP END DO
        ! !$OMP END PARALLEL

        EFF = EminFF                        !Initialization
        do i2 = 1, nIn
            NEtot = 0.d0    ! Initialization
            do izz = 1, Ns

                write(20,*) EFF, NEmatrix(i2, izz), izz   !E and N(E) for each jet cut
                NEtot = NEtot + NEmatrix(i2, izz)   !Sum for all the jet segments
                NEE(izz) = NEE(izz) + NEmatrix(i2, izz) * EFF * (Ein-1d0)

                !Contribution of the segments up to given lengths
                if (z(izz) <= 0.01d0*zmax) NE1 = NEtot
                if (z(izz) <= 0.1d0*zmax) NE2 = NEtot
            enddo
            write(30,*) EFF/eV, EFF**2*NE1, EFF**2*NE2, EFF**2*NEtot
            EFF = EFF*Ein   !Exponential energy increase
        enddo

    else
  
        open(40, file="NEmat.dat")
        open(41, file="NEmat2.dat")
        NEmat = 0d0
        do izz = 2, Ns
            call CalcNE_analytical2(izz, Evec, NEmat)
            ! if (izz < 50 .and. z(izz) > zMin(1000*mec2, z(50))) then
            !     write(40, *) z(izz)/z0, factorIntegrante(z(izz)/z(50), z(50)/z0, 1000d0)
            ! endif
            ! if (izz < 49 .and. z(izz) > zMin(Echar(z(49), 1000*mec2, z(50)), z(49))) then
            !     write(41, *) z(izz)/z0, factorIntegrante(z(izz)/z(49), z(49)/z0  , Echar(z(49), 1000d0*mec2, z(50))/mec2) * &
            !                             factorIntegrante(z(49)/z(50), z(50)/z0, 1000d0)
            ! endif
        enddo
        close(40)
        close(41)

        EFF = EminFF                        !Initialization
        do i2 = 1, nIn
            NEtot = 0.d0    ! Initialization
            do izz = 1, Ns

                write(20,*) EFF, NEmat(izz, i2), izz   !E and N(E) for each jet cut
                NEtot = NEtot + NEmat(izz, i2)   !Sum for all the jet segments
                NEE(izz) = NEE(izz) + NEmat(izz, i2) * EFF * (Ein-1d0)

                !Contribution of the segments up to given lengths
                if (z(izz) <= 0.01d0*zmax) NE1 = NEtot
                if (z(izz) <= 0.1d0*zmax) NE2 = NEtot
            enddo
            write(30,*) EFF/eV, EFF**2*NE1, EFF**2*NE2, EFF**2*NEtot
            EFF = EFF*Ein   !Exponential energy increase
        enddo
    endif

    close(20)
    close(30)
    
    end subroutine Edist
   !============================================================================
   ! Edist (END)
   !============================================================================


   !============================================================================
   ! CalcNE (START)
   !============================================================================
    subroutine CalcNE(Kinj,EFF,Ecut,NE)
   !============================================================================
   ! Calculate the numerical particle distribution by evolving each cut as a 
   ! 1-zone
   !============================================================================
    real(dp), intent(in) :: Kinj,EFF,Ecut
    real(dp), intent(out) :: NE(Ns)
    real(dp) :: tlossFF,E0i,E0f,E0(Ns),backD,dtpast,prod,BB,ee,zz,de,zm, dttFF
    integer :: NN(Ns), jj, izz
    !Initialization of the current energy distribution
    NE = 0d0  !Applies to the whole vector

    !Loop for the injection
    do i = 1, Ns

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
            dtpast = Setdt(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0(i))
            
            !Increase in the energy
            E0(i) = E0(i) - TLosses(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0(i))*dtpast
            ! E0(i) = Echar(z(n), z(i), )

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
        NE(i) = Q_E(Kinj,Ecut,E0(i)) * prod * min(dtFF, tlossFF)
        
        !Correction in N(E) for the first cut when E0 is close to Emax
        if (i == 1 .and. EmaxFF < 2d0*E0(i)) NE(i) = NE(i) * (EmaxFF-E0(i)) * tlossFF/E0(i)     
        
        !if (tlossFF < dtFF .and. E0(i) <= EmaxFF) print *, i,EFF/eV,E0(i)/eV,tlossFF
        NN(i) = n
    enddo

    ! if (EFF > 1*mec2) then
    !     do i=10,Ns
    !         zz = z(i)
    !         ee = EFF
    !         zm = max(zMin(EFF, zz), z0)
    !         do while(zz > zm)
    !             BB = B0/gj * (zz/z0)**(-exp_B)
    !             call Setdt(BB,rj(n),vnorm(n),EphMFF(n),Kappa(n),ee,dtpast)
    !             !dtpast = 0.1 * dabs(ee / TLosses(BB,rj(n),vnorm(n),EphMFF(n),Kappa(n),ee))
    !             !Increase in the energy
    !             de = -TLosses(BB,zz,vj,EphMFF(n),Kappa(n),ee) * dtpast
    !             ee = ee + de
    !             zz = zz - vj * dtpast * gj
    !         enddo
    !             !write(*,*) EFF, i, z(i)/z0, zm/z0, ee/EFF, Echar(zm, EFF, z(i))/EFF
    !             !write(*,*) EFF/mec2, zm/z(i), prod, factorIntegrante(zm/z(i), z(i)/z0, EFF/mec2)
    !     enddo
    ! endif

    return
    end subroutine CalcNE
    !============================================================================
    ! CalcNE (END)
    !============================================================================
  
    
    !============================================================================
    ! CalcNE_analytical (START)
    !============================================================================
    subroutine CalcNE_analytical(EFF, i2, NEmatrix)
    !============================================================================
    ! Calculate the electron energy distribution by analytical methods
    !============================================================================
        integer, intent(in) :: i2
        real(dp), intent(in) :: EFF
        real(dp), intent(inout) :: NEmatrix(nIn, Ns)
        real(dp) :: tlossFF, sum, ee, zz, zm, dttFF, BB, rr, vlat
        integer :: jj, izz
        !Initialization of the current energy distribution

        do izz = 1, Ns
            zm = zMin(EFF, z(izz))
            zz = min(z(izz), zmax_inj)
            jj = max(min(izz, Ninj) - 1, 1)
            sum = 0d0
            do while (zz >= max(z0, zm))
                !print*,'jj=',jj
                if (zz < z(jj)) jj = max(jj - 1, 1)
                ee = Echar(zz, EFF, z(izz))
                BB = BFF_fun(zz)
                rr = r_fun(zz)
                vlat = v_lat(zz)
                tlossFF = ee/abs(TLosses(BB, rr, vlat, 1d0, 1d0, ee))
                dttFF = max(min(0.01 * tlossFF, deltatFF(jj)), 1d-10)
                sum = sum + dQ_Edz(ee, zz, BB, rr, vlat) * dttFF * factorIntegrante(zz/z(izz), z(izz)/z0, EFF/mec2)
                !sum = sum + dQ_Edz(ee, zz, BB, rr, vlat) * dttFF * factorIntegrante2(z(izz)/z0, zz/z0, EFF/mec2)
                zz = zz - vj*gj*dttFF
            enddo
            NEmatrix(i2, izz) = sum * deltaZ(izz)
        enddo
   
    end subroutine CalcNE_analytical
    !============================================================================
    ! CalcNE_analytical (END)
    !============================================================================


    !*************************************************************************************************************
    !   Numerical particle distribution by analytical methods
    !*************************************************************************************************************
    !subroutine CalcNE_analytical(EFF,z,BFF,rj,vnorm,EphMFF,Kappa,NE)
    subroutine CalcNE_analytical2(izz, Evec, NEmat)
        real(dp), intent(in) :: Evec(nIn)
        integer, intent(in) :: izz
        real(dp), intent(inout) :: NEmat(Ns, nIn)
        real(dp) :: NEprev_vec(nIn)
        real(dp) :: tlossFF, sum, ee, zp, zm, dttFF, BB, rr, vlat, ee_prev, Ne_prev
        integer :: jj, count1, count2
        !Initialization of the current energy distribution

        do i2=1, nIn
            NEprev_vec(i2) = NEmat(izz-1, i2) / deltaZ(izz-1)
        enddo
        write(*,*) izz

        count1 = 0
        count2 = 0
        do i2=1, nIn
            zm = zMin(Evec(i2), z(izz))
            zp = z(izz-1)
            if (zm < zp) then
                ee_prev = Echar(zp, Evec(i2), z(izz))
                Ne_prev = Interpol1d(ee_prev, Evec, NEprev_vec)
            else
                Ne_prev = 0d0
                zp = zm
            endif
            sum = 0d0
            !Ne_prev =0d0
            ! zp = max(zm, z0)
            ! do while (zp < min(z(izz-1), zmax_inj))
            !     ee = Echar(zp, Evec(i2), z(izz))
            !     BB = BFF_fun(zp)
            !     rr = r_fun(zp)
            !     vlat = v_lat(zp)
            !     tlossFF = ee/dabs(TLosses(BB, rr, vlat, 1d0, 1d0, ee))
            !     dttFF = min(zp, 1.00001d0*(min(z(izz), zmax_inj)-zp)) / (vj*gj)
            !     !dttFF = max(min(0.1 * tlossFF, dttFF), 1d-10)
            !     sum = sum + dQ_Edz(ee, zp, BB, rr, vlat) * dttFF * factorIntegrante(zp/z(izz), z(izz)/z0, Evec(i2)/mec2)
            !     zp = zp + vj*gj*dttFF
            ! enddo

            do while (zp < min(z(izz), zmax_inj))
                ee = Echar(zp, Evec(i2), z(izz))
                BB = BFF_fun(zp)
                rr = r_fun(zp)
                vlat = v_lat(zp)
                tlossFF = ee/abs(TLosses(BB, rr, vlat, 1d0, 1d0, ee))
                dttFF = min(zp, (min(z(izz), zmax_inj)-zp)) / (vj*gj)
                dttFF = max(min(0.1 * tlossFF, dttFF), 1d-10)
                sum = sum + dQ_Edz(ee, zp, BB, rr, vlat) * dttFF * factorIntegrante(zp/z(izz), z(izz)/z0, Evec(i2)/mec2)
                !sum = sum + dQ_Edz(ee, zp, BB, rr, vlat) * dttFF * factorIntegrante2(z(izz)/z0, zp/z0, Evec(i2)/mec2)
                zp = zp + vj*gj*dttFF
            enddo
           ! write(*,*) factorIntegrante(z(izz-1)/z(izz), z(izz)/z0, Evec(i2)/mec2) * factorIntegrante(1.1*max(z0, zMin(ee_prev, z(izz-1)))/z(izz-1), z(izz-1)/z0, ee_prev/mec2), &
            !factorIntegrante(max(z0, zMin(ee_prev, z(izz-1)))*1.1/z(izz), z(izz)/z0, Evec(i2)/mec2)
            NEmat(izz, i2) = deltaZ(izz) * (sum + Ne_prev * factorIntegrante(z(izz-1)/z(izz), z(izz)/z0, Evec(i2)/mec2))
            !NEmat(izz, i2) = deltaZ(izz) * (sum + Ne_prev * factorIntegrante2(z(izz)/z0, z(izz-1)/z0, Evec(i2)/mec2))
        enddo

        return
    end subroutine CalcNE_analytical2

    !*************************************************************************************************************
    !   Setting of an appropriate time step for the N(E) computation
    !*************************************************************************************************************
    real(dp) function Setdt(B,r_j,v_norm,Eph,Kapp,E)
    real(dp), intent(in) :: B,r_j,v_norm,Eph,Kapp,E
    real(dp) :: tloss
    !Characteristic time for losses
    tloss = abs(E/TLosses(B,r_j,v_norm,Eph,Kapp,E))
    !Set dt according to the characterstic losses and injection times
    Setdt = 0.1d0 * min(tloss,dtFF)
    return
    end function
    

    !************************************************************************************
    !   Determination of the cutoff energy of the particle injection
    !************************************************************************************
    real(dp) function CutE (BFF1,rj1,vnorm1,Eph1,Kappa1)
    real(dp), intent(in) :: BFF1,rj1,vnorm1,Eph1,Kappa1
    real(dp) :: Emax1,Emax2,Eacc,lE1,lE2,lEmid,error,check
    !Injected power (erg/s)
    Eacc = facc*qe*BFF1*c
    !Initialization
    lE1 = log10(1d9*eV)  ! log10(Emin [CGS]) ~ 1 GeV
    lE2 = log10(1d14*eV) ! log10(Emax [CGS]) ~ 100 TeV
    error = 1.d-4        ! error < 10^(-4)
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


    !************************************************************************************
    !   Injection function (only one injection cell)
    !************************************************************************************
    real(dp) function Q_E (K,Ecut,E)
    real(dp), intent(in) :: K,Ecut,E
    !Q(E) = 0 if outside the injection energy limits
    if (E <= EmaxFF .and. E >= EminFF) then
        Q_E = K*E**(-p) * exp(-E/Ecut)
    else
        Q_E = 0d0
    endif
    end function Q_E
    

    !************************************************************************************
    !   Injection function (for multiple injection cells?)
    !************************************************************************************
    real(dp) function Q_Ez (E, iz)
    real(dp), intent(in) :: E
    integer, intent(in) :: iz
    real(dp) :: Ecut, Kinj
    Ecut = CutE(BFF(iz),rj(iz),vnorm(iz),EphMFF(iz),Kappa(iz))!Cutoff energy for Q(E) in the FF
    Kinj = NormKz(Ecut)                                       !Normalization constant for Q(E)
    !Q(E) = 0 if outside the injection energy limits
    if (E <= EmaxFF .and. E >= EminFF) then
        Q_Ez = Kinj*E**(-p) * exp(-E/Ecut)
    else
        Q_Ez = 0d0
    endif
    end function Q_Ez
  

    !************************************************************************************
    !   Injection function (for multiple injection cells)
    !************************************************************************************
    real(dp) function dQ_Edz(E, zz, BB, rr, vlat)
    real(dp), intent(in) :: E, zz, BB, rr, vlat
    real(dp) :: Ecut, Kinj
    Ecut = CutE(BB, rr, vlat, 1d0, 1d0)   !Cutoff energy for Q(E) in the FF
    Kinj = NormKz(Ecut) / (zz * log_zmax_zinj_z0)         !Normalization constant for Q(E)
    !Q(E) = 0 if outside the injection energy limits
    if (E <= EmaxFF .and. E >= EminFF) then
        dQ_Edz = Kinj*E**(-p) * exp(-E/Ecut)
    else
        dQ_Edz = 0d0
    endif
    end function dQ_Edz


    !************************************************************************************
    !   Injection normalization (only one injection cell)
    !************************************************************************************  
    real(dp) function NormK(Ecut)
    real(dp), intent(in) :: Ecut
    real(dp) :: E, Einc, dE, QInt, f, f1, f2, E2
    integer, parameter :: N = 300   !Number of subintervals for the integral
    Einc = (EmaxFF/EminFF)**(1d0/(N-1))  !Increase factor of the energy
    !Initialization
    E = EminFF
    QInt = 0.d0
    ! The integral is solved with the trapezoid method
    f1 = E * Q_E(1.d0,Ecut,E)
    do i = 1, N-1 
        dE = E*(Einc-1d0)                    !Energy interval
        E2 = E*Einc
        f2 = E2 * Q_E(1.d0,Ecut,E2)
        f = (f1 + f2)/2.d0
        QInt = QInt + f*dE
        f1 = f2
        E = E2
    end do
    NormK = LnTFF/QInt    !Normalization constant
    return
    end function NormK


    !************************************************************************************
    !   Injection normalization
    !************************************************************************************  
    real(dp) function NormKz(Ecut)
    real(dp), intent(in) :: Ecut
    real(dp) :: E, Einc, dE, QInt, f, f1, f2, E2
    integer, parameter :: N = 300   !Number of subintervals for the integral
    Einc = (EmaxFF/EminFF)**(1d0/(N-1))  !Increase factor of the energy
    !Initialization
    E = EminFF
    QInt = 0.d0
    ! The integral is solved with the trapezoid method
    f1 = E * Q_E(1.d0,Ecut,E)
    do i = 1, N-1 
        dE = E*(Einc-1d0)                    !Energy interval
        E2 = E*Einc
        f2 = E2 * Q_E(1.d0,Ecut,E2)
        f = (f1 + f2)/2.d0
        QInt = QInt + f*dE
        f1 = f2
        E = E2
    end do
    NormKz = LnTFF/QInt   !Normalization constant

    return
    end function NormKz


    !============================================================================
    ! Calc_Ue (START)
    !============================================================================
    subroutine Calc_Ue()
    !============================================================================
    ! Computation of the electron energy density in each segment (FF)
    !============================================================================
        real(dp) :: Ue, Ne(nIn, Ns), Ee(nIn), Ee_step, dEe, E_tot

        ! Read from a file: Ee(i), Ne(i,l)
        open(70,file='NE_cuts.dat')
        do i = 1, nIn
            do l = 1, Ns
                read(70,*) Ee(i), Ne(i,l), dummy
            enddo
        enddo
        close(70)
        ! Integrate the particle energy distribution of each segment
        Ee_step = Ee(2)/Ee(1)
        open(28, file='U_e.dat')
        do l = 1, Ns
            E_tot = 0.d0
            do i = 1, nIn
                dEe = Ee(i)*(Ee_step-1.d0)
                E_tot = E_tot + Ee(i)*Ne(i,l) * dEe
            enddo
            Ue = E_tot/VFF(l) 
            write(28,*) l, z(l), Ue
        enddo
        close(28)

        return
    end subroutine Calc_Ue
    !============================================================================
    ! Calc_Ue (END)
    !============================================================================


    !*********************************************************************************
    !  Minimum z for a given characteristic curve ending at (zf, Ef)
    !*********************************************************************************
    real(dp) function zMin(E, z_)
    real(dp), intent(in) :: E, z_
    real(dp) :: A, B, D, xi_f, g_f, f, U0

    U0 = (B0/gj)**2 / (8.d0*pi)
    A = 4.d0/3.d0 * sigma_T * (U0/gj) * z0 / mec2
    B = 2.d0/3.d0 * geomIndex
    !B = 2.d0/3.d0 / thetaj
    xi_f = z_/z0
    g_f = E/mec2
    f = 2*exp_B
    D = B + f - 1.d0

    zMin = (1d0 + D/(A*g_f))**(-1d0/D) * z_

    return
    end function


    !*********************************************************************************
    !  Analytical function for the characteristic curve
    !*********************************************************************************
    real(dp) function Echar(z_, Efinal, zfinal)
    real(dp), intent(in) :: z_, Efinal, zfinal
    real(dp) :: A, B, C1, D, xi_final, g_final, xi, f, U0

    U0 = (B0/gj)**2 / (8.d0*pi)
    A = 4.d0/3.d0 * sigma_T * (U0/gj) * z0 / mec2
    B = 2.d0/3.d0 * geomIndex
    !B = 2.d0/3.d0 / thetaj
    xi_final = zfinal/z0
    g_final = Efinal/mec2
    xi = z_/z0
    f = 2*exp_B
    D = B + f - 1.d0
    C1 = xi_final**(-(B+f)) * (xi_final**f * D + A * xi_final * g_final) / (D*g_final)

    Echar = mec2 * (D*xi**f / (C1*D*xi**(B+f) - A*xi))
  
    return
    end function Echar


    !*********************************************************************************
    !  Integrating factor 
    !*********************************************************************************
    real(dp) function factorIntegrante(xp, xi_f, g_f)
    real(dp), intent(in) :: xp, g_f, xi_f
    real(dp) :: A, B, aux, D, f, U0, C1

    !write(*, *) zp/z0, zf/z0, Ef
    U0 = (B0/gj)**2 / (8.d0*pi)
    A = 4.d0/3.d0 * sigma_T * (U0/gj) * z0 / mec2
    B = 2.d0/3.d0 * geomIndex
    f = 2d0*exp_B
    D = B + f - 1.d0
    
    C1 = xi_f**(-(B+f)) * (xi_f**f*D + A*xi_f*g_f) / g_f
    factorIntegrante = xp**(-B) * ( (C1-A*xi_f**(-D)) / (C1-A*xi_f**(-D)*xp**(-D)) )**2

    return
    end function

    !----------------------------------------------------------------------------------

    real(dp) function factorIntegrante2(xm, xp, ym)
    real(dp), intent(in) :: xm, ym, xp
    real(dp) :: A, B, aux, D, f, U0, C1, g

    U0 = (B0/gj)**2 / (8.d0*pi)
    A = 4.d0/3.d0 * sigma_T * (U0/gj) * z0 / mec2
    B = 2.d0/3.d0 * geomIndex
    f = 2d0*exp_B
    D = B + f - 1.d0

    C1 = (xm**(-B - f)*(-xm**f + B*xm**f + f*xm**f + A*xm*ym))/(D*ym)
    D = B - 1d0 + f
    g = B + f
    factorIntegrante2 = (C1*D - A*xm**(1d0-g)) / (C1*D - A*xp**(1d0-g))

    return
    end function

!=======================================================================================

END MODULE edistribution
