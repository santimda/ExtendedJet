MODULE edistribution
    use constants
    use globals
    use losses
    implicit none

contains
    
   !============================================================================
   ! Edist (START)
   !============================================================================
    subroutine Edist(particle)
   !============================================================================
   ! This subroutine calculates the electron energy distribution at each segment
   !============================================================================
    character(1), intent(in) :: particle
    real(dp) :: Kinj, Ecut, Ein, EFF, NE(Ns), NEtot, NE1, NE2, EmaxFF2, EminFF2, restEnergy, Lsy, dEdt
    real(dp) :: NEmatrix(nIn, Ns), NEmat(Ns, nIn), Evec(nIn)
    real(dp) :: NEE(Ns), Nparticles, flux, Eparticles, energyFlux, vol, fluxp, eFluxp, Ljett, magnetization
    integer :: izz

    if (particle == 'e') then
        open(11,file='Emax_FF.dat')
        open(20,file='NE_cuts.dat')
        open(30,file='NE_jet.dat')
        EminFF2 = EminFF
        EmaxFF2 = EmaxFF
        restEnergy = mec2
    else if (particle == 'p') then
        open(20,file='Np_cuts.dat')
        open(30,file='Np_jet.dat')
        open(11, file='Epmax_FF.dat')
        EminFF2 = EpminFF
        EmaxFF2 = EpmaxFF
        restEnergy = mpc2
    else
        write(*,*) "FATAL ERROR: WRONG PARTICLE CLASS!!!"
    endif

    do i = 1,Ns
        !Quantities at each segment
        Ecut = CutE(BFF(i),rj(i),vnorm(i),EphMFF(i),Kappa(i),particle)   !Cutoff energy for Q(E) in the FF
        if (Ecut > EmaxFF2/3.d0) then
            write(*,*) 'EmaxFF in globals.f90 is too low!'
            stop
        else
            write(11,*) i, z(i), Ecut/eV
        endif
    enddo
    close(11)
    
    ! Calculate logarithmically-spaced energy vector
    call vector_log(EminFF2, EmaxFF2, Ein, Evec ) 

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

    if (.true.) then ! FINITE DIFFERENCES
        NEmatrix = 0d0
        call calcNE_paramo(Evec, NEmatrix)
        !call CalcNE_fds(Evec, NEmatrix, 150, 300)
    endif
    
    if (.true.) then

        if (.false.) then ! QUASI-ANALYTIC SOLUTION
            NEmatrix = 0d0

        ! !$OMP PARALLEL SHARED(NEmatrix, z) DEFAULT(Private)
        ! !$OMP DO
            do i2 = 1, nIn
                write(*, *) i2
                call CalcNE_analytical(Evec(i2), i2, NEmatrix, particle)
                EFF = EFF*Ein   !Exponential energy increase
            enddo
        endif

        flux = 0d0
        energyFlux = 0d0

        ! PLOT SOME QUANTITIES
        open(50, file='Quantities.dat')
        do izz = 1, Ns
            Nparticles = 0d0
            Eparticles = 0d0
            Lsy = 0d0
            do i2 = 1, nIn
                Nparticles = Nparticles + NEmatrix(i2, izz) * Evec(i2) * (sqrt(Ein) - 1d0/sqrt(Ein))
                Eparticles = Eparticles + NEmatrix(i2, izz) * Evec(i2)**2 * (sqrt(Ein)-1d0/sqrt(Ein))
                dEdt = 4d0/3d0 * sigma_T * c * (BFF_fun(z(izz))**2d0/(8d0*pi)) * (Evec(i2)/mec2)**2d0
                Lsy = Lsy + NEmatrix(i2, izz) * dEdt * Evec(i2) * (sqrt(Ein) - 1d0/sqrt(Ein))
            enddo
            vol = pi * r_fun(z(izz))**2 * deltaZ(izz)
            Ljett = pi*r_fun(z(izz))**2*lorentzFactor(z(izz))**2*betaJet(z(izz))*c * (Eparticles/(vol*lorentzFactor(z(izz))) &
             + BFF_fun(z(izz))**2/(8d0*pi) + Lsy/(4d0*r_fun(z(izz))**2d0*c))
            magnetization = BFF_fun(z(izz))**2/(8d0*pi*Eparticles/vol)
            write(50, *) z(izz)/R_g, Eparticles/(vol*lorentzFactor(z(izz))), BFF_fun(z(izz))**2/(8d0*pi), Ljett/Lj, magnetization, &
                            (Ljett + Lsy)/Lj
            fluxp = flux
            eFluxp = energyFlux
            flux = Nparticles / vol * lorentzFactor(z(izz)) * betaJet(z(izz)) * c * pi * r_fun(z(izz))**2
            energyFlux = Eparticles / vol * betaJet(z(izz)) * c * pi * r_fun(z(izz))**2
            write(*, *) "z = ", z(izz)/R_g, "Flux = ", flux, "EnergyFlux = ", energyFlux + Lsy, "RatioPart = ", & 
                                fluxp/flux, "RatioE = ", eFluxp/energyFlux
        enddo
        close(50)
        ! !$OMP END DO
        ! !$OMP END PARALLEL

        EFF = EminFF2                        !Initialization
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
            call CalcNE_analytical2(izz, Evec, NEmat,particle)
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
!     subroutine CalcNE(Kinj,EFF,Ecut,NE,particle)
!    !============================================================================
!    ! Calculate the numerical particle distribution by evolving each cut as a 
!    ! 1-zone
!    !============================================================================
!     real(dp), intent(in) :: Kinj,EFF,Ecut
!     real(dp), intent(out) :: NE(Ns)
!     character(1), intent(in) :: particle
!     real(dp) :: tlossFF,E0i,E0f,E0(Ns),backD,dtpast,prod,BB,ee,zz,de,zm, dttFF
!     integer :: NN(Ns), jj, izz
!     !Initialization of the current energy distribution
!     NE = 0d0  !Applies to the whole vector

!     !Loop for the injection
!     do i = 1, Ns

!         !Initial value of the energy for each injection
!         E0(i) = EFF

!         !Initial value of the lower energy limit for the cut-by-cut one-zone computation
!         E0i = EFF
        
!         !Initialization
!         backD = 0d0 !Backwards distance
!         prod = 1d0
!         j = i       !Reference for the product
!         n = i       !Reference segment for B,nf,rj...

!         !Loop for the distance that each population has made
!         do while (backD <= z(i)-z(1))
        
!             !Null contribution to N(E) if E0 > Emax (not physically possible)
!             if (E0(i) > EmaxFF) exit

!             !Change the reference segment when backD is long enough
!             if (backD > (z(i)-z(n))) n = n - 1
            
!             !Selection of the adequate time steps for each position (dtpast already in the FF)
!             dtpast = Setdt(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0(i),particle)
            
!             !Increase in the energy
!             E0(i) = E0(i) - TLosses(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0(i),particle)*dtpast
!             ! E0(i) = Echar(z(n), z(i), )

!             !Condition to find the point where the segment i was at the current position of segment j
!             if (backD > z(i)-z(j)) then

!                 !Upper energy limit at the segment j position
!                 E0f = E0(i)

!                 !Product of the quotient of the losses
!                 prod = prod * TLosses(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0f,particle)/        &
!                 &             TLosses(BFF(n),rj(n),vnorm(n),EphMFF(n),Kappa(n),E0i,particle)

!                 !Lower energy limit at the segment j-1 position
!                 E0i = E0(i)

!                 !Decrease in the segment count
!                 j = j - 1

!             endif

!             !Increase in the backwards distance (LAB frame)
!             backD = backD + vj*dtpast*gj

!         enddo
        
!         !Cooling time for energy E0 and segment i
!         tlossFF = E0(i)/abs(TLosses(BFF(i),rj(i),vnorm(i),EphMFF(i),Kappa(i),E0(i),particle))

!         !Energy distribution for each segment, corrected when tloss' < dt' for a given energy (to avoid 
!         !overestimation at high energies due to the 1st segment being too long)
!         NE(i) = Q_E(Kinj,Ecut,E0(i)) * prod * min(dtFF, tlossFF)
        
!         !Correction in N(E) for the first cut when E0 is close to Emax
!         if (i == 1 .and. EmaxFF < 2d0*E0(i)) NE(i) = NE(i) * (EmaxFF-E0(i)) * tlossFF/E0(i)     
        
!         !if (tlossFF < dtFF .and. E0(i) <= EmaxFF) print *, i,EFF/eV,E0(i)/eV,tlossFF
!         NN(i) = n
!     enddo

!     ! if (EFF > 1*mec2) then
!     !     do i=10,Ns
!     !         zz = z(i)
!     !         ee = EFF
!     !         zm = max(zMin(EFF, zz), z0)
!     !         do while(zz > zm)
!     !             BB = B0/gj * (zz/z0)**(-exp_B)
!     !             call Setdt(BB,rj(n),vnorm(n),EphMFF(n),Kappa(n),ee,dtpast)
!     !             !dtpast = 0.1 * dabs(ee / TLosses(BB,rj(n),vnorm(n),EphMFF(n),Kappa(n),ee))
!     !             !Increase in the energy
!     !             de = -TLosses(BB,zz,vj,EphMFF(n),Kappa(n),ee) * dtpast
!     !             ee = ee + de
!     !             zz = zz - vj * dtpast * gj
!     !         enddo
!     !             !write(*,*) EFF, i, z(i)/z0, zm/z0, ee/EFF, Echar(zm, EFF, z(i))/EFF
!     !             !write(*,*) EFF/mec2, zm/z(i), prod, factorIntegrante(zm/z(i), z(i)/z0, EFF/mec2)
!     !     enddo
!     ! endif

!     return
!     end subroutine CalcNE
    !============================================================================
    ! CalcNE (END)
    !============================================================================
  
    
    !============================================================================
    ! CalcNE_analytical (START)
    !============================================================================
    subroutine CalcNE_analytical(EFF, i3, NEmatrix, particle)
    !============================================================================
    ! Calculate the electron energy distribution by analytical methods
    !============================================================================
        use globals
        integer, intent(in) :: i3
        character(1), intent(in) :: particle
        real(dp), intent(in) :: EFF
        real(dp), intent(inout) :: NEmatrix(nIn, Ns)
        real(dp) :: tlossFF, sum, ee, zz, zm, dttFF, BB, rr, vlat, gammap, restEnergy
        real(dp) :: bgj_, one_gjbz
        integer :: jj, izz
        !Initialization of the current energy distribution

        if ( particle == 'e') then
            restEnergy = mec2
        else if ( particle == 'p') then
            restEnergy = mpc2
        endif

        gammap = EFF / restEnergy
        do izz = 1, Ns
            zm = zMin(EFF, z(izz), particle)
            zz = min(z(izz), zmax_inj)
            jj = max(min(izz, Ninj) - 1, 1)
            sum = 0d0
            !one_gjbz = 1d0 / betaGamma(z(izz))
            one_gjbz = 1d0 / sqrt(lorentzFactor(z(izz))**2 - 1d0)
            !one_gjbz = 1d0 / (sqrt(gjz**2 - 1d0))
            !do while (zz >= max(z0, zm))
            ee = mec2
            do while (zz >= z0 .and. ee < EmaxFF)
                !print*,'jj=',jj
                if (zz < z(jj)) jj = max(jj - 1, 1)
                ee = Echar(zz, EFF, z(izz), particle)
                BB = BFF_fun(zz)
                rr = r_fun(zz)
                vlat = v_lat(zz)
                tlossFF = ee/abs(TLosses(BB, rr, vlat, 1d0, 1d0, ee, particle))
                dttFF = max(min(0.1 * tlossFF, deltatFF(jj)), 1d-15*zz/(betaGamma(zz)*c))
                !bgj_ = betaGamma(zz)
                bgj_ = sqrt(lorentzFactor(zz)**2 - 1d0)
                !sum = sum + dQ_Edz(ee, zz, BB, rr, vlat, particle) * dttFF * (bgj_ * one_gjbz) * &
                !        factorIntegrante(zz/z(izz), z(izz)/z0, gammap, particle)
                sum = sum + dQ_Edz(ee, zz, BB, rr, vlat, particle) * dttFF * (bgj_ * one_gjbz) * & 
                        factorIntegrante(zz/z(izz), z(izz)/z0, gammap, particle)
                zz = zz - bgj_*c*dttFF
                !write(*, *) factorIntegrante(zz/z(izz), z(izz)/z0, gammap, particle)
            enddo
            NEmatrix(i3, izz) = sum * deltaZ(izz)
        enddo
   
    end subroutine CalcNE_analytical
    !============================================================================
    ! CalcNE_analytical (END)
    !============================================================================


    !*************************************************************************************************************
    !   Numerical particle distribution by analytical methods
    !*************************************************************************************************************
    !subroutine CalcNE_analytical(EFF,z,BFF,rj,vnorm,EphMFF,Kappa,NE)
    subroutine CalcNE_analytical2(izz, Evec, NEmat, particle)
        real(dp), intent(in) :: Evec(nIn)
        character(1), intent(in) :: particle
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
            zm = zMin(Evec(i2), z(izz), particle)
            zp = z(izz-1)
            if (zm < zp) then
                ee_prev = Echar(zp, Evec(i2), z(izz), particle)
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
                ee = Echar(zp, Evec(i2), z(izz), particle)
                BB = BFF_fun(zp)
                rr = r_fun(zp)
                vlat = v_lat(zp)
                tlossFF = ee/abs(TLosses(BB, rr, vlat, 1d0, 1d0, ee, particle))
                dttFF = min(zp, (min(z(izz), zmax_inj)-zp)) / (betaGamma(zp)*c)
                dttFF = max(min(0.1 * tlossFF, dttFF), 1d-10)
                sum = sum + dQ_Edz(ee, zp, BB, rr, vlat, particle) * dttFF * factorIntegrante(zp/z(izz), z(izz)/z0, Evec(i2)/mec2, &
                            particle) * betaGamma(zp) / betaGamma(z(izz))
                !sum = sum + dQ_Edz(ee, zp, BB, rr, vlat) * dttFF * factorIntegrante2(z(izz)/z0, zp/z0, Evec(i2)/mec2)
                zp = zp + betaGamma(zp)*c*dttFF
            enddo
           ! write(*,*) factorIntegrante(z(izz-1)/z(izz), z(izz)/z0, Evec(i2)/mec2) * factorIntegrante(1.1*max(z0, zMin(ee_prev, z(izz-1)))/z(izz-1), z(izz-1)/z0, ee_prev/mec2), &
            !factorIntegrante(max(z0, zMin(ee_prev, z(izz-1)))*1.1/z(izz), z(izz)/z0, Evec(i2)/mec2)
            NEmat(izz, i2) = deltaZ(izz) * (sum + Ne_prev * factorIntegrante(z(izz-1)/z(izz), z(izz)/z0, Evec(i2)/mec2, particle))
            !NEmat(izz, i2) = deltaZ(izz) * (sum + Ne_prev * factorIntegrante2(z(izz)/z0, z(izz-1)/z0, Evec(i2)/mec2))
        enddo

        return
    end subroutine CalcNE_analytical2

    !*************************************************************************************************************
    !   Setting of an appropriate time step for the N(E) computation
    !*************************************************************************************************************
    real(dp) function Setdt(B,r_j,v_norm,Eph,Kapp,E,particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: B,r_j,v_norm,Eph,Kapp,E
    real(dp) :: tloss
    !Characteristic time for losses
    tloss = abs(E/TLosses(B,r_j,v_norm,Eph,Kapp,E,particle))
    !Set dt according to the characterstic losses and injection times
    Setdt = 0.1d0 * min(tloss,dtFF)
    return
    end function
    

    !************************************************************************************
    !   Determination of the cutoff energy of the particle injection
    !************************************************************************************
    real(dp) function CutE (BFF1,rj1,vnorm1,Eph1,Kappa1, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: BFF1,rj1,vnorm1,Eph1,Kappa1
    real(dp) :: Emax1,Emax2,Eacc,lE1,lE2,lEmid,error,check
    !Injected power (erg/s)
    Eacc = facc*qe*BFF1*c
    !Initialization
    lE1 = log10(1d9*eV)  ! log10(Emin [CGS]) ~ 1 GeV
    lE2 = log10(1d14*eV) ! log10(Emax [CGS]) ~ 100 TeV
    if (particle == 'p') lE2 = log10(1d20*eV)
    error = 1.d-4        ! error < 10^(-4)
    !Maximum particle energy
    do while (abs(lE2-lE1) > error)
        lEmid = (lE1+lE2)/2d0
        check =  Eacc - abs(TLosses(BFF1,rj1,vnorm1,Eph1,Kappa1,10d0**lEmid, particle))
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
    real(dp) function Q_Ez (E, iz, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: E
    integer, intent(in) :: iz
    real(dp) :: Ecut, Kinj
    Ecut = CutE(BFF(iz),rj(iz),vnorm(iz),EphMFF(iz),Kappa(iz),particle)!Cutoff energy for Q(E) in the FF
    Kinj = NormKz(Ecut,particle)                                       !Normalization constant for Q(E)
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
    real(dp) function dQ_Edz(E, zz, BB, rr, vlat, particle)
    real(dp), intent(in) :: E, zz, BB, rr, vlat
    real(dp) :: Ecut, Kinj, EminFF2, EmaxFF2, intInj
    character(1), intent(in) :: particle
    Ecut = CutE(BB, rr, vlat, 1d0, 1d0, particle)   !Cutoff energy for Q(E) in the FF
    !Kinj = NormKz(Ecut, particle) / (zz * log(zmax_inj/z0))         !Normalization constant for Q(E)
    intInj = z0 / (a_inj-1d0) * (1d0 - (z0/zmax_inj)**(a_inj-1d0))
    Kinj = NormKz(Ecut, particle) * (z0/zz)**a_inj / intInj
    !Q(E) = 0 if outside the injection energy limits
    if (particle == 'e') then
        EminFF2 = EminFF
        EmaxFF2 = EmaxFF
    else if (particle == 'p') then
        EminFF2 = EpminFF
        EmaxFF2 = EpmaxFF
    endif
    !if (E <= EmaxFF2 .and. E >= EminFF2) then
    if (E >= EminFF2 .and. E <= Ecut) then
        dQ_Edz = Kinj*E**(-p) !* exp(-E/Ecut)
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
    real(dp) function NormKz(Ecut, particle)
    real(dp), intent(in) :: Ecut
    character(1), intent(in) :: particle
    real(dp) :: E, Einc, dE, QInt, f, f1, f2, E2, EminFF2, EmaxFF2, LntFF2
    integer, parameter :: N = 300   !Number of subintervals for the integral
    if (particle == 'e') then
        EminFF2 = EminFF
        EmaxFF2 = EmaxFF
        LnTFF2 = LnTFF
    else if (particle == 'p') then
        EminFF2 = EpminFF
        EmaxFF2 = EpmaxFF
        LnTFF2 = LnTpFF
    endif 
    Einc = (EmaxFF2/EminFF2)**(1d0/(N-1))  !Increase factor of the energy
    !Initialization
    E = EminFF2
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
    NormKz = LnTFF2/QInt   !Normalization constant

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
    !  Minimum z for a given characteristic curve ending at (zf, Ef) ! FUNCIONA MAL
    !*********************************************************************************
    real(dp) function zMin(E, z_, particle)
    real(dp), intent(in) :: E, z_
    character(1), intent(in) :: particle
    real(dp) :: A, B, D, xi_f, g_f, f, U0, restEnergy, factor

    restEnergy = mec2
    factor = 1d0
    if (particle == 'p') then
        factor = (mec2/mpc2)**2
        restEnergy = mpc2
    endif

    xi_f = z_/z0
    g_f = E/restEnergy

    ! U0 = (B0/lorentzFactor(z0))**2 / (8.d0*pi)
    ! A = 4.d0/3.d0 * sigma_T * (U0/betaGamma(z0)) * z0 * factor / restEnergy
    ! B = 2.d0/3.d0 * geomIndex
    ! B = 0d0
    ! f = 2*exp_B*geomIndex + 3d0 * gammaIndex
    ! D = B + f - 1d0

    A = z0 * sigma_T * B0**2d0 / (6d0*pi*mec2*lorentzFactor(z0)**3d0)
    f = 2*exp_B*geomIndex + 3d0*gammaIndex
    B = 2d0/3d0 * geomIndex
    B = 0d0
    D = B + f - 1d0

    zMin = (1d0 + D/(A*g_f))**(-1d0/D) * z_

    return
    end function


    !*********************************************************************************
    !  Analytical function for the characteristic curve
    !*********************************************************************************
    real(dp) function Echar(z_, Efinal, zfinal, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: z_, Efinal, zfinal
    real(dp) :: A, B, C1, D, xi_final, g_final, xi, f, U0, factor, restEnergy

    !U0 = (B0/lorentzFactor(z0))**2 / (8.d0*pi)
    factor = 1d0
    restEnergy = mec2
    if (particle == 'p') then
        factor = (mec2/mpc2)**2
        restEnergy = mpc2
    endif

    xi_final = zfinal/z0
    g_final = Efinal/restEnergy
    xi = z_/z0

    A = z0 * sigma_T * B0**2d0 / (6d0*pi*mec2*lorentzFactor(z0)**3d0)
    f = 2*exp_B*geomIndex + 3d0*gammaIndex
    B = 2d0/3d0 * geomIndex
    D = B + f - 1d0

    C1 = xi_final**(-(B+f)) * (xi_final**f * D + A * xi_final * g_final) / (D*g_final)

    Echar = restEnergy * (D*xi**f / (C1*D*xi**(B+f) - A*xi))
    if ((C1*D*xi**(B+f) - A*xi) <= 0d0) then
        Echar = EmaxFF * 100
    endif
  
    return
    end function Echar

    !*********************************************************************************
    !  Analytical function for the characteristic curve
    !*********************************************************************************
    real(dp) function Echar_sy(z_, Efinal, zfinal, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: z_, Efinal, zfinal
    real(dp) :: A, B, C1, D, xi_final, g_final, xi, f, U0, factor, restEnergy

    A = z0 * sigma_T * B0**2d0 / (6d0*pi*mec2*lorentzFactor(z0)**3d0)
    f = 2*exp_B*geomIndex + 3d0*gammaIndex
    xi_final = zfinal/z0
    g_final = Efinal/mec2
    xi = z_/z0
    D = f - 1d0
    
    Echar_sy = mec2 * g_final / (1d0 - A*g_final/D * (1d0/xi**D - 1d0/xi_final**D))
  
    return
    end function Echar_sy

    real(dp) function zMin_sy(E_, z_, particle)
    real(dp), intent(in) :: E_, z_
    character(1), intent(in) :: particle
    real(dp) :: A, B, D, xi_f, g_f, f, U0, restEnergy, factor

    A = z0 * sigma_T * B0**2d0 / (6d0*pi*mec2*lorentzFactor(z0)**3d0)
    f = 2*exp_B*geomIndex + 3d0*gammaIndex
    xi_f = z_/z0
    g_f = E_/mec2
    D = f - 1d0

    zMin_sy = xi_f * (A*g_f / (xi_f**D * D + A*g_f))**(1d0/D) * z0

    return
    end function zMin_sy

    !*********************************************************************************
    !  Integrating factor 
    !*********************************************************************************
    real(dp) function factorIntegrante(xp, xi_f, g_f, particle)
    real(dp), intent(in) :: xp, g_f, xi_f
    real(dp) :: A, B, aux, D, f, U0, C1, factor, restEnergy
    character(1), intent(in) :: particle

    restEnergy = mec2
    factor = 1d0
    if (particle == 'p') then
        factor = (mec2/mpc2)**2
        restEnergy = mpc2
    endif

    A = z0 * sigma_T * B0**2d0 / (6d0*pi*mec2*lorentzFactor(z0)**3d0)
    f = 2*exp_B*geomIndex + 3d0*gammaIndex
    B = 2d0/3d0 * geomIndex
    D = B + f - 1d0
    
    C1 = xi_f**(-(B+f)) * (xi_f**f*D + A*xi_f*g_f) / g_f
    !factorIntegrante = xp**(-B) * ( (C1-A*xi_f**(-D)) / (C1-A*xi_f**(-D)*xp**(-D)) )**2
    factorIntegrante = ( (C1-A*xi_f**(-D)) / (C1-A*xi_f**(-D)*xp**(-D)) )**2

    return
    end function

    !*********************************************************************************
    !  Integrating factor 
    !*********************************************************************************
    real(dp) function factorIntegrante_sy(xp, xi_f, g_f, particle)
    real(dp), intent(in) :: xp, g_f, xi_f
    real(dp) :: A, B, aux, D, f, U0, C1, factor, restEnergy
    character(1), intent(in) :: particle

    A = z0 * sigma_T * B0**2d0 / (6d0*pi*mec2*lorentzFactor(z0)**3d0)
    f = 2*exp_B*geomIndex + 3d0*gammaIndex
    D = f - 1d0

    factorIntegrante_sy = (D*xi_f**(2*D) / (A*g_f*xi_f**D-xp**D*(A*g_f+D*xi_f**D)) * (xp/xi_f)**D)**(-1d0/(A*g_f))

    return
    end function


    !----------------------------------------------------------------------------------

    real(dp) function factorIntegrante2(xm, xp, ym, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: xm, ym, xp
    real(dp) :: A, B, aux, D, f, U0, C1, g, factor, restEnergy

    restEnergy = mec2
    factor = 1d0
    if (particle == 'p') then
        factor = (mec2/mpc2)**2
        restEnergy = mpc2
    endif
    U0 = (B0/lorentzFactor(z0))**2 / (8.d0*pi)
    A = 4.d0/3.d0 * sigma_T * (U0/betaGamma(z0)) * z0 * factor / restEnergy
    B = 2.d0/3.d0 * geomIndex
    B = 0d0
    f = 2d0*exp_B*geomIndex + 3d0*gammaIndex
    g = B + f
    D = g - 1d0

    C1 = xm**(-g) * (D*xm**f + A*xm*ym) / (D*ym)

    factorIntegrante2 = ( (C1*D - A*xm**(-D)) / (C1*D - A*xp**(-D)) )**2 * (xm/xp)**B

    return
    end function

!=======================================================================================


    subroutine TriDiagSys(aa, bb, cc, dd, naux)
        !---------------------------------------------------------------------------
        ! Solves a system with tridiagonal matrix by LU factorization (diag(L) = 1).
        ! a - lower codiagonal (i=2,n)
        ! b - main diagonal (i=1,n)
        ! c - upper codiagonal (i=1,n-1)
        ! d - constant terms (i=1,n); solution on exit
        ! n - order of system.
        !---------------------------------------------------------------------------
        integer, intent(in) :: naux
        real(dp), intent(inout) :: aa(naux), bb(naux), cc(naux), dd(naux)
        integer :: ii
        if (bb(1) == 0d0) then
            write(*, *) "TriDiagSys: singular matrix!"
            return
        endif

        do i=2, naux
            ! factorization
            aa(i) = aa(i)/ bb(i-1)
            bb(i) = bb(i) - aa(i)*cc(i-1)
            if (bb(i) == 0d0) then
                write(*, *) "TriDiagSys: singular matrix!"
                return
            endif
            dd(i) = dd(i) - aa(i)*dd(i-1)
        enddo
        
        dd(naux) = dd(naux) / bb(naux)
        
        ! backward substitution
        do i=naux-1, 1, -1
            dd(i) = (dd(i) - cc(i)*dd(i+1))/bb(i)
        enddo

        return
    end subroutine TriDiagSys

    real(dp) function Bfun(g, z_phys)
        real(dp), intent(in) :: g, z_phys
        real(dp) :: Ee, dgdt
        Ee = g*me*c**2
        dgdt = abs(TLosses(BFF_fun(z_phys), r_fun(z_phys), v_lat(z_phys), 1d0, 1d0, Ee, 'e')) / (me*c**2)
        Bfun = dgdt / (c*betaGamma(z_phys)) * R_g
        return
    end function Bfun

    real(dp) function Rfun(g, z_phys)
        real(dp), intent(in) :: g, z_phys
        real(dp) :: Ee
        Ee = g*me*c**2				
        Rfun = 1d99
        return
    end function Rfun

    real(dp) function Qfun(g, z_phys)
        real(dp), intent(in) :: g, z_phys
        real(dp) :: Ee
        Ee = g*me*c**2
        if (z_phys < zmax_inj) then
            Qfun = dQ_Edz(Ee, z_phys, BFF_fun(z_phys), r_fun(z_phys), v_lat(z_phys), 'e') * R_g * me*c**2
        else
            Qfun = 0d0
        endif
        return
    end function Qfun

    real(dp) function safeLog10(x)
        real(dp), intent(in) :: x
        if (x > 0d0) then
            safeLog10 = log10(x)
        else
            safeLog10 = -50
        endif
        return
    end function safeLog10

    ! :::: This produces an error message ::::
    subroutine an_error(string)
        implicit none
        character(len=*), intent(IN) :: string
        write (*, "('ERROR ', A)") string
        stop 'program terminated by an_error'
    end subroutine an_error


    function assert_eq(nn, string)
        implicit none
        character(len=*), intent(in) :: string
        integer, intent(in), dimension(:) :: nn
        integer :: assert_eq
        if (all(nn(2:) == nn(1))) then
           assert_eq = nn(1)
        else
           write(*,*) 'ERROR: an assert_eq failed with this tag:', string
           stop 'program terminated by assert_eq'
        end if
    end function assert_eq

    subroutine tridag_ser(a,b,c,r,u)
        !  Description:
        !     Solver of a tridiagonal matrix. Based on the code in
        !     "Numberical Recipes".
        implicit none
        real(dp),dimension(:),intent(in) :: a,b,c,r
        real(dp),dimension(:),intent(out) :: u
        real(dp),dimension(size(b)) :: gam
        integer :: n,j
        real(dp) :: bet
        n=assert_eq((/ size(a)+1,size(b),size(c)+1,size(r),size(u) /),'tridag_ser')
        bet=b(1)
        if (bet==0.0d0) call an_error('tridag_ser: error at code stage 1')
        u(1)=r(1)/bet
        do j=2,n
           gam(j)=c(j-1)/bet
           bet=b(j)-a(j-1)*gam(j)
           if (bet==0.0d0) call an_error('tridag_ser: error at code stage 2')
           u(j)=(r(j)-a(j-1)*u(j-1))/bet
        end do
        do j=n-1,1,-1
           u(j)=u(j)-gam(j+1)*u(j+1)
        end do
    end subroutine tridag_ser

    function zeros1D(n, small) result(a)
        implicit none
        integer, intent(in) :: n
        logical, intent(in) :: small
        real(dp), dimension(n) :: a
        if (small) then
           a = 1d-200
        else
           a = 0d0
        end if
    end function zeros1D



    subroutine FP_FinDif_cool(dt, g, nin, nout, gdot, QQ, tesc)
        implicit none
        real(dp), intent(in) :: dt, tesc
        real(dp), intent(in), dimension(:) :: g, nin, gdot, QQ
        real(dp), intent(out), dimension(:) :: nout
        integer :: Ng, Ng1, k
        real(dp), dimension(size(g)) :: dx, dxp2, dxm2, BBp2, BBm2, a, b, c, r
  
        Ng = size(g)
        Ng1 = Ng - 1
  
        dxp2(:Ng1) = g(2:) - g(:Ng1)
        dxp2(Ng) = dxp2(Ng1)
        dxm2(2:) = dxp2(:Ng1)
        dxm2(1) = dxm2(2)
        dx = dsqrt(dxp2 * dxm2)
  
        BBp2(:Ng1) = -0.5d0 * (gdot(2:) + gdot(:Ng1))
        BBm2(2:) = BBp2(:Ng1)
        BBp2(Ng) = -0.5d0 * (gdot(Ng) + gdot(Ng1))
        BBm2(1) = -0.5d0 * (gdot(2) + gdot(1))
  
        r = nin + dt * QQ
        a = zeros1D(Ng,.true.)
        c = -dt * BBp2 / dx
        b = 1d0 + BBm2 * dt / dx + dt / tesc
  
        call tridag_ser(a(2:), b, c(:Ng1), r, nout)
  
    end subroutine FP_FinDif_cool

    subroutine calcNE_paramo(Evec, NEmatrix)

        real(dp), intent(in) :: Evec(nIn)
        real(dp), intent(out) :: NEmatrix(nIn, Ns)
        real(dp), allocatable, dimension(:) :: t, Ntot, gg, dt, dg, zz, Ee
        real(dp), allocatable, dimension(:,:) :: dotg, nn, Qinj, Ddiff
        real(dp) :: Ecut, tlossFF, dttFF, dzz
        integer :: jz, k, m, izz
        real(dp) :: E1, E2, z1, z2, N11, N12, N21, N22, N1, N2, slope1, slope2, slope, Ng
        real(dp) :: g_min, g_max
        
        allocate(t(0:numdt), Ntot(numdt), dt(numdt), dg(numbins), zz(0:numdt))
        allocate(nn(numbins, 0:numdt), dotg(numbins, 0:numdt), gg(numbins), &
                Qinj(numbins, 0:numdt), Ee(numbins), Ddiff(numbins, 0:numdt))
        
        t(0) = 0d0
        zz(0) = z0*1.001
        g_min = EminFF / mec2
        g_max = EmaxFF / mec2
        build_g: do k = 1, numbins
            gg(k) = g_min * (g_max / g_min)**(dble(k - 1) / dble(numbins - 1))
            Ee(k) = gg(k) * mec2
            dotg(k, 0) = -abs(TLosses(BFF_fun(zz(0)), zz(0), 1d0, 1d0, 1d0, Ee(k), 'e')) / mec2
            dotg(k, 0) = dotg(k, 0) / (betaGamma(zz(0))*c)
            if ( k > 1 ) dg(k) = gg(k) - gg(k - 1)
        end do build_g
        dg(1) = dg(2)

        build_Q: do k = 1, numbins
            Qinj(k, 0) = dQ_Edz(Ee(k), zz(0), BFF_fun(zz(0)), r_fun(zz(0)), v_lat(zz(0)), 'e') * mec2
        end do build_Q

        nn(:, 0) = 0d0
        Ddiff(:, 0) = 1d-200

        jz = 0
        z_loop: do while (zz(jz) < zmax)
            Ecut = CutE(BFF_fun(zz(jz)), r_fun(zz(jz)), v_lat(zz(jz)), 1d0, 1d0, 'e')
            tlossFF = Ecut/abs(TLosses(BFF_fun(zz(jz)), zz(jz), 1d0, 1d0, 1d0, Ecut, 'e'))
            dttFF = min(1e4*tlossFF, 1d-3*zz(jz)/(betaGamma(zz(jz))*c))
            dzz = dttFF * betaGamma(zz(jz)) * c
            
            jz = jz + 1
            zz(jz) = zz(jz-1) + dzz

            call FP_FinDif_cool(dzz, &
            &                   gg, &
            &                   nn(:, jz - 1), &
            &                   nn(:, jz), &
            &                   dotg(:, jz - 1), &
            &                   Qinj(:, jz - 1), &
            &                   1d200)

            build_dotg: do k = 1, numbins
                dotg(k, jz) = -abs(TLosses(BFF_fun(zz(jz)), zz(jz), 1d0, 1d0, 1d0, Ee(k), 'e'))/mec2!, &
                                     !nn(k, jz) / (betaGamma(zz(jz)) * c * pi*r_fun(zz(jz))**2) )) / mec2
                dotg(k, jz) = dotg(k, jz) / (betaGamma(zz(jz))*c)
            end do build_dotg

            if (zz(jz) < zmax_inj) then
                build_Qin: do k = 1, numbins
                    Qinj(k, jz) = dQ_Edz(Ee(k), zz(jz), BFF_fun(zz(jz)), r_fun(zz(jz)), v_lat(zz(jz)), 'e') * mec2
                end do build_Qin
            else
                Qinj(:, jz) = 0d0
            end if

            !write(*, *) jz, zz(jz) / R_g

        end do z_loop

        ! INTERPOLATION
        m = 1
        do i2=1, nIn
            do while (Ee(m) < Evec(i2))
                m = m + 1
            enddo
            m = max(1, m-1)
            m = min(numbins, m)
            E1 = Ee(m)
            E2 = Ee(m+1)
            jz = 0
            do izz=1, Ns
                do while (zz(jz) < z(izz))
                    jz = jz + 1
                end do
                jz = max(0, jz-1)
                jz = min(numdt-1, jz)
                z1 = zz(jz)
                z2 = zz(jz+1)
                N11 = safeLog10(nn(m, jz))
                N12 = safeLog10(nn(m, jz+1))
                N21 = safeLog10(nn(m+1, jz))
                N22 = safeLog10(nn(m+1, jz+1))
                slope1 = (N12-N11)/log10(z2/z1)
                slope2 = (N22-N21)/log10(z2/z1)
                N1 = N11 + slope1 * log10(z(izz)/z1)
                N2 = N21 + slope2 * log10(z(izz)/z1)
                slope = (N2-N1)/log10(E2/E1)
                !NEmatrix(i2, izz) = nn(m, jz)
                NEmatrix(i2, izz) = N1 + slope * log10(Evec(i2)/E1)
                NEmatrix(i2, izz) = 10**NEmatrix(i2, izz) * deltaZ(izz) / (betaGamma(z(izz)) * c * mec2)
            enddo
        enddo

    end subroutine


    subroutine calcNE_onezone(Evec, NEmatrix)
        real(dp), intent(in) :: Evec(nIn)
        real(dp), intent(out) :: NEmatrix(nIn, Ns)
        real(dp) :: tcool
        integer i2, izz

        write(*,*)"ONEZONE"
        do i2=1, nIn
            do izz=1, Ns
                tcool = Evec(i2)/abs(TLosses(BFF_fun(z(izz)), r_fun(z(izz)), v_lat(z(izz)), 1d0, 1d0, Evec(i2), 'e'))
                NEmatrix(i2, izz) = dQ_Edz(Evec(i2), z(izz), BFF_fun(z(izz)), r_fun(z(izz)), v_lat(z(izz)), 'e') * tcool * &
                 deltaZ(izz)
            enddo
        enddo

        return
    end subroutine calcNE_onezone

    subroutine CalcNE_fds(Evec, NEmatrix, Mg, Nz)
        
        ! We define a new mesh of points
        real(dp), intent(in) :: Evec(nIn)
        real(dp), intent(out) :: NEmatrix(nIn, Ns)
        integer, intent(in) :: Mg, Nz
        real(dp) :: z_vec(1000000)
        real(dp) :: paso_z
        real(dp) :: delta_z(Nz)
        real(dp) :: g_au(Mg+3), delta_g_m_au(Mg+2), delta_g(Mg+1)
        real(dp) :: g_min, g_max, paso_g
        real(dp) :: aa(Mg+1), bb(Mg+1), cc(Mg+1)
        real(dp) :: dd(Mg+1), NNee(1000000, Mg+1) !NNee(Nz, Mg+1)
        real(dp) :: Bm_minusHalf, Bm_plusHalf, Qm, Rm
        real(dp) :: E1, E2, z1, z2, N11, N12, N21, N22, N1, N2, slope1, slope2, slope
        real(dp) :: zz, dzz, dttFF, tlossFF, Ecut, z_phys
        integer :: jz, jg, m, izz, i2

        g_min = EminFF / (me*c**2)
        g_max = EmaxFF / (me*c**2)
        
        paso_g = (g_max/g_min)**(1d0/Mg)
        delta_g = 1d0
        !delta_g_m_au = 1d0
        g_au = 1d0
        g_au(1) = g_min/paso_g

        do jg=2, Mg+3
            g_au(jg) = g_au(jg-1) * paso_g
        enddo
        !do jg=1, Mg+2
        !    delta_g_m_au(jg) = g_au(jg+1) - g_au(jg)
        !enddo
        do jg=1, Mg+1
            delta_g(jg) = 0.5*(g_au(jg+2) - g_au(jg))
        enddo
            
        ! z_vec = 1d0
        ! z_vec(1) = z0 / R_g
        ! paso_z = (zmax/z0)**(1d0/(Nz-1))
        ! delta_z = 0d0
        ! do jz=2, Nz+1
        !     z_vec(jz) = z_vec(jz-1) * paso_z
        ! enddo
        ! do jz=1, Nz
        !     delta_z(jz) = z_vec(jz+1) - z_vec(jz)
        ! enddo
    

        dd = 0d0
        NNee = 0d0

        ! do jz=1, Nz
        !     aa = 0d0
        !     bb = 0d0
        !     cc = 0d0
        !     do m=1, Mg+1
        !         Bm_minusHalf = 0.5 * (Bfun(g_au(m), z_vec(jz)) + Bfun(g_au(m+1), z_vec(jz)))
        !         Bm_plusHalf = 0.5 * (Bfun(g_au(m+1), z_vec(jz)) + Bfun(g_au(m+2), z_vec(jz)))

        !         Qm = Qfun(g_au(m+1), z_vec(jz))
        !         !Rm = Rfun(g_au(m+1), z_vec(jz))
                
        !         !if (m <= Mg) then
        !             cc(m) = - delta_z(jz) * Bm_plusHalf / delta_g(m)
        !         !end if
        !         bb(m) = 1.0 + delta_z(jz) * Bm_minusHalf / delta_g(m) !+ delta_z(jz)/Rm
        !         dd(m) = Qm*delta_z(jz) + dd(m)
        !     enddo

        !     call TriDiagSys(aa, bb, cc, dd, Mg)
        !     do m=1, Mg+1
        !         NNee(jz, m) = dd(m) / (betaGamma(z_vec(jz)*R_g)*c * me*c**2)
        !     enddo
        ! enddo

        
        jz = 1
        zz = z0 / R_g
        do while (zz*R_g < zmax)
            z_phys = zz*R_g
            Ecut = CutE(BFF_fun(z_phys), r_fun(z_phys), v_lat(z_phys), 1d0, 1d0, 'e')
            tlossFF = Ecut/abs(TLosses(BFF_fun(z_phys), r_fun(z_phys), v_lat(z_phys), 1d0, 1d0, Ecut, 'e'))
            dttFF = max(min(10*tlossFF, z_phys/(betaGamma(z_phys)*c)), dttFF_min)
            dzz = dttFF * betaGamma(z_phys) * c / R_g
            aa = 0d0
            bb = 0d0
            cc = 0d0
            do m=1, Mg+1
                
                Bm_minusHalf = 0.5 * (Bfun(g_au(m), z_phys) + Bfun(g_au(m+1), z_phys))
                Bm_plusHalf = 0.5 * (Bfun(g_au(m+1), z_phys) + Bfun(g_au(m+2), z_phys))
                Qm = Qfun(g_au(m+1), z_phys)
                
                if (m <= Mg) then
                    cc(m) = - dzz * Bm_plusHalf / delta_g(m)
                end if
                bb(m) = 1.0 + dzz * Bm_minusHalf / delta_g(m) !+ delta_z(jz)/Rm
                dd(m) = Qm*dzz + dd(m)
            enddo

            call TriDiagSys(aa, bb, cc, dd, Mg)
            do m=1, Mg+1
                NNee(jz, m) = dd(m) / (betaGamma(z_phys)*c * me*c**2)
            enddo
            z_vec(jz) = z_phys
            zz = zz + dzz
            write(*, *) "jz = ", jz, "zz = ", zz
            jz = jz+1
        enddo
        
        m = 1
        do i2=1, nIn
            do while (g_au(m)*me*c**2 < Evec(i2))
                m = m + 1
            enddo
            m = max(1, m-1)
            m = min(Mg+2, m)
            E1 = g_au(m) * me*c**2
            E2 = g_au(m+1) * me*c**2
            jz = 1
            do izz=1, Ns
                do while (z_vec(jz) < z(izz))
                    jz = jz + 1
                end do
                jz = max(1, jz-1)
                jz = min(Nz-1, jz)
                z1 = z_vec(jz)
                z2 = z_vec(jz+1)
                N11 = safeLog10(NNee(jz, m))
                N12 = safeLog10(NNee(jz+1, m))
                N21 = safeLog10(NNee(jz, m+1))
                N22 = safeLog10(NNee(jz+1, m+1))
                slope1 = (N12-N11)/log10(z2/z1)
                slope2 = (N22-N21)/log10(z2/z1)
                N1 = N11 + slope1 * log10(z(izz)/z1)
                N2 = N21 + slope2 * log10(z(izz)/z1)
                slope = (N2-N1)/log10(E2/E1)
                NEmatrix(i2, izz) = N1 + slope * log10(Evec(i2)/E1)
                NEmatrix(i2, izz) = 10**NEmatrix(i2, izz) * deltaZ(izz) ! / (betaGamma(z(izz))*c * me*c**2)
            enddo
        enddo
        
        return
    end subroutine CalcNE_fds

END MODULE edistribution
