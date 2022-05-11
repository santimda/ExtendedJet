MODULE edistribution
    use constants
    use mathematics
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
    real(dp) :: Ecut, Ein, EFF, NE(Ns), NEtot, NE1, NE2, EmaxFF2, EminFF2, restEnergy, Lsy, dEdt
    real(dp) :: NEmatrix(nIn, Ns), Evec(nIn)
    real(dp) :: NEE(Ns), Nparticles, flux, Eparticles, energyFlux, vol, fluxp, eFluxp, Ljett, magnetization
    real(dp) :: sqrt_factor, dEdt_factor
    integer :: izz
    logical :: print_quantities = .False.

    if (particle == 'e') then
        open(11,file='Emax_FF.dat')
        open(20,file='NE_cuts.dat')
        open(30,file='NE_jet.dat')
        EminFF2 = EminFF
        EmaxFF2 = EmaxFF
        restEnergy = mec2
    else if (particle == 'p') then
        open(11,file='Epmax_FF.dat')
        open(20,file='Np_cuts.dat')
        open(30,file='Np_jet.dat')
        EminFF2 = EpminFF
        EmaxFF2 = EpmaxFF
        restEnergy = mpc2
    else
        print *, "FATAL ERROR: WRONG PARTICLE CLASS!!!"
        stop
    endif
        
    ! Calculate the cutoff energy for particle injection at each segment
    do i = 1, Ns
        Ecut = CutE(BFF(i),rj(i),EphMFF(i),Kappa(i),particle)   !Cutoff energy for Q(E) in the FF
        if (Ecut > EmaxFF2/3.d0) then
            print *, 'EmaxFF in globals.f90 is too low!'
            stop
        endif
        write(11,*) i, z(i), Ecut/eV
    enddo
    close(11)
    
    ! Build a logarithmically-spaced energy vector
    call vector_log(EminFF2, EmaxFF2, Ein, Evec ) 

    NEE = 0d0
    NEmatrix = 0d0
    call calcNE_paramo(Evec, NEmatrix, particle)

    if (print_quantities) then
        flux = 0d0
        energyFlux = 0d0
        sqrt_factor = sqrt(Ein) - 1d0/sqrt(Ein)
        open(50, file='Quantities.dat')
        do izz = 1, Ns
            Nparticles = 0d0
            Eparticles = 0d0
            Lsy = 0d0
            dEdt_factor = 4d0/3d0 * sigma_T * c * (BFF_fun(z(izz))**2/(8d0*pi))
            do i2 = 1, nIn
                Nparticles = Nparticles + NEmatrix(i2, izz) * Evec(i2) * sqrt_factor
                Eparticles = Eparticles + NEmatrix(i2, izz) * Evec(i2)**2 * sqrt_factor
                dEdt = dEdt_factor * (Evec(i2)/mec2)**2
                Lsy = Lsy + NEmatrix(i2, izz) * dEdt * Evec(i2) * sqrt_factor
            enddo
            vol = pi * r_fun(z(izz))**2 * deltaZ(izz)
            Ljett = pi*r_fun(z(izz))**2*lorentzFactor(z(izz))**2*betaJet(z(izz))*c * (Eparticles/(vol*lorentzFactor(z(izz))) &
             + BFF_fun(z(izz))**2/(8d0*pi) + Lsy/(4d0*r_fun(z(izz))**2*c))
            magnetization = BFF_fun(z(izz))**2/(8d0*pi*Eparticles/vol)
            write(50, *) z(izz)/R_g, Eparticles/(vol*lorentzFactor(z(izz))), BFF_fun(z(izz))**2/(8d0*pi), Ljett/Lj, magnetization, &
                            (Ljett + Lsy)/Lj
            fluxp = flux
            eFluxp = energyFlux
            flux = Nparticles / vol * lorentzFactor(z(izz)) * betaJet(z(izz)) * c * pi * r_fun(z(izz))**2
            energyFlux = Eparticles / vol * betaJet(z(izz)) * c * pi * r_fun(z(izz))**2
            print *, "z = ", z(izz)/R_g, "Flux = ", flux, "EnergyFlux = ", energyFlux + Lsy, "RatioPart = ", & 
                               fluxp/flux, "RatioE = ", eFluxp/energyFlux
        enddo
        close(50)
    endif

    EFF = EminFF2       !Initialization
    do i2 = 1, nIn
        NEtot = 0.d0    !Initialization
        do izz = 1, Ns
            write(20,*) EFF, NEmatrix(i2, izz), izz     !E and N(E) for each jet cut
            NEtot = NEtot + NEmatrix(i2, izz)           !Sum for all the jet segments
            NEE(izz) = NEE(izz) + NEmatrix(i2, izz) * EFF * (Ein-1d0)

            !Contribution of the segments up to given lengths
            if (z(izz) <= 0.01d0*zmax) NE1 = NEtot
            if (z(izz) <= 0.1d0*zmax) NE2 = NEtot
        enddo
        write(30,*) EFF/eV, EFF**2*NE1, EFF**2*NE2, EFF**2*NEtot
        EFF = EFF*Ein   !Increase the energy of the loop
    enddo

    close(20)
    close(30)
    
    end subroutine Edist
        

    !************************************************************************************
    !   Determination of the cutoff energy of the particle injection
    !************************************************************************************
    real(dp) function CutE (BFF1,rj1,Eph1,Kappa1, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: BFF1,rj1,Eph1,Kappa1
    real(dp) :: Emax1,Emax2,Eacc,lE1,lE2,lEmid,check
    real(dp), parameter :: error_tolerance = 1.d-3 ! log(1 + DE/E) ~ |DE|/E <= error_tolerance
    
    !Initialization
    lE1 = log10(1d9*eV)  ! log10(Emin [CGS]) ~ 1 GeV
    lE2 = log10(1d14*eV) ! log10(Emax [CGS]) ~ 100 TeV
    if (particle == 'p') lE2 = log10(1d20*eV) !Maximum particle energy    

    Eacc = facc*qe*BFF1*c ! energy gain: dE/dt_acc [erg/s]
    do while (abs(lE2-lE1) > error_tolerance)
        lEmid = (lE1+lE2)/2d0
        check =  Eacc - abs(TLosses(BFF1,rj1,Eph1,Kappa1,10d0**lEmid, particle))
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
    !   Injection function (at a specific injection cell)
    !************************************************************************************
    real(dp) function Q_E (K, Ecut, E, EminFF2, EmaxFF2)
    real(dp), intent(in) :: K, Ecut, E, EminFF2, EmaxFF2
    !Q(E) = 0 if outside the injection energy limits
    if (E <= EmaxFF2 .and. E >= EminFF2) then
        Q_E = K*(E/EminFF2)**(-p) * exp(-E/Ecut)
    else
        Q_E = 0d0
    endif
    end function Q_E

    !************************************************************************************
    !   Injection normalization
    !************************************************************************************  
    real(dp) function Norm_q0p(Ecut, zz, particle)
    real(dp), intent(in) :: Ecut, zz
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
    Einc = (EmaxFF2/EminFF2)**(1d0/(N-1d0))  !Increase factor of the energy
    !Initialization
    E = EminFF2
    QInt = 0d0
    ! The integral is solved with the trapezoid method
    f1 = E * Q_E(1d0, Ecut, E, EminFF2, EmaxFF2)
    do i = 1, N-1 
        dE = E*(Einc-1d0)                    !Energy interval
        E2 = E*Einc
        f2 = E2 * Q_E(1d0, Ecut, E2, EminFF2, EmaxFF2)
        f = (f1 + f2)/2d0
        QInt = QInt + f*dE
        f1 = f2
        E = E2
    end do
    Norm_q0p = dLe_dz(zz)/(QInt*lorentzFactor(zz)*pi*r_fun(zz)**2)   !Normalization constant
    end function Norm_q0p

    !************************************************************************************
    !   Injection function (for multiple injection cells)
    !************************************************************************************
    real(dp) function qp_E(E, zz, BB, rr, particle)
    real(dp), intent(in) :: E, zz, BB, rr
    real(dp) :: Ecut, Kinj, EminFF2, EmaxFF2, intInj, q0p
    character(1), intent(in) :: particle
    Ecut = CutE(BB, rr, 1d0, 1d0, particle)   !Cutoff energy for Q(E) in the FF

    !Kinj = NormKz(Ecut, particle) / (zz * log(zmax_inj/z0))         !Normalization constant for Q(E)
    !intInj = z0 / (a_inj-1d0) * (1d0 - (z0/zmax_inj)**(a_inj-1d0))
    !Kinj = NormKz(Ecut, particle) * (z0/zz)**a_inj / intInj
    q0p = Norm_q0p(Ecut, zz, particle)
    !intInj = z0 / (a_inj - 1d0 - geomIndex) * (1d0 - (z0/zmax_inj)**(a_inj - 1d0 - geomIndex))
    !Kinj = NormKz(Ecut, particle) * (r0/r_fun(zz))**a_inj / intInj
    if (particle == 'e') then
        EminFF2 = EminFF
        EmaxFF2 = EmaxFF
    else if (particle == 'p') then
        EminFF2 = EpminFF
        EmaxFF2 = EpmaxFF
    endif
    if (E >= EminFF2 .and. E <= EmaxFF2) then
        qp_E = q0p*(E/EminFF2)**(-p) * exp(-E/Ecut)
    else
        qp_E = 0d0
    endif
    end function qp_E

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

    end subroutine Calc_Ue
    !============================================================================
    ! Calc_Ue (END)
    !============================================================================

    subroutine FP_FinDif_cool(dt, g, nin, nout, gdot, QQ, tesc)
        implicit none
        real(dp), intent(in) :: dt, tesc
        real(dp), intent(in), dimension(:) :: g, nin, gdot, QQ
        real(dp), intent(out), dimension(:) :: nout
        integer :: Ng, Ng1
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

    !============================================================================
    ! calcNE_paramo (START)
    !============================================================================
    subroutine calcNE_paramo(Evec, NEmatrix, particle)

        real(dp), intent(in) :: Evec(nIn)
        character(1), intent(in) :: particle
        real(dp), intent(out) :: NEmatrix(nIn, Ns)
        real(dp), allocatable, dimension(:) :: t, Ntot, gg, dt, dg, zz, Ee
        real(dp), allocatable, dimension(:,:) :: dotg, nn, Qinj, Ddiff
        real(dp) :: Ecut, tlossFF, dttFF, dzz, EminFF2, EmaxFF2, mc2
        integer :: jz, k, m, izz
        real(dp) :: E1, E2, z1, z2, N11, N12, N21, N22, N1, N2, slope1, slope2, slope
        real(dp) :: g_min, g_max
        
        allocate(t(0:numdt), Ntot(numdt), dt(numdt), dg(numbins), zz(0:numdt))
        allocate(nn(numbins, 0:numdt), dotg(numbins, 0:numdt), gg(numbins), &
                Qinj(numbins, 0:numdt), Ee(numbins), Ddiff(numbins, 0:numdt))
        
        if (particle == 'e') then
            EminFF2 = EminFF
            EmaxFF2 = EmaxFF
            mc2 = mec2
        else if (particle == 'p') then
            EminFF2 = EpminFF
            EmaxFF2 = EpmaxFF
            mc2 = mpc2
        endif

        t(0) = 0d0
        zz(0) = z0*1.001
        g_min = EminFF / mc2
        g_max = EmaxFF / mc2
        build_g: do k = 1, numbins
            gg(k) = g_min * (g_max / g_min)**(dble(k - 1) / dble(numbins - 1))
            Ee(k) = gg(k) * mc2
            dotg(k, 0) = -abs(TLosses(BFF_fun(zz(0)), zz(0), 1d0, 1d0, Ee(k), particle)) / mc2
            dotg(k, 0) = dotg(k, 0) / (betaGamma(zz(0)) *c)
            if ( k > 1 ) dg(k) = gg(k) - gg(k - 1)
        end do build_g
        dg(1) = dg(2)

        build_Q: do k = 1, numbins
            !Qinj(k, 0) = dQ_Edz(Ee(k), zz(0), BFF_fun(zz(0)), r_fun(zz(0)), v_lat(zz(0)), 'e') * mec2
            Qinj(k, 0) = qp_E(Ee(k), zz(0), BFF_fun(zz(0)), r_fun(zz(0)), particle) * mc2 * pi * r_fun(zz(0))**2
        end do build_Q

        nn(:, 0) = 0d0
        Ddiff(:, 0) = 1d-200

        jz = 0
        z_loop: do while (zz(jz) < zmax)
            Ecut = CutE(BFF_fun(zz(jz)), r_fun(zz(jz)), 1d0, 1d0, particle)
            tlossFF = Ecut/abs(TLosses(BFF_fun(zz(jz)), zz(jz), 1d0, 1d0, Ecut, particle))
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
                dotg(k, jz) = -abs(TLosses(BFF_fun(zz(jz)), zz(jz), 1d0, 1d0, Ee(k), particle))/mc2!, &
                                     !nn(k, jz) / (betaGamma(zz(jz)) * c * pi*r_fun(zz(jz))**2) )) / mec2
                dotg(k, jz) = dotg(k, jz) / (betaGamma(zz(jz))*c)
            end do build_dotg

            if (zz(jz) < zmax_inj) then
                build_Qin: do k = 1, numbins
                    !Qinj(k, jz) = dQ_Edz(Ee(k), zz(jz), BFF_fun(zz(jz)), r_fun(zz(jz)), v_lat(zz(jz)), 'e') * mec2
                    Qinj(k, jz) = qp_E(Ee(k), zz(jz), BFF_fun(zz(jz)), r_fun(zz(jz)), particle) * mc2 * pi * r_fun(zz(jz))**2
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
                NEmatrix(i2, izz) = 10**NEmatrix(i2, izz) * deltaZ(izz) / (betaGamma(z(izz)) * c * mc2)
            enddo
        enddo

    end subroutine
    !============================================================================
    ! calcNE_paramo (END)
    !============================================================================
    

END MODULE edistribution