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
    real(dp) :: Ecut, Ein, NEtot, EmaxFF2, EminFF2, restEnergy, Lsy, dEdt
    real(dp) :: NEmatrix(nIn, Ns), Evec(nIn)
    real(dp) :: Nparticles, flux, Eparticles, energyFlux, vol, fluxp, eFluxp, Ljett, magnetization
    real(dp) :: sqrt_factor, dEdt_factor
    integer :: izz
    logical :: print_quantities = .False.

    if (particle == 'e') then
        EminFF2 = EminFF
        EmaxFF2 = EmaxFF
        restEnergy = mec2
    else if (particle == 'p') then
        EminFF2 = EpminFF
        EmaxFF2 = EpmaxFF
        restEnergy = mpc2
    else
        print *, "FATAL ERROR: WRONG PARTICLE CLASS!!!"
        stop
    endif
        
    open(11,file='E'//particle//'_max_FF.dat')
    open(20,file='N'//particle//'_cuts.dat')

    ! Calculate the cutoff energy for particle injection at each segment
    do i = 1, Ns
        Ecut = CutE(BFF(i),z(i),EphMFF(i),Kappa(i),particle)   !Cutoff energy for Q(E) in the FF
        if (Ecut > EmaxFF2/2.d0) then
            Ecut = EmaxFF2/2d0
            !print *, 'EmaxFF in system_parameters.f90 is too low!'
            !stop
        endif
        write(11,*) i, z(i), Ecut/eV
    enddo
    close(11)
    
    ! Build a logarithmically-spaced energy vector
    call vector_log(EminFF2, EmaxFF2, Ein, Evec )

    NEmatrix = 0d0
    !call calcNE_onezone(Evec, NEmatrix, particle)
    call calcNE_paramo(Evec, NEmatrix, particle)

    NEmatrix(:,Ns) = 0d0

    if (print_quantities) then
        flux = 0d0
        energyFlux = 0d0
        sqrt_factor = sqrt(Ein) - 1d0/sqrt(Ein)
        open(50, file='Quantities.dat')
        do izz = 1, Ns
            Nparticles = 0d0
            Eparticles = 0d0
            Lsy = 0d0
            dEdt_factor = 4d0/3d0 * sigma_T * c * (BFF(izz)**2/(8d0*pi))
            do i2 = 1, nIn
                Nparticles = Nparticles + NEmatrix(i2, izz) * Evec(i2) * sqrt_factor
                Eparticles = Eparticles + NEmatrix(i2, izz) * Evec(i2)**2 * sqrt_factor
                dEdt = dEdt_factor * (Evec(i2)/mec2)**2
                Lsy = Lsy + NEmatrix(i2, izz) * dEdt * Evec(i2) * sqrt_factor
            enddo
            vol = pi*(rj_out(izz)**2-rj_in(izz)**2) * deltaZ(izz)
            Ljett = pi*(rj_out(izz)**2-rj_in(izz)**2)*lorentzFactor(z(izz))**2*betaJet(z(izz))*c * (Eparticles/(vol*lorentzFactor(z(izz))) &
             + BFF(izz)**2/(8d0*pi) + Lsy/(4d0*(rj_out(izz)**2-rj_in(izz)**2)*c))
            magnetization = BFF(izz)**2/(8d0*pi*Eparticles/vol)
            write(50, *) z(izz)/R_g, Eparticles/(vol*lorentzFactor(z(izz))), BFF(izz)**2/(8d0*pi), Ljett/Lj, magnetization, &
                            (Ljett + Lsy)/Lj
            fluxp = flux
            eFluxp = energyFlux
            flux = Nparticles / vol * lorentzFactor(z(izz)) * betaJet(z(izz)) * c * pi*(rj_out(izz)**2-rj_in(izz)**2)
            energyFlux = Eparticles / vol * betaJet(z(izz)) * c * pi*(rj_out(izz)**2-rj_in(izz)**2)
            print *, "z = ", z(izz)/R_g, "Flux = ", flux, "EnergyFlux = ", energyFlux + Lsy, "RatioPart = ", & 
                               fluxp/flux, "RatioE = ", eFluxp/energyFlux
        enddo
        close(50)
    endif

    ! Save the particle energy distribution (JUST MAKE THIS A SHARED OBJECT?)
    do i2 = 1, nIn
        do izz = 1, Ns
            write(20,*) Evec(i2), NEmatrix(i2, izz), izz     !E and N(E) for each jet cut
        enddo
    enddo
    close(20)

    ! Calculate the integrated distribution along different z-ranges and for the total jet
    call N_bin(Evec,NEmatrix,particle)

    end subroutine Edist
        

    !************************************************************************************
    !   Determination of the cutoff energy of the particle injection (t_acc=t_loss or Hillas)
    !************************************************************************************
    real(dp) function CutE (BFF1, z1, Eph1, Kappa1, particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: BFF1,z1,Eph1,Kappa1
    real(dp) :: Emax1,Emax2,Eacc,lE1,lE2,lEmid,check
    real(dp), parameter :: error_tolerance = 1.d-3 ! log(1 + DE/E) ~ |DE|/E <= error_tolerance
    
    !Initialization
    lE1 = log10(1d9*eV)  ! log10(Emin [CGS]) ~ 1 GeV
    lE2 = log10(2d14*eV) ! log10(Emax [CGS]) ~ 200 TeV
    if (particle == 'p') lE2 = log10(1d20*eV) !Maximum particle energy    

    Eacc = facc*qe*BFF1*c ! energy gain: dE/dt_acc [erg/s]
    do while (abs(lE2-lE1) > error_tolerance)
        lEmid = (lE1+lE2)/2d0
        check =  Eacc - abs(TLosses(BFF1,z1,Eph1,Kappa1,10d0**lEmid, particle))
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
    Emax2 = r_fun_out(z1)*sqrt(3d0*qe*BFF1*Eacc/(2d0*c))
    !Ecut obtained by the minimum of the values above
    CutE = min(Emax1,Emax2)
    end function CutE

    !************************************************************************************
    !   Determination of the break energy of the particle injection  (CONDITION???)
    !************************************************************************************
    real(dp) function BreakE(BFF1,zzz,deltaZZ,particle)
        real(dp), intent(in) :: BFF1, zzz, deltaZZ
        character(1), intent(in) :: particle
        real(dp) :: mc2

        if (particle == 'e') then
            mc2 = mec2
        else
            mc2 = mpc2
        endif
        BreakE = 6*pi*mc2**2 * betaGamma(zzz) / (sigma_T*BFF1**2*deltaZZ*(mec2/mc2)**2)
    end function BreakE

    !************************************************************************************
    !   Injection function (at a specific injection cell)
    !************************************************************************************
    real(dp) function Q_E (K, Ecut, E, EminFF2, EmaxFF2)
        implicit none
        real(dp), intent(in) :: K, Ecut, E, EminFF2, EmaxFF2

        !Q(E) = 0 if outside the injection energy limits
        if (E > EmaxFF2 .or. E < EminFF2) then
            Q_E = 0d0
            return
        end if

        ! If E/Ecut is very small, skip exp()
        Q_E = K * (E / EminFF2)**(-p)
        if (E/Ecut > 5d-4) then
            Q_E = Q_E * exp(-E/Ecut)
        end if
    end function Q_E
    !************************************************************************************
    !   Injection normalization (with parallelization)
    !************************************************************************************  
    real(dp) function Norm_q0p(Ecut, zz, particle)
        real(dp), intent(in) :: Ecut, zz
        character(1), intent(in) :: particle
        real(dp) :: EminFF2, EmaxFF2, Einc, QInt, dE
        real(dp), allocatable :: Egrid(:), f(:)
        integer, parameter :: N = 120

        ! choose ranges depending on particle
        if (particle == 'e') then
            EminFF2 = EminFF
            EmaxFF2 = EmaxFF
        else
            EminFF2 = EpminFF
            EmaxFF2 = EpmaxFF
        endif

        ! allocate
        allocate(Egrid(N))
        allocate(f(N))

        ! geometric grid
        Einc = (EmaxFF2/EminFF2)**(1d0/(N-1d0))  !Increase factor of the energy
        Egrid(1) = EminFF2
        do i = 2, N
            Egrid(i) = Egrid(i-1) * Einc
        end do

        ! compute integrand values in parallel (thread-safe if Q_E is)
        !$omp parallel do default(shared) private(i) schedule(static)
        do i = 1, N
            f(i) = Egrid(i) * Q_E(1d0, Ecut, Egrid(i), EminFF2, EmaxFF2)
        end do
        !$omp end parallel do

        ! trapezoidal integration
        QInt = 0d0
        !$omp parallel do default(shared) private(i,dE) reduction(+:QInt) schedule(static)
        do i = 1, N-1
            dE = Egrid(i) * (Einc - 1d0)
            QInt = QInt + 0.5d0*(f(i) + f(i+1)) * dE
        end do
        !$omp end parallel do

        Norm_q0p = dLe_dz(zz) / ( QInt * lorentzFactor(zz) * pi * ( r_fun_out(zz)**2 - r_fun_in(zz)**2 ) )

        deallocate(Egrid)
        deallocate(f)
    end function Norm_q0p


    !************************************************************************************
    !   Injection function (for multiple injection cells)
    !************************************************************************************
    real(dp) function qp_E(E, zz, BB, particle)
    real(dp), intent(in) :: E, zz, BB
    real(dp) :: Ecut, EminFF2, EmaxFF2, q0p, fact
    character(1), intent(in) :: particle
    Ecut = CutE(BB, zz, 1d0, 1d0, particle)   !Cutoff energy for Q(E) in the FF
    q0p = Norm_q0p(Ecut, zz, particle)        !Normalization constant
    if (particle == 'e') then
        EminFF2 = EminFF
        EmaxFF2 = EmaxFF
    else if (particle == 'p') then
        EminFF2 = EpminFF
        EmaxFF2 = EpmaxFF
    endif
    EminFF2 = 1.1*mec2 !MAX(EminFF2, (p-2d0)/(p-1d0)*mpc2*(1d0+sigma_B(zz)))
    !fact = (mpc2-mec2) * (1.0 - exp(-zz/(z0*1000)))
    !EminFF2 = MAX(EminFF2, (p-2d0)/(p-1d0)*(1d0+sigma_B(zz))* (mec2 + fact))
    if (E >= EminFF2 .and. E <= EmaxFF2) then
        qp_E = q0p*(E/EminFF2)**(-p) * exp(-E/Ecut)
    else
        qp_E = 0d0
    endif
    end function qp_E

    !============================================================================
    ! Bin the particle energy distribution at different z ranges
    !============================================================================
    subroutine N_bin(E,N_tot,particle)
        character(LEN=*), intent(in) :: particle
        real(dp), intent(in) :: E(nIn),N_tot(nIn,Ns)
        integer, parameter :: N_bins = 6 ! For the plot dist_particle_bins
        real(dp) :: N_bin_i(N_bins,nIn), N_total_i
        real(dp) :: zin, zvec(N_bins+1) 

        N_bin_i = 0.d0

        ! Build a logarithmically-spaced vector of z intervals
        call vector_log(z(1), z(Ns), zin, zvec ) 
        open(80,file='z_bin.dat') ! Store the z values
        do i = 1, N_bins+1
            write(80,*) zvec(i), zvec(i)/z0
        enddo
        close(80)

        open(81,file='N'//particle//'_total.dat')
        do i = 1, nIn
            N_total_i = 0.d0
            do l = 1, Ns
                ! The bin "i" is for z(i) < z < z(i+1)
                if (z(l) >= zvec(1) .and. z(l) < zvec(2)) then
                    N_bin_i(1,i) = N_bin_i(1,i) + N_tot(i,l) * VFF(l)
                elseif (z(l) >= zvec(2) .and. z(l) < zvec(3)) then
                    N_bin_i(2,i) = N_bin_i(2,i) + N_tot(i,l) * VFF(l)
                elseif (z(l) >= zvec(3) .and. z(l) < zvec(4)) then
                    N_bin_i(3,i) = N_bin_i(3,i) + N_tot(i,l) * VFF(l)
                elseif (z(l) >= zvec(4) .and. z(l) < zvec(5)) then
                    N_bin_i(4,i) = N_bin_i(4,i) + N_tot(i,l) * VFF(l)
                elseif (z(l) >= zvec(5) .and. z(l) < zvec(6)) then
                    N_bin_i(5,i) = N_bin_i(5,i) + N_tot(i,l) * VFF(l)
                elseif (z(l) >= zvec(6) .and. z(l) < zvec(7)) then
                    N_bin_i(6,i) = N_bin_i(6,i) + N_tot(i,l) * VFF(l)
                end if  
                N_total_i = N_total_i + N_tot(i,l) * VFF(l)
            end do
            write(81,*) E(i)/eV, N_total_i 
        end do
        close(81)

        open(82,file='N'//particle//'_bin.dat')
        do j = 1,N_bins
            do i = 1, nIn
                write(82,*) j, E(i)/eV, N_bin_i(j,i)
            end do
        end do
        close(82)
    end subroutine N_bin

    !============================================================================
    ! Computation of the electron energy density in each segment (FF)
    !============================================================================
    subroutine Calc_Ue()
        real(dp) :: Ue, Ne(nIn, Ns), Ee(nIn), Ee_step, dEe, E_tot

        ! Read from a file: Ee(i), Ne(i,l)
        open(70,file='Ne_cuts.dat')
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

    subroutine FP_FinDif_cool(dt, g, n_in, n_out, gdot, QQ, tesc)
        implicit none
        real(dp), intent(in) :: dt, tesc
        real(dp), intent(in), dimension(:) :: g, n_in, gdot, QQ
        real(dp), intent(out), dimension(:) :: n_out
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
  
        r = n_in + dt * QQ
        a = zeros1D(Ng,.true.)
        c = -dt * BBp2 / dx
        b = 1d0 + BBm2 * dt / dx + dt / tesc
        call tridag_ser(a(2:), b, c(:Ng1), r, n_out)
  
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
        real(dp) :: Ecut, tlossFF, dttFF, dzz, EminFF2, EmaxFF2, mc2, Bz0, Bjz, const_z0, const_zjz
        integer :: jz, k, m, izz
        real(dp) :: E1, E2, z1, z2, N11, N12, N21, N22, N1, N2, slope1, slope2, slope
        real(dp) :: g_min, g_max, Ebr, dz_aux
        
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
        g_min = EminFF2 / mc2
        g_max = EmaxFF2 / mc2
        !$OMP PARALLEL DO PRIVATE(k)
        build_g: do k = 1, numbins
            gg(k) = g_min * (g_max / g_min)**(dble(k - 1) / dble(numbins - 1))
            Ee(k) = gg(k) * mc2
            Bz0 = BFF_fun(zz(0))
            dotg(k, 0) = -abs(TLosses(Bz0, zz(0), 1d0, 1d0, Ee(k), particle)) / mc2             ! dotg'
            dotg(k, 0) = dotg(k, 0) / (betaGamma(zz(0)) *c)                                                ! dotg' / (betaGamma c)
            if ( k > 1 ) dg(k) = gg(k) - gg(k - 1)
        end do build_g
        !$OMP END PARALLEL DO
        dg(1) = dg(2)

        const_z0 = mc2 * pi *(r_fun_out(zz(0))**2-r_fun_in(zz(0))**2)
        !$OMP PARALLEL DO PRIVATE(k)
        build_Q: do k = 1, numbins
            Qinj(k, 0) = qp_E(Ee(k), zz(0), Bz0, particle) * const_z0        ! dQ'(g', z)/dz
        end do build_Q
        !$OMP END PARALLEL DO

        nn(:, 0) = 0d0
        Ddiff(:, 0) = 1d-200

        jz = 0
        z_loop: do while (zz(jz) < zmax*1.01)
            dz_aux = zz(jz)*(sqrt(dlogz)-1d0/sqrt(dlogz))
            Bjz = BFF_fun(zz(jz))
            Ebr = BreakE(Bjz, zz(jz), dz_aux, particle)
            Ecut = CutE(Bjz, zz(jz), 1d0, 1d0, particle)
            Ebr = min(Ebr, Ecut)
            tlossFF = Ebr/abs(TLosses(Bjz, zz(jz), 1d0, 1d0, Ebr, particle))
            dttFF = min(0.01*tlossFF, dz_aux/(betaGamma(zz(jz))*c))
            !if ((1d0*tlossFF, 1d-2*zz(jz)/(betaGamma(zz(jz))*c))
            dzz = dttFF * betaGamma(zz(jz)) * c
            !write(*,*) Ebr, Ecut, Bjz, sigma_B(zz(jz))
            jz = jz + 1
            zz(jz) = zz(jz-1) + dzz

            call FP_FinDif_cool(dzz, &
            &                   gg, &
            &                   nn(:, jz - 1), &
            &                   nn(:, jz), &
            &                   dotg(:, jz - 1), &
            &                   Qinj(:, jz - 1), &
            &                   1d200)

            !write(*,*) nn(:,jz)
            !$OMP PARALLEL DO PRIVATE(k)
            build_dotg: do k = 1, numbins
                dotg(k, jz) = -abs(TLosses(Bjz, zz(jz), 1d0, 1d0, Ee(k), particle))/mc2!, &
                                     !nn(k, jz) / (betaGamma(zz(jz)) * c * pi*r_fun(zz(jz))**2) )) / mec2
                dotg(k, jz) = dotg(k, jz) / (betaGamma(zz(jz))*c)
            end do build_dotg
            !$OMP END PARALLEL DO

            const_zjz = mc2 * pi * (r_fun_out(zz(jz))**2-r_fun_in(zz(jz))**2)* exp(-(zz(jz)/zmax_inj))
            !if (zz(jz) < zmax_inj) then
            !$OMP PARALLEL DO PRIVATE(k)
            build_Qin: do k = 1, numbins
              Qinj(k, jz) = qp_E(Ee(k), zz(jz), Bjz, particle) * const_zjz
            end do build_Qin
            !$OMP END PARALLEL DO
            !else
            !    Qinj(:, jz) = 0d0
            !end if

        end do z_loop

        ! INTERPOLATION
        !$OMP PARALLEL DO PRIVATE(i2, m, E1, E2, jz, izz, z1, z2, N11, N12, N21, N22, N1, N2, slope1, slope2, slope)
        do i2=1, nIn
            m = 1
            do while (Ee(m) < Evec(i2))
                m = m + 1
            enddo
            m = max(1, m-1)
            m = min(numbins, m)
            E1 = Ee(m)
            E2 = Ee(m+1)
            jz = 0
            do izz=2, Ns
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
                slope1 = (N12-N11)/safeLog10(z2/z1)
                slope2 = (N22-N21)/safeLog10(z2/z1)
                N1 = N11 + slope1 * safeLog10(z(izz)/z1)
                N2 = N21 + slope2 * safeLog10(z(izz)/z1)
                slope = (N2-N1)/safeLog10(E2/E1)
                !NEmatrix(i2, izz) = nn(m, jz)
                NEmatrix(i2, izz) = N1 + slope * safeLog10(Evec(i2)/E1)
                NEmatrix(i2, izz) = 10**NEmatrix(i2, izz) / ((pi * (rj_out(izz)**2-rj_in(izz)**2)) * betaGamma(z(izz)) * c * mc2)
            enddo
        enddo
        !$OMP END PARALLEL DO

    end subroutine
    !============================================================================
    ! calcNE_paramo (END)
    !============================================================================
    

    !============================================================================
    ! ONEZONE
    !============================================================================
    subroutine calcNE_onezone(Evec, NEmatrix, particle)

        real(dp), intent(in) :: Evec(nIn)
        character(1), intent(in) :: particle
        real(dp), intent(out) :: NEmatrix(nIn, Ns)
        real(dp) :: t_cool, t_adv, Q
        integer :: i, j

        do i=1,nIn
            do j=1,Ns
                Q = qp_E(Evec(i), z(j), BFF_fun(z(j)), particle) * exp(-(z(j)/zmax_inj))
                t_cool = Evec(i)/abs(TLosses(BFF_fun(z(j)), z(j), 1d0, 1d0, Evec(i), particle))
                t_adv = z(j)/(betaGamma(z(j))*c)
                NEmatrix(i, j) = Q * min(t_cool, t_adv)
            end do
        end do
    end subroutine

END MODULE edistribution
