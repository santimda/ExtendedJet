MODULE synrad
    use globals
    use constants
    implicit none

contains

    !============================================================================
    ! SEDsyn (START)
    !============================================================================
    subroutine SEDsyn(jet, particle)
    !============================================================================
    ! This subroutine calculates the synchrotron SED 
    ! The input is the '' (jet) or 'CJ' (counter-jet). 
    !============================================================================
    character(LEN=*), intent(in) :: jet ! 'CJ' or 'jet'
    character(1), intent(in) :: particle
    real(dp) :: E(nIn), NE(nIn,Ns), NEcut(nIn), Ein, Eph, dEph, EphFF(Ns), dEphFF(Ns)
    real(dp) :: Lumcut(Ns), LumcutFF(Ns), L1, L2, L3, L4, L5, L6, Lum, LumFF, LBolFF
    real(dp) :: LumcutA(Ns), LumcutAFF(Ns), LumA, LumAFF, LBolAFF, z_proj, sjo
    real(dp) :: j_at_z, alpha_at_z
    logical, parameter :: bin_SEDsyn = .True.

    !N(E) of each cut
    open(30,file='N'//particle//'_cuts.dat')
    open(40,file='SED'//particle//'_syn_cuts'//jet//'.dat')
    open(45,file='n_ph_syn_cuts_FF_'//jet//'.dat')
    open(50,file='SED'//particle//'_syn_jet'//jet//'.dat')
    do i = 1, nIn
        do j = 1, Ns
            read(30,*) E(i), NE(i,j), dummy
            NE(i,j) = NE(i,j) * deltaZ(j)
        enddo
    enddo
    close(30)
    
    open(122,file='Em_Abs.dat')
    Ein = (Emax_sy/Emin_sy)**(1d0/(NSED_sy-1))  !Energy increase factor
    Eph = Emin_sy                               !Initialization
    sjo = sin(jo)

    !Initialization of the bolometric luminosities
    LBolFF = 0d0
    LBolAFF = 0d0

    do l = 1, NSED_sy

        !Initialization of the total luminosity for each photon energy
        Lum = 0d0           !Unabsorbed (de)boosted luminosity for Eph
        LumFF = 0d0         !In the FF
        LumA = 0d0          !Absorbed (de)boosted luminosity for Eph
        LumAFF = 0d0        !In the FF
        LumcutAFF = 0d0
        LumcutFF = 0d0
        LumcutA = 0d0       !Absorbed radiation for each cut

        dEph = Eph*(Ein-1d0)    !Energy interval

        !Loop for the segments
        do j = 2, Ns
            EphFF(j) = Eph/delta_obs(j)     !Photon energy in the FF
            dEphFF(j) = dEph/delta_obs(j)   !Energy interval in the FF
            NEcut = NE(:,j)                 !N(E) for the segment j

            !Unabsorbed synchrotron luminosity for each segment
            LumcutFF(j) = SyncLum(BFF(j),EphFF(j),E,NEcut,particle)

            !Emissivity and Absortion at height j
            !j_at_z=jsync(EphFF(j), BFF(j), NEcut, E)
            !alpha_at_z=alpha_sync(EphFF(j), BFF(j), NEcut, E)
            
            !write(122,*)Eph, j_at_z, alpha_at_z, BFF(j), NEcut

            Lumcut(j) = boost(j)*LumcutFF(j)      !(de)Boosted
            !Unabsorbed (de)boosted luminosity for Eph
            Lum = Lum + Lumcut(j)
            LumFF = LumFF + LumcutFF(j)*dEphFF(j)
        enddo
        
        !Unabsorbed bolometric luminosity
        LBolFF = LBolFF + LumFF
                
        !Absorbed radiation for each cut
        call SSA(BFF,rj,sjoFF,VFF,EphFF,E,NE,LumcutFF,LumcutAFF,particle)

        do j = 1, Ns
            NEcut = NE(:,j)
            !call SSA(BFF,rj,sjoFF,VFF,EphFF,E,NE,LumcutFF,LumcutAFF,particle)
            LumcutA(j) = boost(j)*LumcutAFF(j)    !(de)Boosted
            write(40,*) Eph/eV, Eph*Lumcut(j), Eph*LumcutA(j), j
            write(45,*) EphFF(j), LumcutAFF(j)/VFF(j) * rj(j)/c / EphFF(j)
            
            !Absorbed (de)boosted luminosity for Eph
            LumA = LumA + LumcutA(j)
            LumAFF = LumAFF + LumcutAFF(j)*dEphFF(j)
            
            !Jet and c-jet contribution of the segments up to given lengths
            z_proj = z(j) * sjo
            if (z_proj <= 7.3*R_g) L1 = LumA
            if (z_proj <= 2.6d1*R_g) L2 = LumA
            if (z_proj <= 6.5d1*R_g) L3 = LumA
            if (z_proj <= 1.95d2*R_g) L4 = LumA
            if (z_proj <= 6.5d2*R_g) L5 = LumA
            if (z_proj <= 3d5*R_g) L6 = LumA
            write(*,*) l, j, Ns
        enddo

        write(50,*) Eph/eV, Eph*L1, Eph*L2, Eph*L3, Eph*L4, Eph*L5, Eph*L6, Eph*Lum, Eph*LumA
        
        !Absorbed bolometric luminosity
        LBolAFF = LBolAFF + LumAFF
        
        Eph = Eph*Ein       !Increase in the photon energy

    enddo
    
    !write(*,'(3X,A,ES10.3E2,1X,A)') 'Jet unabsorbed synchrotron luminosity (FF):',LBolFF,'erg/s'
    write(*,'(3X,A,ES10.3E2,1X,A)') 'Jet absorbed synchrotron luminosity (FF):',LBolAFF,'erg/s'

    close(40)
    close(122)
    close(45)
    close(50)
    
    ! Calculate the integrated SED for different z-ranges and for the total jet
    if (bin_SEDsyn .and. jet=='') call SEDsyn_bin(jet,particle)
    
    end subroutine SEDsyn
    !============================================================================
    ! SEDsyn (END)
    !============================================================================
    

    !************************************************************************************
    !   Synchrotron luminosity per photon energy unit
    !************************************************************************************
    real(dp) function SyncLum(B,Eph,E,NE,particle)
    character(1), intent(in) :: particle
    real(dp), intent(in) :: B, Eph, E(nIn), NE(nIn)
    real(dp) :: dE, Ec, x, Psync, E_step, factor
    real(dp), parameter :: cte_Ec = sqrt(3d0/2d0)*(qe*h)/(2d0*pi*me**3*c**5)
    real(dp), parameter :: cte_sy = 1.85d0*sqrt(2d0)*qe**3/(h*me*c**2) 
    real(dp) :: f, f1, f2
    !Initialization
    SyncLum = 0d0
    factor = 1d0
    if (particle == 'p') factor = mec2/mpc2
    E_step = E(2)/E(1)

    ! The integral is solved with the trapezoid method
    Ec = cte_Ec*B*E(1)**2 * factor**3
    x = Eph/Ec
    !The approximation is only valid for 0.1 < x < 10
    if ((x < 1d-1) .or. (x > 1d1)) then
        Psync = 0d0
    else    
        Psync = cte_sy * factor * B * x**(1d0/3d0)*exp(-x)
    endif   
    f1 = Psync * NE(1)
    do i = 1, nIn-1 
        dE = E(i)*(E_step - 1.d0)     !Energy interval
        Ec = cte_Ec*B*E(i+1)**2 * factor**3
        x = Eph/Ec
        !The approximation is only valid for 0.1 < x < 10
        if ((x < 1d-1) .or. (x > 1d1)) then
            Psync = 0d0
        else    
            Psync = cte_sy * factor * B * x**(1d0/3d0)*exp(-x)
        endif   
        f2 = Psync * NE(i+1)
        f = (f1 + f2)/2d0
        SyncLum = SyncLum + f*dE
        f1 = f2
    end do

    end function SyncLum
    
    !****************************************************************************
    ! Synchrotron emitted power for a specific photon energy
    !****************************************************************************
    subroutine PSynchrotron(E_ph, E, B, resultado)
    implicit none
        real(dp), intent(in) :: E_ph, B
        real(dp), dimension(:), intent(in) :: E
        real(dp) :: K
        real(dp), dimension(:), allocatable :: x, R
        real(dp), dimension(size(E)) :: Ec
        real(dp), dimension(size(E)), intent(out) :: resultado
        integer :: j, n
    n = size(E)
    ! Define constants
    K = sqrt(3.0d0) * qe**3 / (h * me * c**2)
    ! Allocate arrays
    allocate(x(n))
    allocate(R(n))
    do j = 1, n
        ! Calculate critical energy for each element of E
        Ec(j) = (3.0d0 / (4.0d0 * pi)) * (qe * h * B) / (me**3 * c**5) * E(j)**2
        ! Calculate x and R for each element of E
        x(j) = E_ph / Ec(j)
        R(j) = (1.808 * x(j)**(1.0d0 / 3.0d0) / sqrt(1.0d0 + 3.4 * x(j)**(2.0d0 / 3.0d0))) * &
                ((1.0d0 + 2.21 * x(j)**(2.0d0 / 3.0d0) + 0.347 * x(j)**(4.0d0 / 3.0d0)) / &
                (1.0d0 + 1.353 * x(j)**(2.0d0 / 3.0d0) + 0.217 * x(j)**(4.0d0 / 3.0d0))) * exp(-x(j))
    end do

    ! Compute the emitted power for each element of E
    resultado = K * B * R
    ! Deallocate arrays
    deallocate(x)
    deallocate(R)

    end subroutine PSynchrotron

    !****************************************************************************
    ! Synchrotron emissivity for a specific photon energy
    !****************************************************************************
    function jsync(E_ph, B, N, E) result(integral)
    implicit none
    integer :: nsize  ! Number of elements in arrays
    real(dp), intent(in) :: E_ph, B
    real(dp), dimension(:), intent(in) :: E
    real(dp), dimension(:), intent(in) :: N
    real(dp) :: integral
    real(dp) , dimension(:), allocatable :: integ
    integer :: i
    nsize=size(E)
    allocate(integ(nsize))
    ! The values inside the integral
    call PSynchrotron(E_ph, E, B, integ)
    integ=integ*N
    ! The integral is solved with the trapezoid method
    integral = 0.0d0
    
    do i = 1, nsize-1
        integral = integral + 0.5d0 * (E(i+1) - E(i)) * (integ(i) + integ(i+1))
        
    end do
    integral=integral/(4d0*pi)
    deallocate(integ)

    end function jsync

    !****************************************************************************
    ! Synchrotron absorption for a specific photon energy
    !****************************************************************************
    function alpha_sync(E_ph,B, N, E) result(alpha)
    implicit none
    real(dp), dimension(:) :: E, N
    real(dp), dimension(:), ALLOCATABLE :: log_eVec, log_neVec
    real(dp), intent(in) :: B, E_ph
    real(dp) :: logne
    real(dp) :: alpha

    allocate(log_eVec(size(E)), log_neVec(size(N)))
    log_eVec = log(E)
    log_neVec = log(N)

    alpha = quad(integrand, log(E(1)*1.1d0), log(E(size(E))/1.1d0), E_ph, B, log_eVec, log_neVec)
    deallocate(log_eVec, log_neVec)

    end function alpha_sync

    !Function that calculates the interpolation
    function log_interp1d(x, xp, yp)
    real(dp) :: log_interp1d, x
    real(dp), dimension(:) :: xp, yp
    integer :: i, j, n
    real(dp) :: t
    
    n = size(xp)
    log_interp1d = 0.0  ! Initialize the output array
    if (x .lt. xp(1)) THEN 
        log_interp1d = yp(1) 
    ELSE if (x .gt. xp(n)) THEN 
        log_interp1d = yp(n) 
    ELSE 
        j = 1 
        do WHILE (x .gt. xp(j)) 
            j = j + 1 
        end do 
        t = (x - xp(j-1)) / (xp(j) - xp(j-1)) 
        log_interp1d = yp(j-1) + t * (yp(j) - yp(j-1)) 
    end if
    end function log_interp1d

    !Funtion that calculates the term inside the integral
    function integrand(loge, E_ph, bField, log_eVec,log_neVec)
    
    real(dp) :: df
    real(dp), dimension(:) :: log_eVec, log_neVec
    real(dp) :: loge
    real(dp) :: E_e, integrand
    real(dp), dimension(1) :: E_e_int 
    real(dp), intent(in) :: E_ph, bField
    real(dp), dimension(1) :: p

    E_e=exp(loge)
    E_e_int=E_e
    call compute_derivative(log_interp1d, loge, 1d-4, df, log_eVec,log_neVec)
    if (df.lt. 0. .and. df.gt.-10.) then
        call PSynchrotron(E_ph, E_e_int, bField, p)
        integrand=E_e * p(1) * E_e * (exp(log_interp1d(loge,log_eVec,log_neVec))/(E_e*E_e)) * (df-2)
    else
        integrand=0.0
    end if

    end function integrand

    !Integration function
    function quad(f, a, b, E_ph, bField , log_eVec,log_neVec) result(q)
        implicit none
        real(dp), intent(in) :: a, b, E_ph, bField
        interface
            real(8) function f(x, E_ph, bField, log_eVec, log_neVec)
            real(8) :: x
            real(8), dimension(:) :: log_eVec, log_neVec
            real(8), intent(in) :: E_ph, bField
        end function f
        end interface
        real(dp), dimension(:) :: log_eVec, log_neVec
        real(dp) :: q
        real(dp) :: h, x
        integer :: n, i

        n = 100000 ! Number of intervals, you can adjust this
        h = (b - a) / n
        q = 0.0
        do i = 0, n
            x = a + i * h
            q = q + f(x, E_ph, bField, log_eVec, log_neVec)
        end do

        q = q * h

    end function quad

    !Subroutines that calculates the derivative  d(log(N))/d(log(E)) 
    subroutine compute_derivative(f, x, h, df, log_eVec, log_neVec)
    interface
        real(8) function f(x,log_eVec, log_neVec)
        real(8) :: x
        real(8), dimension(:) :: log_eVec, log_neVec
      end function f
    end interface
    real(dp), intent(in) :: x, h
    real(dp), dimension(:) :: log_eVec, log_neVec
    real(dp),intent(out) :: df
    df = (f(x + h,log_eVec,log_neVec) - f(x - h,log_eVec,log_neVec)) / (2.0 * h)
    end subroutine compute_derivative

    !============================================================================
    ! SSA (start)
    !============================================================================
    subroutine SSA (B,r_j,sjo,V,Eph,E,NE,Lum,LumA,particle)
    !============================================================================
    ! Luminosity corrected for synchrotron self-absorption
    !============================================================================
    character, intent(in) :: particle
    real(dp), intent(in) :: B(Ns),r_j(Ns),sjo(Ns),V(Ns),Eph(Ns),E(nIn),NE(nIn,Ns),Lum(Ns)
    real(dp), intent(out) :: LumA(Ns)
    real(dp) :: Absor_coef,Emmision_coef,Source_fun,depth,tau,Ap, mass

    !Loop for the segments
    do j = 2, Ns
        !Absorption and emission coefficients, and source function, for each segment (Pacholczyck, 1970)
        mass = me
        if (particle == 'p') mass = mpc2/c**2
        Absor_coef = 2.62d0*qe**3*h**2*B(j)/(4d0*pi*Eph(j)**2*mass*V(j)) * EintSSA(B(j),Eph(j),E,NE(:,j),particle)
        Emmision_coef = Lum(j)/(4d0*pi*V(j))      !Per stereoradian
        Source_fun = Emmision_coef/max(Absor_coef,1.d-30) !!!! CHECK !!!!

        !Length that photons in the observer's direction must travel inside the segment (in the FF)
        depth = r_j(j)/sjo(j)
        !Optical depth for SSA
        tau = Absor_coef*depth / delta_obs(j)
        !tau = Absor_coef*min(depth,z(j)/gj)!/delta_obs(j)
        !Estimation of the projected area, perpendicular to the observer's direction, obtained by equating
        !the unabsorbed Lum to its value for optically thin emission: 4*pi*Ap*Source_fun*tau
        !Ap = V(j)/min(depth,z(j)/gj)
        Ap = V(j)*delta_obs(j)/depth

        !Luminosity corrected for synchrotron self-absorption
        if (tau < 0.01d0) then  !!!! CHECK !!!!
            LumA(j) = Lum(j)
        else
            LumA(j) = 4d0*pi*Ap*Source_fun*(1d0-exp(-tau))
        endif
    enddo

    end subroutine SSA
    !============================================================================
    ! SSA (END)
    !============================================================================


    !************************************************************************************
    !   Integral for the particle energy in SSA
    !************************************************************************************
    real(dp) function EintSSA(B,Eph,E,NE,particle)
    real(dp), intent(in) :: B, Eph, E(nIn), NE(nIn)
    character(1), intent(in) :: particle
    real(dp) :: dE, Ec, x, F, E_step, factor
    real(dp) :: f1, f2, sum
    real(dp), parameter :: cte_Ec = sqrt(3d0/2d0)*(qe*h)/(2d0*pi*me**3*c**5)
    !Initialization
    EIntSSA = 0d0   
    factor = 1d0
    if (particle == 'p') factor = mec2/mpc2
    E_step = E(2)/E(1)

    ! The integral is solved with the trapezoid method
    Ec = cte_Ec*B*E(1)**2 * factor**3
    x = Eph/Ec
    !The approximation is only valid for 0.1 < x < 10
    if ((x < 1d-1) .or. (x > 1d1)) then
        f1 = 0d0
    else    
        f1 = NE(1)/E(1)*exp(-x)*x**(1d0/3d0)*(x + 2d0/3d0)
    endif   
    do i = 1, nIn-1 
        dE = E(i)*(E_step - 1.d0)     !Energy interval
        Ec = cte_Ec*B*E(i+1)**2 * factor**3
        x = Eph/Ec
        !The approximation is only valid for 0.1 < x < 10
        if (x < 1d-3 .or. x > 1d1) then
            f2 = 0d0 !!!!CHECK THIS!!!!, troublesome in first cell?
        else
            f2 = NE(i+1)/E(i+1)*exp(-x)*x**(1d0/3d0)*(x + 2d0/3d0)
        endif      
        f = (f1 + f2)/2d0
        EintSSA = EintSSA + f*dE
        f1 = f2
    end do

    end function EintSSA


    !============================================================================
    ! Bin the particle energy distribution at different z ranges
    !============================================================================
    subroutine SEDsyn_bin(jet,particle)
        character(LEN=*), intent(in) :: jet,particle
        integer, parameter :: n_bins = 6 ! For the plot SEDsyn_bins
        real(dp) :: SED_bin_i(n_bins,NSED_sy), SED_total_i
        real(dp) :: SED_cut(NSED_sy,Ns), Esyn(NSED_sy)
        real(dp) :: zin, zvec(n_bins+1) 

        ! Read the synchrotron SED at each location
        open(40,file='SED'//particle//'_syn_cuts'//jet//'.dat')
        write(*,*) Ns, NSED_sy
        do i = 1, NSED_sy            
            do l = 1, Ns
                !The file contains Eph/eV, Eph*Lumcut(j), Eph*LumcutA(j), j
                read(40,*) Esyn(i), SED_cut(i,l)
            enddo
        enddo

        ! Build a logarithmically-spaced vector of z intervals
        call vector_log(z(1), z(Ns), zin, zvec ) 
        open(80,file='z_bin.dat') ! Store the z values
        do i = 1, n_bins+1
            write(80,*) zvec(i), zvec(i)/z0
        enddo
        close(80)

        open(81,file='SEDsyn'//particle//'_total.dat')
        SED_bin_i = 0.d0
        do i = 1, NSED_sy
            SED_total_i = 0.d0
            do l = 1, Ns
                ! The bin "i" is for z(i) < z < z(i+1)
                if (z(l) >= zvec(1) .and. z(l) < zvec(2)) then
                    SED_bin_i(1,i) = SED_bin_i(1,i) + SED_cut(i,l)
                elseif (z(l) >= zvec(2) .and. z(l) < zvec(3)) then
                    SED_bin_i(2,i) = SED_bin_i(2,i) + SED_cut(i,l)
                elseif (z(l) >= zvec(3) .and. z(l) < zvec(4)) then
                    SED_bin_i(3,i) = SED_bin_i(3,i) + SED_cut(i,l)
                elseif (z(l) >= zvec(4) .and. z(l) < zvec(5)) then
                    SED_bin_i(4,i) = SED_bin_i(4,i) + SED_cut(i,l)
                elseif (z(l) >= zvec(5) .and. z(l) < zvec(6)) then
                    SED_bin_i(5,i) = SED_bin_i(5,i) + SED_cut(i,l)
                elseif (z(l) >= zvec(6) .and. z(l) < zvec(7)) then
                    SED_bin_i(6,i) = SED_bin_i(6,i) + SED_cut(i,l)
                end if  
                SED_total_i = SED_total_i + SED_cut(i,l)
            end do
            write(81,*) Esyn(i), SED_total_i 
        end do
        close(81)

        open(82,file='SEDsyn'//particle//'_bin.dat')
        do j = 1,N_bins
            do i = 1, NSED_sy
                write(82,*) j, Esyn(i), SED_bin_i(j,i)
            end do
        end do
        close(82)
    end subroutine SEDsyn_bin

    !============================================================================
    ! SEDsyn (END)
    !============================================================================
    subroutine Calc_Usyn(jet)
    !============================================================================
    ! Computation of the synchrotron energy density in each segment
    ! The input is the '' (jet) or 'CJ' (counter-jet). 
    !============================================================================
        character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet 
        real(dp) :: L_tot(Ns), L_abs(Ns), Lun(NSED_sy, Ns), Lsy(NSED_sy, Ns)
        real(dp) :: Esy(NSED_sy), Esy_step, dEsy
        real(dp) :: U_sy_FF

        ! Read from a file: Esy(i,l), E*Lsy_unabs(i,l), E*Lsy_abs(i,l), l
        open(27, file='SEDe_syn_cuts_ord'//jet//'.dat')
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
            write(29,*) l, z(l)/R_g, L_tot(l), L_abs(l)/LnTFF
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
    subroutine SEDsyn_data(E_data, S_data, jet, particle)
    !============================================================================
    ! This routine calculates the intrinsic synchrotron emission in the fluid frame
    ! at the energies E_data [cgs] (vector), and returns the specific flux at each 
    ! of those energies [cgs] (also a vector). The number of points in E_data has
    ! to be smaller than m_data_max defined in globals.f90.
    ! The input is the '' (jet) or 'CJ' (counter-jet). 
    !============================================================================
    character(1), intent(in) :: particle
    real(dp), intent(in) :: E_data(:)
    real(dp), intent(out) :: S_data(m_data_max)
    character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet
    real(dp), allocatable :: Eph(:), LumA(:)
    real(dp) :: E(nIn), NE(nIn,Ns), NEcut(nIn), EphFF(Ns)
    real(dp) :: LumcutFF(Ns), LumcutAFF(Ns), LumcutA(Ns)
    real(dp) :: dil
    integer :: m_data

    !N(E) of each cut
    open(30,file='Ne_cuts.dat')
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
            NEcut = NE(:,j)                 !N(E) for the segment j
            !Unabsorbed synchrotron luminosity for each segment
            LumcutFF(j) = SyncLum(BFF(j),EphFF(j),E,NEcut,particle)
        enddo
        
        !Absorbed radiation for each cut
        call SSA(BFF,rj,sjoFF,VFF,EphFF,E,NE,LumcutFF,LumcutAFF,particle)
        
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
    
    end subroutine SEDsyn_data
    !============================================================================
    ! SEDsyn_data (END)
    !============================================================================
    
    !============================================================================
    ! SEDsyn_data (START)
    !============================================================================
    subroutine SEDsyn_z_data(E_data, z_data, S_data, jet, particle)
    !============================================================================
    ! This routine calculates the intrinsic synchrotron emission in the fluid frame
    ! at the energies E_data [cgs] (vector), and returns the specific flux at each 
    ! of those energies [cgs] (also a vector) integrated up to the height z_data. 
    ! The number of points in E_data has 
    ! to be smaller than m_data_max defined in globals.f90.
    ! The input is the '' (jet) or 'CJ' (counter-jet). 
    !============================================================================
    character(1), intent(in) :: particle
    real(dp), intent(in) :: E_data(:)
    real(dp), intent(in) :: z_data(:)
    real(dp), intent(out) :: S_data(m_data_max)
    character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet
    real(dp), allocatable :: Eph(:), LumA(:)
    real(dp) :: E(nIn), NE(nIn,Ns), NEcut(nIn), EphFF(Ns)
    real(dp) :: LumcutFF(Ns), LumcutAFF(Ns), LumcutA(Ns)
    real(dp) :: dil
    integer :: m_data

    !N(E) of each cut
    open(30,file='Ne_cuts.dat')
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
            NEcut = NE(:,j)                 !N(E) for the segment j
            !Unabsorbed synchrotron luminosity for each segment
            LumcutFF(j) = SyncLum(BFF(j),EphFF(j),E,NEcut,particle)
        enddo
        
        !Absorbed radiation for each cut
        call SSA(BFF,rj,sjoFF,VFF,EphFF,E,NE,LumcutFF,LumcutAFF,particle)
        
        LumA(l) = 0d0 ! Total luminosity at each photon energy
        j = 1
        do while ( (abs(z(j)) <= z_data(l) .and. j <= Ns ) .or. j==1) 
            LumcutA(j) = boost(j)*LumcutAFF(j)  !(de)Boosted
            LumA(l) = LumA(l) + LumcutA(j)      !Absorbed (de)boosted luminosity for Eph
            j = j + 1
        enddo
    enddo

    ! Define the output vector for the first m_data positions
    S_data = 0.d0
    dil = 4.d0*pi*dist**2
    do l = 1, m_data
      S_data(l) = LumA(l)/dil
    enddo
    
    deallocate( Eph, LumA )
    
    end subroutine SEDsyn_z_data
    !============================================================================
    ! SEDsyn_z_data (END)
    !============================================================================

!================================================================================


END MODULE synrad
