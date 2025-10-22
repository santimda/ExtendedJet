MODULE initialize
    use globals
    use constants
    use controls
    implicit none
    
    contains 

    !============================================================================
    ! initialize_allocate (START)
    !============================================================================
    subroutine initialize_allocate(sigma_value, Lj_value, z0_value, a_inj_value)
    !============================================================================
    ! allocate and initialize vectors 
    !============================================================================
        real(dp), intent(in) :: sigma_value, Lj_value, z0_value, a_inj_value
        real(dp) :: zx
        logical :: fexist

        dlogz = (zmax/z0)**(1d0/(Ns-1d0))

        allocate( z(Ns), rj_out(Ns), rj_in(Ns), gj(Ns), bj(Ns) )
        allocate( BFF(Ns), vnorm(Ns), VFF(Ns), EphMFF(Ns), nfFF(Ns), Kappa(Ns), sjoFF(NS) ) 
        allocate( deltaZ(Ns), deltatFF(Ns) )
        allocate( delta_obs(Ns), boost(Ns) )

        ! Initialize values given by free parameters:
        sigma_0 = sigma_value !Update the value of sigma_0 (default value defined in system_parameters.f90)
        Lj = Lj_value         !Update the value of Lj (default value defined in system_parameters.f90)
        !z0 = z0_value         !Update the value of z0 (default value defined in system_parameters.f90)
        a_inj = a_inj_value   !Update the value of z0 (default value defined in system_parameters.f90)

        LnT = Lj*fnT          !Available nonthermal luminosity
        LnTFF = LnT!/gj0**2   !Available nonthermal luminosity in the FF
        
        mu = gj0 * (sigma_0 + 1d0)
        Mdot_j = Lj / (c**2 * gj0 * (sigma_0 + 1d0))
        max_enthalpy = 1d0 + LnT / (lorentzFactor(zmax_inj)*Mdot_j * c**2)
        !B0 = sqrt(fB*8d0*Lj/(betaJet(z0)*c*r0**2)) / lorentzFactor(z0)
        
        !dl = (zmax-z0)/Ns    !Segment length
        !B0 = sqrt(fB*8d0*Lj/(betaJet(z0)*c*r0**2)) / lorentzFactor(z0)
        
        cte_inj = LnTFF/( (zmax_inj**(-a_inj+1.) - z0**(-a_inj+1.))/(-a_inj+1.) )

        !write(*,'(A,Es9.2)') 'B0[G] = ', BFF_fun(z0)
        !write(*,'(A,Es9.2)') 'Lj[erg/s] = ', Lj

        ! Initialize arrays
        if (.not.jet_calc) then 
            ! Initialize using the values from a previous run (if the file exists)
            inquire(file='quant_FF.dat',exist=fexist)
            if (.not. fexist) then
                print *,'You first have to run trajectory and evolution modules at least once!'
                stop
            else
                open(30,file='quant_FF.dat')  
                do i = 1, Ns
                    !Quantities at each segment
                    read(30,*) z(i), BFF(i), rj_out(i), rj_in(i), vnorm(i), Kappa(i), sjoFF(i), VFF(i), EphMFF(i), nfFF(i), deltaZ(i), deltatFF(i)
                enddo
                close(30)
            endif
        endif

    end subroutine initialize_allocate
    !============================================================================
    ! initialize_allocate (END)
    !============================================================================
    

    !============================================================================
    ! initialize_deallocate (START)
    !============================================================================
    subroutine initialize_deallocate
    !============================================================================
    ! deallocate vectors 
    !============================================================================

        deallocate( z, rj_out, rj_in, gj, bj )
        deallocate( BFF, vnorm, VFF, EphMFF, nfFF, Kappa, sjoFF )
        deallocate( deltaZ, deltatFF )
        deallocate( delta_obs, boost )

    end subroutine initialize_deallocate
    !============================================================================
    ! initialize_deallocate (END)
    !============================================================================


END MODULE initialize
