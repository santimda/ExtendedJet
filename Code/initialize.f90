MODULE initialize
    use globals
    use constants
    use controls
    implicit none
    
    contains 

    !============================================================================
    ! initialize_allocate (START)
    !============================================================================
    subroutine initialize_allocate(sigma_value, zmax_value, Lj_value, a_inj_value) !,z0_value
    !============================================================================
    ! allocate and initialize vectors 
    !============================================================================
        real(dp), intent(in) :: sigma_value, zmax_value, Lj_value, a_inj_value !,z0_value
        integer :: i
        logical :: fexist

        ! Initialize values given by free parameters:
        sigma_0 = sigma_value !Update the value of sigma_0 (default value defined in system_parameters.f90)
        Lj = Lj_value         !Update the value of Lj (default value defined in system_parameters.f90)
        a_inj = a_inj_value   !Update the value of a_inj (default value defined in system_parameters.f90)
        zmax_rg = zmax_value/R_g   !zmax is derived from zmax_rg (R_g is set by read_namelist/derive_parameters)

        ! Recompute every quantity tied to the parameters above (rout_0, mu,
        ! Mdot_j, LnT, LnTFF, max_enthalpy, cte_inj, dlogz, ...)
        call derive_parameters()

        allocate( z(Ns), rj_out(Ns), rj_in(Ns), gj(Ns), bj(Ns) )
        allocate( BFF(Ns), vnorm(Ns), VFF(Ns), EphMFF(Ns), nfFF(Ns), Kappa(Ns), sjoFF(NS) )
        allocate( deltaZ(Ns), deltatFF(Ns) )
        allocate( delta_obs(Ns), boost(Ns) )

        ! Initialize arrays
        if (.not.jet_calc) then 
            ! Initialize using the values from a previous run (if the file exists)
            inquire(file=trim(opath('quant_FF.dat')),exist=fexist)
            if (.not. fexist) then
                print *,'You first have to run trajectory and evolution modules at least once!'
                stop
            else
                open(30,file=trim(opath('quant_FF.dat')))
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
