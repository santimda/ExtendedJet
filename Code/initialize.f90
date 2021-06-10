MODULE initialize
    use globals
    use constants
    use controls
    implicit none
    
    contains 

    !============================================================================
    ! initialize_allocate (START)
    !============================================================================
    subroutine initialize_allocate(fB_value, Lj_value)
    !============================================================================
    ! allocate and initialize vectors 
    !============================================================================
        real(dp), intent(in) :: fB_value, Lj_value
        logical :: fexist

        allocate( z(Ns), rj(Ns) )
        allocate( BFF(Ns), vnorm(Ns), VFF(Ns), EphMFF(Ns), nfFF(Ns), Kappa(Ns), sjoFF(NS) ) 
        allocate( deltaZ(Ns), deltatFF(Ns) )

        ! Initialize values given by free parameters:
   
        Lj = Lj_value 
        fnT = 1.d-2              !Nonthermal energy fraction
        LnT = Lj*fnT             !Available nonthermal luminosity
        LnTFF = LnT!/gj**2       !Available nonthermal luminosity in the FF

        cte_inj = LnTFF/( (zmax_inj**(-a_inj+1.) - z0**(-a_inj+1.))/(-a_inj+1.) )
        
        fB = fB_value
        B0 = sqrt(fB*8d0*Lj/(vj*(thetaj*z0)**2))

        ! Initialize arrays
        if (.not.jet_calc) then 
            ! Initialize using the values from a previous run (if the file exists)
            inquire(file='quant_FF.dat',exist=fexist)
            if (.not. fexist) then
                print*,'You first have to run trajectory and evolution modules at least once!'
                stop
            else
                open(30,file='quant_FF.dat')  
                do i = 1, Ns
                    !Quantities at each segment
                    read(30,*) z(i), BFF(i), rj(i), vnorm(i), Kappa(i), sjoFF(i), VFF(i), EphMFF(i), nfFF(i), deltaZ(i), deltatFF(i)
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

        deallocate( z, rj )
        deallocate( BFF, vnorm, VFF, EphMFF, nfFF, Kappa, sjoFF )
        deallocate( deltaZ, deltatFF )

    end subroutine initialize_deallocate
    !============================================================================
    ! initialize_deallocate (END)
    !============================================================================


END MODULE initialize
