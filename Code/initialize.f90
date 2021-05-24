MODULE initialize
    use globals
    use controls
    implicit none
    
    contains 

    !============================================================================
    ! initialize_allocate (START)
    !============================================================================
    subroutine initialize_allocate
    !============================================================================
    ! allocate and initialize vectors 
    !============================================================================
        logical :: fexist

        allocate( z(Ns), rj(Ns) )
        allocate( BFF(Ns), vnorm(Ns), VFF(Ns), EphMFF(Ns), nfFF(Ns), Kappa(Ns), sjoFF(NS) ) 
        allocate( deltaZ(Ns), deltatFF(Ns) )

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
