program test_synchrotron
    use synchrotron
    implicit none
    real(8) :: E_ph, B, result_psyn, result_jsync, result_alpha, E_max, E_min, factor
    real(8), dimension(100) :: E, N_e
    integer :: k


    
    ! Define test values
    E_ph = 1.0d-16  ! Sample photon energy
    B = 1.0d0      ! Sample magnetic field
    ! Define range for geometrically spaced energy values
    E_min = 1.63742116e-06
    E_max = 8.18710578e-03
    factor = (E_max / E_min) ** (1.0d0 / (size(E) - 1))

    ! Generate geometrically spaced energy array
    E(1) = E_min
    do k = 2, size(E)
        E(k) = E(k-1) * factor
    end do
    N_e = E ** (-2.8d0)

    ! Call PSynchrotron function
    result_psyn = PSynchrotron(E_ph, E(1), B)
    print *, "PSynchrotron result: ", result_psyn

    ! Call jsync function
    result_jsync = jsync(E_ph, B, N_e, E)
    print *, "jsync result: ", result_jsync

    ! Call alpha_sync function
    result_alpha = alpha_sync(E_ph, B, N_e, E)
    print *, "alpha_sync result: ", result_alpha

end program test_synchrotron

