module synchrotron
    implicit none
    real(8), parameter :: pi = 3.141592653589793
    real(8), parameter :: qe = 4.80320425d-10  ! Electron charge [esu]
    real(8), parameter :: h = 6.62607015d-27  ! Planck constant [erg·s]
    real(8), parameter :: c = 2.99792458d10   ! Speed of light [cm/s]
    real(8), parameter :: me = 9.10938356d-28 ! Electron mass [g]
    
contains
    real(8) function PSynchrotron(E_ph, E_e, B) result(resultado)
    implicit none
        real(8), intent(in) :: E_ph, B, E_e
        real(8) :: K
        real(8) :: x, R
        real(8) :: Ec
    
    ! Define constants
    K = sqrt(3.0d0) * qe**3 / (h * me * c**2)
    
    Ec = (3.0d0 / (4.0d0 * pi)) * (qe * h * B) / (me**3 * c**5) * E_e**2
        ! Calculate x and R for each element of E
        x = E_ph / Ec
        R = (1.808 * x**(1.0d0 / 3.0d0) / sqrt(1.0d0 + 3.4 * x**(2.0d0 / 3.0d0))) * &
                ((1.0d0 + 2.21 * x**(2.0d0 / 3.0d0) + 0.347 * x**(4.0d0 / 3.0d0)) / &
                (1.0d0 + 1.353 * x**(2.0d0 / 3.0d0) + 0.217 * x**(4.0d0 / 3.0d0))) * exp(-x)

    ! Compute the emitted power for each element of E
    resultado = K * B * R

    end function PSynchrotron

    !****************************************************************************
    ! Synchrotron emissivity for a specific photon energy
    !****************************************************************************
    function jsync(E_ph, B, N, E) result(integral)
    implicit none
    integer :: nsize  ! Number of elements in arrays
    real(8), intent(in) :: E_ph, B
    real(8), dimension(:), intent(in) :: E
    real(8), dimension(:), intent(in) :: N
    real(8) :: integral
    real(8) , dimension(:), allocatable :: integ
    integer :: i
    nsize=size(E)
    allocate(integ(nsize))
    ! The values inside the integral
    do i=1, nsize
        integ(i)=PSynchrotron(E_ph, E(i), B)
    enddo
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
    real(8), dimension(:) :: E, N
    real(8), dimension(:), ALLOCATABLE :: log_eVec, log_neVec
    real(8), intent(in) :: B, E_ph
    real(8) :: logne
    real(8) :: alpha

    allocate(log_eVec(size(E)), log_neVec(size(N)))
    log_eVec = log(E)
    log_neVec = log(N)

    alpha = quad(integrand, log(E(1)*1.1d0), log(E(size(E))/1.1d0), E_ph, B, log_eVec, log_neVec)
    alpha = -h**3*c**2/ (8.*pi*E_ph*E_ph)*alpha
    deallocate(log_eVec, log_neVec)

    end function alpha_sync

    !Function that calculates the interpolation
    real(8) function log_interp1d(x, xp, yp)
    real(8) :: x
    real(8), dimension(:) :: xp, yp
    integer :: i, j, n
    real(8) :: t
    
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
    
    real(8) :: df
    real(8), dimension(:) :: log_eVec, log_neVec
    real(8) :: loge
    real(8) :: E_e, integrand
    real(8), dimension(1) :: E_e_int 
    real(8), intent(in) :: E_ph, bField
    real(8), dimension(1) :: p

    E_e=exp(loge)
    
    call compute_derivative(log_interp1d, loge, 1d-4, df, log_eVec,log_neVec)
    if (df.lt. 0. .and. df.gt.-10.) then
        integrand=E_e * PSynchrotron(E_ph, E_e, bField) * E_e * (exp(log_interp1d(loge,log_eVec,log_neVec))/(E_e*E_e)) * (df-2)
    else
        integrand=0.0
    end if

    end function integrand

    !Integration function
    function quad(f, a, b, E_ph, bField , log_eVec,log_neVec) result(q)
        implicit none
        real(8), intent(in) :: a, b, E_ph, bField
        interface
            real(8) function f(x, E_ph, bField, log_eVec, log_neVec)
            real(8) :: x
            real(8), dimension(:) :: log_eVec, log_neVec
            real(8), intent(in) :: E_ph, bField
        end function f
        end interface
        real(8), dimension(:) :: log_eVec, log_neVec
        real(8) :: q
        real(8) :: h, x
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
    real(8), intent(in) :: x, h
    interface
        real(8) function f(x,log_eVec, log_neVec)
        real(8) :: x
        real(8), dimension(:) :: log_eVec, log_neVec
      end function f
    end interface
    real(8), dimension(:) :: log_eVec, log_neVec
    real(8),intent(out) :: df
    df = (f(x + h,log_eVec,log_neVec) - f(x - h,log_eVec,log_neVec)) / (2.0 * h)
    end subroutine compute_derivative
end module