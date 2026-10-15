module numerical_utils
    use constants, only: dp
    implicit none
    private

    public :: vector_log, integrate
    public :: locate, Interpol1d
    public :: safeLog10, safeLog, tridag_ser, zeros1D

contains

    subroutine vector_log(Vmin, Vmax, Vint, Vvec)
        real(dp), intent(in) :: Vmin, Vmax
        real(dp), intent(out) :: Vint
        real(dp), intent(out) :: Vvec(:)
        integer :: i

        Vint = (Vmax/Vmin)**(1.d0/real(size(Vvec)-1, dp))
        Vvec(1) = Vmin
        do i = 2, size(Vvec)
            Vvec(i) = Vvec(i-1)*Vint
        end do
    end subroutine vector_log


    subroutine integrate(N, H, FI, S)
        integer, intent(in) :: N
        real(dp), intent(in) :: H, FI(N)
        real(dp), intent(out) :: S
        real(dp) :: S0, S1, S2
        integer :: i

        S0 = 0.d0
        S1 = 0.d0
        S2 = 0.d0
        do i = 2, N-1, 2
            S1 = S1 + FI(i-1)
            S0 = S0 + FI(i)
            S2 = S2 + FI(i+1)
        end do
        S = H*(S1 + 4.d0*S0 + S2)/3.d0
        if (mod(N, 2) == 0) S = S + H*(5.d0*FI(N) + 8.d0*FI(N-1) - FI(N-2))/12.d0
    end subroutine integrate


    integer function locate(xx, x)
        real(dp), intent(in) :: xx(:), x
        integer :: jl, jm, ju, n

        n = size(xx)
        jl = 0
        ju = n + 1
        do while (ju-jl > 1)
            jm = (ju + jl)/2
            if ((xx(n) >= xx(1)) .eqv. (x >= xx(jm))) then
                jl = jm
            else
                ju = jm
            end if
        end do

        if (x == xx(1)) then
            locate = 1
        else if (x == xx(n)) then
            locate = n - 1
        else
            locate = jl
        end if
    end function locate


    real(dp) function Interpol1d(xp, x, v)
        real(dp), intent(in) :: xp, x(:), v(:)
        integer :: j_x, n
        real(dp) :: slope

        n = size(v)
        j_x = locate(x, xp)
        if (j_x == 0) then
            slope = (v(2) - v(1))/(x(2) - x(1))
            Interpol1d = max(v(1) + slope*(xp - x(1)), 0.d0)
        else if (j_x == n) then
            slope = (v(n) - v(n-1))/(x(n) - x(n-1))
            Interpol1d = max(v(n) + slope*(xp - x(n)), 0.d0)
        else if (v(j_x) > 0.d0 .and. v(j_x+1) > 0.d0) then
            slope = log(v(j_x+1)/v(j_x))/log(x(j_x+1)/x(j_x))
            Interpol1d = v(j_x)*(xp/x(j_x))**slope
        else
            slope = (v(j_x+1) - v(j_x))/(x(j_x+1) - x(j_x))
            Interpol1d = v(j_x) + slope*(xp - x(j_x))
        end if
    end function Interpol1d


    real(dp) function safeLog10(x)
        real(dp), intent(in) :: x

        if (x > 0.d0) then
            safeLog10 = log10(x)
        else
            safeLog10 = -500.d0
        end if
    end function safeLog10


    real(dp) function safeLog(x)
        real(dp), intent(in) :: x

        if (x > 0.d0) then
            safeLog = log(x)
        else
            safeLog = -500.d0
        end if
    end function safeLog


    subroutine an_error(message)
        character(len=*), intent(in) :: message

        write (*, "('ERROR ', A)") trim(message)
        error stop 'Program terminated by an_error'
    end subroutine an_error


    integer function assert_eq(nn, message)
        integer, intent(in) :: nn(:)
        character(len=*), intent(in) :: message

        if (all(nn(2:) == nn(1))) then
            assert_eq = nn(1)
        else
            call an_error('assert_eq failed with tag: '//trim(message))
        end if
    end function assert_eq


    subroutine tridag_ser(a, b, c, r, u)
        real(dp), intent(in) :: a(:), b(:), c(:), r(:)
        real(dp), intent(out) :: u(:)
        real(dp) :: gam(size(b)), bet
        integer :: n, j

        n = size(b)
        if (size(a) /= n-1 .or. size(c) /= n-1 .or. size(r) /= n .or. size(u) /= n) then
            call an_error('tridag_ser: inconsistent array sizes')
        end if
        bet = b(1)
        if (bet == 0.d0) call an_error('tridag_ser: error at code stage 1')
        u(1) = r(1)/bet
        do j = 2, n
            gam(j) = c(j-1)/bet
            bet = b(j) - a(j-1)*gam(j)
            if (bet == 0.d0) call an_error('tridag_ser: error at code stage 2')
            u(j) = (r(j) - a(j-1)*u(j-1))/bet
        end do
        do j = n-1, 1, -1
            u(j) = u(j) - gam(j+1)*u(j+1)
        end do
    end subroutine tridag_ser


    function zeros1D(n, small) result(values)
        integer, intent(in) :: n
        logical, intent(in) :: small
        real(dp) :: values(n)

        if (small) then
            values = 1.d-180
        else
            values = 0.d0
        end if
    end function zeros1D

end module numerical_utils
