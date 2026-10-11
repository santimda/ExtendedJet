program test_numerical_utilities
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use constants, only: dp
    use numerical_utils, only: vector_log, integrate, locate, Interpol1d, &
                               safeLog10, safeLog, tridag_ser, zeros1D
    implicit none

    real(dp), parameter :: tolerance = 1.d-12
    real(dp) :: values(5), expected_values(5), ratio
    real(dp) :: x_odd(5), f_odd(5), integral_odd
    real(dp) :: x_even(4), f_even(4), integral_even
    real(dp) :: increasing(4), decreasing(4)
    real(dp) :: x_interp(2), v_interp(2)
    real(dp) :: lower(2), diagonal(3), upper(2), rhs(3), solution(3)
    real(dp) :: zero_values(4), small_values(4)

    call vector_log(1.d0, 1.d4, ratio, values)
    expected_values = (/ 1.d0, 10.d0, 100.d0, 1000.d0, 10000.d0 /)
    call assert_array_close('vector_log values', values, expected_values, tolerance)
    call assert_close('vector_log ratio', ratio, 10.d0, tolerance)

    x_odd = (/ 0.d0, 1.d0, 2.d0, 3.d0, 4.d0 /)
    f_odd = x_odd**2
    call integrate(5, 1.d0, f_odd, integral_odd)
    call assert_close('integrate odd points', integral_odd, 64.d0/3.d0, tolerance)

    x_even = (/ 0.d0, 1.d0, 2.d0, 3.d0 /)
    f_even = x_even**2
    call integrate(4, 1.d0, f_even, integral_even)
    call assert_close('integrate even points', integral_even, 9.d0, tolerance)

    increasing = (/ 1.d0, 2.d0, 4.d0, 8.d0 /)
    call assert_equal('locate below increasing range', locate(increasing, 0.5d0), 0)
    call assert_equal('locate increasing interval', locate(increasing, 3.d0), 2)
    call assert_equal('locate increasing first endpoint', locate(increasing, 1.d0), 1)
    call assert_equal('locate increasing last endpoint', locate(increasing, 8.d0), 3)
    call assert_equal('locate above increasing range', locate(increasing, 9.d0), 4)

    decreasing = (/ 8.d0, 4.d0, 2.d0, 1.d0 /)
    call assert_equal('locate decreasing interval', locate(decreasing, 3.d0), 2)
    call assert_equal('locate above decreasing range', locate(decreasing, 9.d0), 0)
    call assert_equal('locate below decreasing range', locate(decreasing, 0.5d0), 4)

    x_interp = (/ 1.d0, 10.d0 /)
    v_interp = (/ 10.d0, 1000.d0 /)
    call assert_close('log interpolation', Interpol1d(sqrt(10.d0), x_interp, v_interp), &
                      100.d0, tolerance)
    x_interp = (/ 1.d0, 2.d0 /)
    v_interp = (/ 0.d0, 10.d0 /)
    call assert_close('linear interpolation', Interpol1d(1.5d0, x_interp, v_interp), &
                      5.d0, tolerance)
    v_interp = (/ 10.d0, 20.d0 /)
    call assert_close('non-negative extrapolation', Interpol1d(0.d0, x_interp, v_interp), &
                      0.d0, tolerance)

    call assert_close('safeLog10 positive', safeLog10(100.d0), 2.d0, tolerance)
    call assert_close('safeLog10 zero', safeLog10(0.d0), -500.d0, tolerance)
    call assert_close('safeLog positive', safeLog(exp(2.d0)), 2.d0, tolerance)
    call assert_close('safeLog negative', safeLog(-1.d0), -500.d0, tolerance)

    lower = -1.d0
    diagonal = 2.d0
    upper = -1.d0
    rhs = (/ 1.d0, 0.d0, 1.d0 /)
    call tridag_ser(lower, diagonal, upper, rhs, solution)
    call assert_array_close('tridiagonal solver', solution, (/ 1.d0, 1.d0, 1.d0 /), tolerance)

    zero_values = zeros1D(4, .false.)
    small_values = zeros1D(4, .true.)
    call assert_array_close('zeros1D zeros', zero_values, (/ 0.d0, 0.d0, 0.d0, 0.d0 /), tolerance)
    call assert_array_close('zeros1D small', small_values, &
                            (/ 1.d-180, 1.d-180, 1.d-180, 1.d-180 /), 0.d0)

    print *, 'Numerical utility tests passed.'

contains

    subroutine assert_close(name, actual, expected, allowed_error)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, allowed_error

        if (.not. ieee_is_finite(actual) .or. abs(actual-expected) > allowed_error) then
            print *, 'FAILED: ', trim(name)
            print *, 'actual, expected = ', actual, expected
            error stop 1
        end if
    end subroutine assert_close


    subroutine assert_array_close(name, actual, expected, allowed_error)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual(:), expected(:), allowed_error

        call assert_close(name, maxval(abs(actual-expected)), 0.d0, allowed_error)
    end subroutine assert_array_close


    subroutine assert_equal(name, actual, expected)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, 'FAILED: ', trim(name)
            print *, 'actual, expected = ', actual, expected
            error stop 1
        end if
    end subroutine assert_equal

end program test_numerical_utilities
