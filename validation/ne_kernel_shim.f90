! Test-only f2py shim: explicit-size interface to edistribution::FP_FinDif_cool, whose assumed-shape
! intent(out) argument f2py cannot size. Built by validation/ne_validation.py together with the production
! sources, so the kernel under test is the one in Code/edistribution.f90 (not a copy).
subroutine fp_steps(n, nsteps, dt, g, n_in, gdot, qq, tesc, n_out)
    use edistribution, only: FP_FinDif_cool
    implicit none
    integer, intent(in) :: n, nsteps
    double precision, intent(in) :: dt, tesc
    double precision, intent(in) :: g(n), n_in(n), gdot(n), qq(n)
    double precision, intent(out) :: n_out(n)
    double precision :: work(n)
    integer :: i

    work = n_in
    do i = 1, nsteps
        call FP_FinDif_cool(dt, g, work, n_out, gdot, qq, tesc)
        work = n_out
    end do
end subroutine fp_steps
