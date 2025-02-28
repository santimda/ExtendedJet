module mathematics
    use constants
contains

    real(dp) function safeLog10(x)
        real(dp), intent(in) :: x
        if (x > 0d0) then
            safeLog10 = log10(x)
        else
            safeLog10 = -500
        endif
        return
    end function safeLog10

    real(dp) function safeLog(x)
        real(dp), intent(in) :: x
        if (x > 0d0) then
            safeLog = log(x)
        else
            safeLog = -500
        endif
        return
    end function safeLog

    ! :::: This produces an error message ::::
    subroutine an_error(string)
        implicit none
        character(len=*), intent(IN) :: string
        write (*, "('ERROR ', A)") string
        stop 'program terminated by an_error'
    end subroutine an_error

    function assert_eq(nn, string)
        implicit none
        character(len=*), intent(in) :: string
        integer, intent(in), dimension(:) :: nn
        integer :: assert_eq
        if (all(nn(2:) == nn(1))) then
           assert_eq = nn(1)
        else
           write(*,*) 'ERROR: an assert_eq failed with this tag:', string
           stop 'program terminated by assert_eq'
        end if
    end function assert_eq

    subroutine tridag_ser(a,b,c,r,u)
        !  Description:
        !     Solver of a tridiagonal matrix. Based on the code in
        !     "Numberical Recipes".
        implicit none
        real(dp),dimension(:),intent(in) :: a,b,c,r
        real(dp),dimension(:),intent(out) :: u
        real(dp),dimension(size(b)) :: gam
        integer :: n,j
        real(dp) :: bet
        n=assert_eq((/ size(a)+1,size(b),size(c)+1,size(r),size(u) /),'tridag_ser')
        bet=b(1)
        if (bet==0.0d0) call an_error('tridag_ser: error at code stage 1')
        u(1)=r(1)/bet
        do j=2,n
           gam(j)=c(j-1)/bet
           bet=b(j)-a(j-1)*gam(j)
           if (bet==0.0d0) call an_error('tridag_ser: error at code stage 2')
           u(j)=(r(j)-a(j-1)*u(j-1))/bet
        end do
        do j=n-1,1,-1
           u(j)=u(j)-gam(j+1)*u(j+1)
        end do
    end subroutine tridag_ser

    function zeros1D(n, small) result(a)
        implicit none
        integer, intent(in) :: n
        logical, intent(in) :: small
        real(dp), dimension(n) :: a
        if (small) then
           a = 1d-200
        else
           a = 0d0
        end if
    end function zeros1D

end module
