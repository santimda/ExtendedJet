MODULE map
  use globals
  implicit none
  private
  public :: map_sy

  contains

  !============================================================================
  ! map_sy (START)
  !============================================================================
  ! This subroutine transforms Lsy(z) into an emission map Lsy(x,y) 
  ! The input is the frequency in GHz at which to calculate the map
  !============================================================================

  subroutine map_sy(nu)

    real(dp), intent(in) :: nu ! [in GHz]
    real(dp) :: Esy(nIn), Lsy(nIn, Ns), dummy
    real(dp) :: E_nu, L_nu(Ns), log_Eint
    integer :: i_nu
    allocate( z(Ns), rj(Ns) )


    ! Read from a file: Esy(i,l), E*Lsy_unabs(i,l), E*Lsy_abs(i,l), l
    open(27, file='SEDsyn_cuts_ord.dat')
    do l = 1, Ns 
      do i = 1, nIn
        read(27,*) Esy(i), dummy, Lsy(i,l), dummy
        Lsy(i,l) = Lsy(i,l)/Esy(i)
      enddo
    enddo
    close(27)

    E_nu = h*nu*1.d9/eV ! Convert nu to the same units as Esy
    log_Eint = log(Esy(2)/Esy(1))
    i_nu = int(log(E_nu/Esy(1))/log_Eint) + 1 

    do l = 1, Ns
      L_nu(l) = Lsy(i_nu,l)
    enddo

    ! TO BE COMPLETED: 
    ! from z(l), r(l), L_nu(l) ---> (x, y, S_nu(x,y))
    ! UNITS: [x,y] = arcsec, [S_nu] = mJy ????
    ! STORE IT IN map_sy.dat 

    deallocate( z, rj )

  end subroutine map_sy

  !============================================================================
  ! map_sy (END)
  !============================================================================

  !============================================================================
  ! map_sy_conv (START)
  !============================================================================
  ! This subroutine convolves the emission map Lsy(x,y) wit the 
  ! response of a gaussian beam to calculate a synthetic map
  !============================================================================

  subroutine map_sy_conv()

    real(dp), parameter :: f_Ny = 0.6
    integer, parameter :: ml_tot = 2*Ns   ! Cantidad total de datos ????
    !real(dp), parameter :: coc = sqrt(8.d0*log(2.d0))       ! difference between sigma and theta!
    !real(dp), parameter :: theta_x = sigma_x*coc, theta_y = sigma_y*coc ! HPBW = Resolution
    real(dp), parameter :: sigma_x = 58.d0, sigma_y = 57.d0  ! ????
    real(dp), dimension(:), allocatable :: x_c, y_c
    real(dp), dimension(:,:), allocatable :: F_map_abs
    real(dp) :: x_map(ml_tot), y_map(ml_tot), F_abs(ml_tot)
    real(dp) :: ancho_x, ancho_y, D_x, D_y, sjo
    real(dp) :: x_c0, y_c0, d_pesada2, suma, D_f, atenua, sum, mJy_beam
    real(dp) :: A_beam, A_cell, dummy
    integer :: m_x, m_y

    A_beam = 2*pi*sigma_x*sigma_y     ! = Gaussian normalization factor 
    A_cell = f_Ny**2*sigma_x*sigma_y  ! cell sample size
    
    sjo = abs(sin(jo))  ! CHECK IF SINE OF COSINE

    !-------------------------------------------------------
    ! Inicializo los datos de la región donde convolucionaré

    ancho_x = zmax * thetaj * 1.2d0 ! [mas](extra margin)
    ancho_y = zmax * sjo * 1.2d0    ! [mas](extra margin)

    D_x = f_Ny*sigma_x
    D_y = f_Ny*sigma_y
    m_x = int(ancho_x/D_x) ! Number of ellipses in the x-direction
    m_y = int(ancho_y/D_y) ! Number of ellipses in the y-direction
    x_c0 = -ancho_x/2.d0   ! Center of the first ellipse (bottom left corner)
    y_c0 = -ancho_y/4.d0
    allocate( x_c(m_x), y_c(m_y), F_map_abs(m_x,m_y)  )
    !    write(*,*) m_x,m_y


    !------------------------------------------------------
    ! Data reading

    sum = 0.d0
    open(13, file='map_sy.dat')
    do l = 1, ml_tot
      read(13,*) x_map(l), y_map(l), F_abs(l)
      sum = sum + F_abs(l)
    end do
    write(*,*) 'Pre-conv: S(GHz)=',sum,'mJy'
    close(13)

   !------------------------------------------

    do i = 1,m_x
      x_c(i) = x_c0 + (i-0.5)*D_x  ! x-coord of the center of the i-th ellipse

      do j = 1,m_y
        y_c(j) = y_c0 + (j-0.5)*D_y  ! y-coord of the center of the (x_i,y_j) ellipse

        suma = 0.d0  ! simulate noise with random number?
        do l = 1, ml_tot
          d_pesada2 = ((x_map(l)-x_c(i))**2)/(2.d0*sigma_x**2) + ((y_map(l)-y_c(j))**2)/(2.d0*sigma_y**2) ! added extra sigma
          if (d_pesada2 < 6.d0) then   
            ! Gaussian beam -> /(sigma*sqrt(2*pi)), both in x and y
            atenua = exp( -d_pesada2 ) !/ norm
            D_F = F_abs(l) * atenua   
          else
            D_F = 0.d0 ! flux drops < 1e-3 farther from there, no need to sum
          end if
          suma = suma + D_f
        end do

        F_map_abs(i,j) = suma

      end do

    end do

    !--------------------------------------------------------------
    ! Write data and do consistency check

    mJy_beam = A_cell/A_beam ! Transform mJy/beam -> mJy (for checks)

    sum = 0.d0
    open (17, file='map_sy_conv.dat')
    do i = 1, m_x
      do j = 1, m_y
        write(17,61) x_c(i), y_c(j), F_map_abs(i,j)
        sum = sum + F_map_abs(i,j)*mJy_beam
        !print*, sum, F_map_abs(i,j), mJy_beam
      end do
      write(17,*) 
    end do
    write(*,*) 'Post-conv: S(GHz)=',sum,'mJy'
    close(17)

    61  format(F11.4,X,F11.4,X,E11.4)

    deallocate( x_c, y_c, F_map_abs )

  end subroutine map_sy_conv

  !============================================================================
  ! map_sy_conv (END)
  !============================================================================


END MODULE map
