MODULE trajectory
  use globals
  implicit none
  private
  public :: trajectory_calc

  contains

  !============================================================================
  ! trajectory_calc (START)
  !============================================================================
  ! This subroutine calculates a vertical jet trajectory. The input is the '' (jet) 
  ! or 'cj' (counter-jet). 
  !============================================================================

  subroutine trajectory_calc(jet)

    character(LEN=*), intent(in) :: jet ! 'cj' or 'jet'
    real(dp) :: SO(3),Px,Py,Pz,Ptot,y0,z0
    allocate (x(Ns),y(Ns),z(Ns),rj(Ns),lenj(Ns),d(Ns),jw(Ns),wo(Ns),jo(Ns),isOccult(Ns))

    write(*,'(X,A,X,I4)') 'Number of segments:',Ns

    !System-observer vector (unitary)
    SO = (/sin(inc*pi/180d0)*cos(azi*pi/180d0),sin(inc*pi/180)*sin(azi*pi/180d0),cos(inc*pi/180d0)/)

    !Initial conditions
    x(1) = a*cos(orbPhase*pi/180d0)
    y(1) = a*sin(orbPhase*pi/180d0)
    z(1) = zrec
    Px = 0d0      !Blob momentum projection in the x direction
    Py = 0d0      !Blob momentum projection in the y direction
    Pz = Pdotj*dt !Blob momentum projection in the z direction
    lenj(1) = zrec

    do i = 1, Ns

      Ptot = sqrt(Px**2+Py**2+Pz**2)        !Total momentum of the jet      
      d(i) = sqrt(x(i)**2+y(i)**2+z(i)**2)  !Distance to the star
      rj(i) = lenj(i)*tan(thetaj)           !Jet radius     
      jo(i) = acos((Px*SO(1)+Py*SO(2)+Pz*SO(3))/Ptot)     !Jet-observer angle

      !Jet trajectory
      if (i /= Ns) then
        x(i+1) = x(i)
        y(i+1) = y(i)
        if (jet == 'jet') then
          z(i+1) = z(i) + vj*dt
        elseif (jet == 'cj') then
          z(i+1) = z(i) - vj*dt
        endif
        lenj(i+1) = lenj(i) + vj*dt
      endif

      !Coordinates at x=0 of a photon emitted at (x(i),y(i),z(i)) and directed towards the observer (valid
      !only for SO(1) /= 0)
      y0 = y(i) - SO(2)/SO(1)*x(i)
      z0 = z(i) - SO(3)/SO(1)*x(i)

    enddo

    open(10, file='trajectory'//jet//'.dat')
    do i = 1, Ns
      write(10,*) x(i),y(i),z(i),rj(i),lenj(i),d(i),jo(i)
    enddo
    close(10)

    deallocate(x,y,z,rj,lenj,d,jo)
    return
  end subroutine trajectory_calc


  !============================================================================
  ! trajectory_calc (END)
  !============================================================================


END MODULE trajectory
