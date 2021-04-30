MODULE opac_gg
  use controls
  use globals
  implicit none
  
contains

  !============================================================================
  ! opac_gg_calc (jet)
  !============================================================================

  subroutine opac_gg_calc ()
    real*8 :: dummy,ang_int,ang_min
    real*8 :: Eg(NSED_IC),LIC_unabs(NSED_IC),LIC_abs(NSED_IC)
    real*8 :: LIC_tot_abs(NSED_IC),LIC_tot_unabs(NSED_IC),LIC1(NSED_IC),LIC2(NSED_IC)
    integer :: iEf,iang
    allocate(tau_gg(mE_gg,m_ang_gg),lenj(Ns),wo(Ns),d(Ns)) 

    if (remake_table_gg) then
    ! Calculate a matrix with the values of tau at a distance of 1 Rstar
    ! for different wind-observer angles and gamma-ray photon energies Eg,
    ! such that: tau(Eg,ang)

       write(*,'(3X,A)') 'Calculate gamma-gamma opacity'
       call tau_calc(tau_gg)

       open (88,file='opac_gg.dat')
       do i=1,mE_gg
          write(88,*) (tau_gg(i,j), j=1,m_ang_gg)
       enddo
       close(88)
       
    else
       open (88,file='opac_gg.dat')
       do i=1,mE_gg
          read(88,*) (tau_gg(i,j), j=1,m_ang_gg)
       enddo
       close(88)
    endif
    Ef_int_gg = (Ef_max_gg/Ef_min_gg)**(1.d0/float(mE_gg-1))
    
    ang_int = pi/float(m_ang_gg)
    ang_min = ang_int/2.d0
    
    
    !-------------------------------------------------------------------
    ! READ QUANTITIES
    !-------------------------------------------------------------------
    
    open(30,file='Quant_FF.dat')
    open(11,file='Trajectory.dat')
    do j = 1,Ns
       read(30,*) lenj(j),dummy,dummy,dummy,dummy,dummy,dummy,dummy,wo(j)
       read(11,*) dummy,dummy,dummy,dummy,dummy,d(j)
    enddo
    close(30)
    close(11)

    !-------------------------------------------------------------------
    ! CORRECT IC SED
    !-------------------------------------------------------------------
    
    LIC_tot_unabs = 0.d0
    LIC_tot_abs = 0.d0
    open(40,file='SEDIC_cuts.dat')
    open(41,file='SEDIC_abs_cuts.dat')
    do j = 1,Ns
       
       do i = 1,NSED_IC
          read(40,*) Eg(i),LIC_unabs(i),dummy
       enddo
       Eg = Eg*eV
       if (j < Ns) read(40,*)
       
       iang = int(m_ang_gg*abs(wo(j)-ang_min)/pi)+1
       do i=1,NSED_IC
          LIC_tot_unabs(i) = LIC_tot_unabs(i)+LIC_unabs(i)
          if (Ef_min_gg < Eg(i) .and. Eg(i) < Ef_max_gg) then
             iEf = int(log(Eg(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
             LIC_abs(i) = LIC_unabs(i) * exp( -(Rstar/d(j))*tau_gg(iEf,iang) )
             LIC_tot_abs(i) = LIC_tot_abs(i) + LIC_abs(i)
          else
             LIC_abs(i) = LIC_unabs(i)
             LIC_tot_abs(i) = LIC_tot_unabs(i)
          endif
          
       	  !Jet contribution of the segments up to given lengths
       	  if (lenj(j)-lenj(1) <= 0.1d0*a) LIC1(i) = LIC_tot_abs(i)
       	  if (lenj(j)-lenj(1) <= a) LIC2(i) = LIC_tot_abs(i)
       	  

          write(41,*) Eg(i)/eV,LIC_unabs(i),LIC_abs(i),j
       enddo
       write(41,*)
    enddo
	close(40)
    close(41)

    open(42,file='SEDIC_abs_jet.dat')
    do i = 1,NSED_IC
       write(42,*) Eg(i)/eV,LIC1(i),LIC2(i),LIC_tot_unabs(i),LIC_tot_abs(i)     	  
    enddo
    close(42)
    
    deallocate(tau_gg,lenj,wo,d)
  END subroutine opac_gg_calc
  
  !============================================================================
  ! opac_gg_calc (counter-jet)
  !============================================================================

  subroutine opac_gg_calcCJ ()
    real*8 :: dummy,ang_int,ang_min
    real*8 :: Eg(NSED_IC),LIC_unabs(NSED_IC),LIC_abs(NSED_IC)
    real*8 :: LIC_tot_abs(NSED_IC),LIC_tot_unabs(NSED_IC),LIC1(NSED_IC),LIC2(NSED_IC)
    integer :: iEf,iang
    allocate(tau_gg(mE_gg,m_ang_gg),lenj(Ns),wo(Ns),d(Ns)) 

    open (88,file='opac_gg.dat')
    do i=1,mE_gg
    	read(88,*) (tau_gg(i,j), j=1,m_ang_gg)
    enddo
    close(88)
    
    Ef_int_gg = (Ef_max_gg/Ef_min_gg)**(1.d0/float(mE_gg-1))    
    ang_int = pi/float(m_ang_gg)
    ang_min = ang_int/2.d0    
    
    !-------------------------------------------------------------------
    ! READ QUANTITIES
    !-------------------------------------------------------------------
    
    open(30,file='Quant_FFCJ.dat')
    open(11,file='TrajectoryCJ.dat')
    do j = 1,Ns
       read(30,*) lenj(j),dummy,dummy,dummy,dummy,dummy,dummy,dummy,wo(j)
       read(11,*) dummy,dummy,dummy,dummy,dummy,d(j)
    enddo
    close(30)
    close(11)

    !-------------------------------------------------------------------
    ! CORRECT IC SED
    !-------------------------------------------------------------------
    
    LIC_tot_unabs = 0.d0
    LIC_tot_abs = 0.d0
    open(40,file='SEDIC_cutsCJ.dat')
    open(41,file='SEDIC_abs_cutsCJ.dat')
    do j = 1,Ns
       
       do i = 1,NSED_IC
          read(40,*) Eg(i),LIC_unabs(i),dummy
       enddo
       Eg = Eg*eV
       if (j < Ns) read(40,*)
       
       iang = int(m_ang_gg*abs(wo(j)-ang_min)/pi)+1
       do i=1,NSED_IC
          LIC_tot_unabs(i) = LIC_tot_unabs(i)+LIC_unabs(i)
          if (Ef_min_gg < Eg(i) .and. Eg(i) < Ef_max_gg) then
             iEf = int(log(Eg(i)/Ef_min_gg)/log(Ef_int_gg)) + 1
             LIC_abs(i) = LIC_unabs(i) * exp( -(Rstar/d(j))*tau_gg(iEf,iang) )
             LIC_tot_abs(i) = LIC_tot_abs(i) + LIC_abs(i)
          else
             LIC_abs(i) = LIC_unabs(i)
             LIC_tot_abs(i) = LIC_tot_unabs(i)
          endif
          
       	  !Jet contribution of the segments up to given lengths
       	  if (lenj(j)-lenj(1) <= 0.1d0*a) LIC1(i) = LIC_tot_abs(i)
       	  if (lenj(j)-lenj(1) <= a) LIC2(i) = LIC_tot_abs(i)       	  

          write(41,*) Eg(i)/eV,LIC_unabs(i),LIC_abs(i),j
       enddo
       write(41,*)
    enddo
    close(40)
    close(41)

    open(42,file='SEDIC_abs_CJ.dat')
    do i = 1,NSED_IC
       write(42,*) Eg(i)/eV,LIC1(i),LIC2(i),LIC_tot_unabs(i),LIC_tot_abs(i) 	  
    enddo
    close(42)
    
    deallocate(tau_gg,lenj,wo,d)
  END subroutine opac_gg_calcCJ


  !--------------------------------------------------------------------
  !
  !     Subroutines
  !
  !--------------------------------------------------------------------

  subroutine tau_calc(tau_table)
    integer, parameter :: dp = kind(1.d0)
    integer, parameter :: ir_max = 100
    integer, parameter :: ie0_max = 20
    real(dp),intent(out) :: tau_table(mE_gg,m_ang_gg)
    integer :: iang,ie,ir,ie0    
    real(dp) :: Cnumc2,tau,Ef(mE_gg),dang,ang,xi,yi,dr
    real(dp) :: angi,mu,nu,Bb,Be,dEps_int,EpsM,Eps0,Eps0_min_factor,Eps0_min,Eps0_max,Eps_int
    real(dp) :: cte_sigma_gg,Nbmc2,Nst,xe,sigma_gg,sum
    real(dp) :: f(ir_max),cang,sang

    ! Ef_min_gg and Ef_max_gg are defined in globals, and also mE_gg
    ! tau_table is where the opacity matrix is stored and returned

    Cnumc2 = 1.23d20
    cte_sigma_gg = (0.5d0*pi*(qe**2/mec2)**2)
    EpsM = EphM/mec2

    call vector_log(Ef_min_gg,Ef_max_gg,Ef_int_gg,Ef)
    Ef = Ef/mec2 ! Adimensional

    dang = pi/m_ang_gg
    dr = Rstar/15.d0   ! step for integration R_final ~ Rstar + Rstar*ir_max/15

    do ie=1,mE_gg       ! loop in gamma-ray photon energy: Ef(ie)
       ang = dang/2.d0  ! initial angle

       Eps0_min_factor = 1.d0/Ef(ie)/0.5d0
       do iang=1,m_ang_gg ! loop in interaction angle: ang
          cang = cos(ang)
          sang = sin(ang)
          xi = Rstar        ! initial position
          yi = 0.d0
          tau = 0.d0

          do ir=1,ir_max  ! loop in gamma-ray photon trayectory: xi,yi
             sum = 0.d0    
             xi = dr*cang+xi
             yi = dr*sang+yi
             mu = (xi*cang+yi*sang)/sqrt(xi**2+yi**2)
             angi=acos(mu)! interaction angle along the path integral; i.e. acos e_*e_gamma
             !             cangi = cos(angi)
             Eps0_min = Eps0_min_factor/(1.d0-mu)
             Eps0_max = 20.d0*EpsM
             Eps_int = (Eps0_max/Eps0_min)**(1.d0/float(ie0_max))
             if (Eps0_min >= Eps0_max) goto 0111 !no absorption
             Eps0 = Eps0_min*Eps_int**0.5d0  

             do ie0=1,ie0_max ! loop in stellar photon energies, Eps0
                dEps_int = Eps0*(Eps_int-1.d0)

                ! Densidad de fotones target
                nu = Eps0*Cnumc2
                Bb = (2.d0*h*nu**3/c**2)/(exp(h*nu/(kB*Tstar))-1.d0)
                nbmc2 = Bb*Cnumc2/(nu*h)
                nst = (nbmc2*pi/c)*(Rstar**2/(xi**2+yi**2)) 

                ! sigma gamma-gamma (Gould & Schreder 1967)
                xe = Ef(ie)*Eps0*0.5d0*(1.d0-mu)
                Be = (1.d0-1./xe)**0.5d0
                sigma_gg = cte_sigma_gg*(1.d0-Be**2) &
                     *((3.d0-Be**4)*log((1+Be)/(1-Be))-2.d0*Be*(2.d0-Be**2))
                if (xe < 1.d0) sigma_gg=0.d0

                ! Opacidad: tau(Efe,ang)
                sum = sum+sigma_gg*nst*(1.d0-mu)*dEps_int
                ! tau=tau+sigma_gg*nst*dr*(1.d0-cangi)*deint
                Eps0 = Eps0*Eps_int
             enddo
0111         continue
             f(ir)=sum
          enddo
          call integrate(ir_max,dr,f,tau)
          tau_table(ie,iang)=tau

          ang=ang+dang

       enddo

    enddo

    return

  end subroutine tau_calc

  !======================================================================

END MODULE opac_gg
