MODULE opac_ff
  !This module calculates the f-f absorption and emission of the unshocked SW 
  use globals
  use controls
  implicit none
  public :: opac_ff_calc

contains

  !============================================================================
  ! opac_ff_calc (jet)
  !============================================================================

  subroutine opac_ff_calc ()
    real*8 :: dummy,ang_int,ang_min
    real*8 :: Eg(NSED_sy),Lsy_noabs(NSED_sy),Lsy_unabs(NSED_sy),Lsy_abs(NSED_sy)
    real*8 :: Lsy_tot_abs(NSED_sy),Lsy_tot_noabs(NSED_sy),Lsy_tot_unabs(NSED_sy),Lsy1(NSED_sy),Lsy2(NSED_sy)
    real*8 :: Ef_int_ff
    integer :: iEf,iang
    allocate(tau_ff(mE_ff,m_ang_ff),lenj(Ns),wo(Ns),d(Ns)) 

    if (remake_table_ff) then
       ! Calculate a matrix tau(Eg,ang) with the values of tau at a distance "a"
       ! for different wind-observer angles and emitted photon energies Eg

       write(*,'(3X,A)') 'Calculate free-free opacity'
       call tau_ff_calc(tau_ff)
       
       open (88,file='opac_ff.dat')
       do i=1,mE_ff
          write(88,*) (tau_ff(i,j), j=1,m_ang_ff)
       enddo
       close(88)

    else
       
       open (88,file='opac_ff.dat')
       do i=1,mE_ff
          read(88,*) (tau_ff(i,j), j=1,m_ang_ff)
       enddo
       close(88)
       
    endif

    Ef_int_ff = (Ef_max_ff/Ef_min_ff)**(1.d0/float(mE_ff-1))
    ang_int = pi/float(m_ang_ff)
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
    ! CORRECT SYNCHROTRON SED
    !-------------------------------------------------------------------
    
    Lsy_tot_unabs = 0.d0
    Lsy_tot_abs = 0.d0
    Lsy_tot_noabs = 0.d0
    open(40,file='SEDsyn_cuts_ord.dat')
    open(41,file='SEDsyn_abs_cuts.dat')
    do j = 1,Ns

       do i = 1,NSED_sy
          read(40,*) Eg(i),Lsy_noabs(i),Lsy_unabs(i) ! Lsy_unabs is already corrected for SSA
       enddo
       Eg = Eg*eV
       if (j < Ns) read(40,*)

       iang = int(m_ang_ff*abs(wo(j)-ang_min)/pi)+1
       do i=1,NSED_sy
          Lsy_tot_unabs(i) = Lsy_tot_unabs(i)+Lsy_unabs(i)
          Lsy_tot_noabs(i) = Lsy_tot_noabs(i)+Lsy_noabs(i)
          if (Ef_min_ff < Eg(i) .and. Eg(i) < Ef_max_ff) then
             iEf = int(log(Eg(i)/Ef_min_ff)/log(Ef_int_ff)) + 1
             Lsy_abs(i) = Lsy_unabs(i) * exp( -(a/d(j))**3*tau_ff(iEf,iang) )
             Lsy_tot_abs(i) = Lsy_tot_abs(i) + Lsy_abs(i)
          elseif  (Eg(i) <= Ef_min_ff) then
             Lsy_abs(i) = 0d0 ! totally absorbed
             Lsy_tot_abs(i) = 0d0 ! totally absorbed
          else
             Lsy_abs(i) = Lsy_unabs(i)
             Lsy_tot_abs(i) = Lsy_tot_unabs(i)
          endif

       	  !Jet contribution of the segments up to given lengths
       	  if (lenj(j)-lenj(1) <= 0.1d0*a) Lsy1(i) = Lsy_tot_abs(i)
       	  if (lenj(j)-lenj(1) <= a) Lsy2(i) = Lsy_tot_abs(i)

          write(41,*) Eg(i)/eV,Lsy_noabs(i),Lsy_unabs(i),Lsy_abs(i),j
       enddo
       write(41,*)
    enddo
    close(40)
    close(41)

    open(42,file='SEDsyn_abs_jet.dat')
    do i = 1,NSED_sy
       write(42,*) Eg(i)/eV,Lsy1(i),Lsy2(i),Lsy_tot_noabs(i),Lsy_tot_unabs(i),Lsy_tot_abs(i)   	  
    enddo
    close(42)

    deallocate(tau_ff,lenj,wo,d)
  END subroutine opac_ff_calc

  
  !============================================================================
  ! opac_ff_calc (counter-jet)
  !============================================================================

  subroutine opac_ff_calcCJ ()
    real*8 :: dummy,ang_int,ang_min
    real*8 :: Eg(NSED_sy),Lsy_noabs(NSED_sy),Lsy_unabs(NSED_sy),Lsy_abs(NSED_sy)
    real*8 :: Lsy_tot_abs(NSED_sy),Lsy_tot_noabs(NSED_sy),Lsy_tot_unabs(NSED_sy),Lsy1(NSED_sy),Lsy2(NSED_sy)
    real*8 :: Ef_int_ff
    integer :: iEf,iang
    allocate(tau_ff(mE_ff,m_ang_ff),lenj(Ns),wo(Ns),d(Ns)) 

    open (88,file='opac_ff.dat')
    do i=1,mE_ff
    	read(88,*) (tau_ff(i,j), j=1,m_ang_ff)
    enddo
    close(88)

    Ef_int_ff = (Ef_max_ff/Ef_min_ff)**(1.d0/float(mE_ff-1))
    ang_int = pi/float(m_ang_ff)
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
    ! CORRECT SYNCHROTRON SED
    !-------------------------------------------------------------------
    
    Lsy_tot_unabs = 0.d0
    Lsy_tot_abs = 0.d0
    Lsy_tot_noabs = 0.d0
    open(40,file='SEDsyn_cuts_ordCJ.dat')
    open(41,file='SEDsyn_abs_cutsCJ.dat')
    do j = 1,Ns

       do i = 1,NSED_sy
          read(40,*) Eg(i),Lsy_noabs(i),Lsy_unabs(i) ! Lsy_unabs is already corrected for SSA
       enddo
       Eg = Eg*eV
       if (j < Ns) read(40,*)

       iang = int(m_ang_ff*abs(wo(j)-ang_min)/pi)+1
       do i=1,NSED_sy
          Lsy_tot_unabs(i) = Lsy_tot_unabs(i)+Lsy_unabs(i)
          Lsy_tot_noabs(i) = Lsy_tot_noabs(i)+Lsy_noabs(i)
          if (Ef_min_ff < Eg(i) .and. Eg(i) < Ef_max_ff) then
             iEf = int(log(Eg(i)/Ef_min_ff)/log(Ef_int_ff)) + 1
             Lsy_abs(i) = Lsy_unabs(i) * exp( -(a/d(j))**3*tau_ff(iEf,iang) )
             Lsy_tot_abs(i) = Lsy_tot_abs(i) + Lsy_abs(i)
          elseif  (Eg(i) < Ef_min_ff) then
             Lsy_abs(i) = 0d0 ! totally absorbed
             Lsy_tot_abs(i) = 0d0 ! totally absorbed
          else
             Lsy_abs(i) = Lsy_unabs(i)
             Lsy_tot_abs(i) = Lsy_tot_unabs(i)
           endif

       	  !Jet contribution of the segments up to given lengths
       	  if (lenj(j)-lenj(1) <= 0.1d0*a) Lsy1(i) = Lsy_tot_abs(i)
       	  if (lenj(j)-lenj(1) <= a) Lsy2(i) = Lsy_tot_abs(i)
       	  
          write(41,*) Eg(i)/eV,Lsy_noabs(i),Lsy_unabs(i),Lsy_abs(i),j
       enddo
       write(41,*)
    enddo
    close(40)
    close(41)
    
    open(42,file='SEDsyn_abs_CJ.dat')
    do i = 1,NSED_sy
       write(42,*) Eg(i)/eV,Lsy1(i),Lsy2(i),Lsy_tot_noabs(i),Lsy_tot_unabs(i),Lsy_tot_abs(i)   	  
    enddo
    close(42)

    deallocate(tau_ff,lenj,wo,d)
  END subroutine opac_ff_calcCJ

  
  !============================================================================
  ! tau_ff_calc 
  !============================================================================

  !     The outputs of the program are:
  !     1) opac_ff.dat, which gives a table with nu(Hz) and opacity tau_ff(nu,ang)
  !     2) windemiff.dat, which gives a table with nu(Hz) and the emission S_nu in mJy

  subroutine tau_ff_calc(tau_table)
    integer, parameter :: dp = kind(1.d0)
    integer, parameter :: ir_max = 300
    real(dp), parameter:: T_w = 0.4d0*Tstar  ! assume constant T_w (valid far from the star)
    real(dp),intent(out) :: tau_table(mE_ff,m_ang_ff)
    real(dp) :: Mdot_eff = Mdot/sqrt(filling) ! effective mass-loss rate (account for clumping)
    real(dp) :: nu(mE_ff),nu_int,g_ff(mE_ff),S(mE_ff)
    real(dp) :: alpha_ff,cte_alpha_ff,cte_nu,dr,A_ff,v_w,ne,ni
    real(dp) :: tau,dang,ang,cang,sang,xi,yi,ri
    real(dp) :: f(ir_max)
    integer :: iang,ie,ir

    ! Frequencies at which tau_ff is calculated (calculate gaunt factor)
    open(14,file='nu.dat')
    call vector_log(Ef_min_ff/h,Ef_max_ff/h,nu_int,nu) ! Calculate opacities for 100 MHz? - 10 THz
    do i=1,mE_ff       
       g_ff(i) = 9.77d0*(1.d0+0.13d0*log10(T_w**1.5d0/Zq_w/nu(i))) ! Approx.
       write(14,*) nu(i)
    enddo
    close(14)     
    
    cte_alpha_ff = 3.7d8*Zq_w**2*T_w**(-0.5d0)
    
    dang = pi/m_ang_ff
    dr = Rstar/16.d0    ! step for integration: R_final ~ a + Rstar*ir_max/16*cos(ang)
    do ie=1,mE_ff       ! loop in emitted photon energy: Ef(ie) = h*nu(ie)
       ang = dang/2.d0  ! initial angle
       cte_nu = nu(ie)**(-3)*(1.d0-exp(-h*nu(ie)/kB/T_w))*g_ff(ie)

       do iang=1,m_ang_ff ! loop in interaction angle: ang
          cang = cos(ang)
          sang = sin(ang)
          xi = a        ! initial position (taken > Rstar for simplicity)
          yi = 0.d0
          tau = 0.d0
          do ir=1,ir_max  ! loop in emitted photon trayectory: xi,yi
             xi = dr*cang+xi
             yi = dr*sang+yi
             ri = sqrt(xi**2 + yi**2)
             if (ri > Rstar) v_w = 1d5 + vW_inf*(1d0 - Rstar/ri)**0.8d0 !vW = vW_0 + vW_inf*(1-R/r)^0.8
             A_ff = Mdot_eff/(4.d0*pi*v_w*mu_w*mp)
             ni = A_ff/ri**2/(gamma_w+1.d0)  ! mui = mu*(g_e+1)
             ne = gamma_w*ni                      
             alpha_ff = cte_alpha_ff*ne*ni*cte_nu !1/cm
             f(ir) = alpha_ff ! The photon energies are the same
             ! tau = tau + alpha_ff*dr
          enddo
          call integrate(ir_max,dr,f,tau)
          
		  tau_table(ie,iang)=tau

          ang=ang+dang

       enddo

    enddo


    !------------------------------------------------------------------------------
    !     
    ! The total radio flux from the wind is given by Eq. 8 from Wright 1975:

    open(19,file='windemiff.dat')
    S = 0.d0
    i=1
    do while (i < mE_ff .and. nu(i) < 1.d-1*kB*T_w/h )! formula valid for h*nu << k*T      
       S(i) = (23.2d0/(dist/kpc)**2) * ( ((Mdot_eff/Msun_yr)/&
            (vW_inf/1.d5)/mu_w)**2 *nu(i)*gamma_w*g_ff(i)*Zq_w**2 )**(2.d0/3.d0) ! Jy = 1d-26 W/m/Hz
       if (nu(i) > 1.d-3*kB*T_w/h) then
          S(i) = S(i)*exp(-h*nu(i)/kB/T_w)      ! Invented... 
       endif

       write(19,*) nu(i),S(i)*1.d3 ! Stored in mJy
       i=i+1
    enddo
    close(19)

    !
    !-------------------------------------------------------------------------------
    
    return

  END subroutine tau_ff_calc

  !===============================================================================

END MODULE opac_ff
