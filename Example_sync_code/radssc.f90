module radssc
   use constants
   use global

   implicit none
   private
   public :: radssc_run, radssc_data_run

   contains

   !============================================================================
   ! radssc_run
   !============================================================================
   ! This routine calculates the intrinsic SSC emission in the fluid frame
   !-----------------------------------------------------------------------------
   subroutine radssc_run(L_SSC_tot)  
      implicit none
      real(dp), intent(out) :: L_SSC_tot
      real(dp) :: Lic(mE_f), Eg(mE_f), Eg_int, dEg, f, q_gamma, P_g, Lic_tot
      real(dp) :: cte_ic, dotNsy(m_nu), nph(m_nu), Eph(m_nu), dEph, Eint_ph
      integer :: i2

      interface
      function f(E_g, E_e, E_f) 
         real(kind(1.d0)), intent(in) :: E_g, E_e, E_f
      end function f
      end interface

      allocate( N_e(mE_e), E_e(mE_e), Lsy(m_nu) )

      !------------------------------------------------------------------------
      !  Reading information here  

      ! N_e(i) = number of e- with energy E_e(i) per unit energy
      open (14,file='dist_e.dat')
      do i = 1, mE_e
         read(14,*) E_e(i), N_e(i)
         E_e(i) = E_e(i)*eV    !   [erg]
      enddo
      close(14) 
      Eint_e = E_e(2)/E_e(1)

      ! Lsy(i) = luminosity emitted with energy Eph(i) per unit energy
      open (14,file='radsyn.dat')
      do i = 1, m_nu
         read(14,*) Eph(i), Lsy(i)
         Eph(i) = Eph(i)*eV ! Convert from eV to cgs
         ! dotNsy(i) = photon rate emitted with energy Eph(i) per unit energy
         dotNsy(i) = Lsy(i)/Eph(i)
      enddo
      close(14) 
      Eint_ph = Eph(2)/Eph(1)
      nph = dotNsy/(4.d0*pi*R**2*c) * 2.24d0  ! Same as naima

      !===========================================================
      !       Inverse Compton (isotropic)
      !===========================================================

      cte_ic = 8.d0*pi*re**2*c

      call vector_log(max(1.d3*Eph(1),1.d1*eV), E_e(mE_e), Eg_int, Eg)

      open(17,file='radssc.dat')
      Lic_tot = 0.d0
      do i = 1, mE_f    ! E_f(emitidos)
         q_gamma = 0.d0
         do j = 1, mE_e ! E_e
            dE_e = E_e(j)*(Eint_e-1.d0)
            P_g = 0.d0
            do i2 = 1, m_nu
               dEph = Eph(i2)*(Eint_ph-1.d0)  
               P_g = P_g + cte_ic * f(Eg(i), E_e(j), Eph(i2)) * nph(i2)*dEph
            enddo
            q_gamma = q_gamma + N_e(j)*P_g*dE_e
         enddo
         Lic(i) = Eg(i)*q_gamma 
         write (17,*) Eg(i)/eV, Lic(i)

         dEg = Eg(i)*(Eg_int-1.d0)
         Lic_tot = Lic_tot + Lic(i)*dEg
      enddo
      L_SSC_tot = Lic_tot
      close(17)

      deallocate ( N_e, E_e, Lsy )

   END subroutine radssc_run


   !=============================================================================
   ! radssc_data_run
   !=============================================================================
   ! This routine calculates the intrinsic SSC emission in the fluid frame
   ! at the energies E_data, and returns the flux at each of those energies
   !-----------------------------------------------------------------------------
   subroutine radssc_data_run(E_data, S_data)
      implicit none
      real(dp), intent(in) :: E_data(:)
      real(dp), intent(out) :: S_data(m_data_max)
      real(dp), allocatable :: Eg(:), Lic(:)
      real(dp) :: f, q_gamma, P_g, a
      real(dp) :: cte_ic, dotNsy(m_nu), nph(m_nu), Eph(m_nu), dEph, Eint_ph
      integer :: i2, m_data

      interface
      function f(E_g, E_e, E_f) 
         real(kind(1.d0)), intent(in) :: E_g, E_e, E_f
      end function f
      end interface

      allocate( N_e(mE_e), E_e(mE_e), Lsy(m_nu) )
      
      !------------------------------------------------------------------------
      !  Reading information here  

      ! N_e(i) = number of e- with energy E_e(i) per unit energy
      open (14,file='dist_e.dat')
      do i = 1, mE_e
         read(14,*) E_e(i), N_e(i)
         E_e(i) = E_e(i)*eV    !   [erg]
      enddo
      close(14) 
      Eint_e = E_e(2)/E_e(1)

      ! Lsy(i) = luminosity emitted with energy Eph(i) per unit energy
      open (14,file='radsyn.dat')
      do i = 1, m_nu
         read(14,*) Eph(i), Lsy(i)
         Eph(i) = Eph(i)*eV ! Convert from eV to cgs
         ! dotNsy(i) = photon rate emitted with energy Eph(i) per unit energy
         dotNsy(i) = Lsy(i)/Eph(i)
      enddo
      close(14) 
      Eint_ph = Eph(2)/Eph(1)
      nph = dotNsy/(4.d0*pi*R**2*c) * 2.24d0  ! Same as naima

      !------------------------------------------------------------------------
      !  Radiative calculations start here  

      m_data = SIZE(E_data)
      if (m_data > m_data_max) then
         print*, 'Too many points to fit, make m_data_max larger in global.f90!'
         stop
      endif
      allocate( Eg(m_data), Lic(m_data) )
      Eg = E_data

      cte_ic = 8.d0*pi*re**2*c

      do i = 1, m_data    ! E_f(emitidos)
         q_gamma = 0.d0
         do j = 1, mE_e ! E_e
            dE_e = E_e(j)*(Eint_e-1.d0)
            P_g = 0.d0
            do i2 = 1, m_nu
               dEph = Eph(i2)*(Eint_ph-1.d0)  
               P_g = P_g + cte_ic * f(Eg(i), E_e(j), Eph(i2)) * nph(i2)*dEph
            enddo
            q_gamma = q_gamma + N_e(j)*P_g*dE_e
         enddo
         Lic(i) = Eg(i)*q_gamma 
         !write (*,*) Eg(i)/eV, Lic(i)
      enddo

      ! Define the output vector for the first m_data positions
      S_data = 0.d0
      a = 4.d0*pi*d**2
      do i = 1, m_data
         S_data(i) = Lic(i)/a
      enddo

      deallocate( Eg, Lic )
      deallocate( N_e, E_e, Lsy )

   END subroutine radssc_data_run


   !==============================

END module radssc

!====================================================================================
real(dp) function f(Eg, Ee, Eph)     !f from Eq. 5.86 
use constants
real(dp) :: Eg, Ee, Eph, z, b_, t_  !Eg = gamma photon, Ee = electron, Eph = target
z = Eg/Ee                      
b_ = (4.d0*Eph*Ee)/(mec2**2)    
t_ = z/(1.-z)
if( ((Eph/Ee) <= z) .and. (z <= (b_/(1.d0+b_))) ) then    ! Eq. 5.87
   f = ( 1.d0 + 0.5*z*t_ + t_/b_ - 2.d0*(t_/b_)**2 - 0.5*(z/b_)*t_**2 - 2.d0*(t_/b_)*log(b_/t_) )/(b_*Ee)   !1/erg
else
   f = 0.d0
endif

end function f

