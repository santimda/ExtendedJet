module radssc
   use constants
   use global
   use controles
   use radsyn

   implicit none
   private
   public :: radssc_run, radssc_data_run

   contains

   !============================================================================
   ! radssc_run
   !============================================================================
   ! This routine calculates the intrinsic SSC emission in the fluid frame
   !-----------------------------------------------------------------------------
   subroutine radssc_run(E_sy, L_sy, E_ic, L_ic)  
      implicit none
      real(dp), intent(in) :: E_sy(mE_sy), L_sy(mE_sy)
      real(dp), intent(out) :: E_ic(mE_ph), L_ic(mE_ph)
      real(dp) :: E_ic_int, dE_ic, f, q_gamma, P_g, cte_ic
      real(dp) :: dotNsy(mE_sy), n_sy(mE_sy), dE_sy, Eint_sy
      integer :: i2

      interface
      function f(E_g, E_e, E_f) 
         real(kind(1.d0)), intent(in) :: E_g, E_e, E_f
      end function f
      end interface

      !------------------------------------------------------------------------
      !  Reading information here  

      if (.not. dist_calc ) then
         ! N_e(i) = number of e- with energy E_e(i) per unit energy
         open (14,file='dist_e.dat')
         do i = 1, mE_e
            read(14,*) E_e(i), N_e(i)
            E_e(i) = E_e(i)*eV    !   [erg]
         enddo
         close(14) 
      endif
      Eint_e = E_e(2)/E_e(1)

      !------------------------------------------------------------------------
      ! Obtain the synchrotron SED 

      do i = 1, mE_sy
         ! dotNsy(i) = photon rate emitted with energy E_sy(i) per unit energy
         dotNsy(i) = L_sy(i)/E_sy(i)
      enddo
      Eint_sy = E_sy(2)/E_sy(1)
      n_sy = dotNsy/(4.d0*pi*R**2*c) * 2.24d0  ! Same as naima

      !===========================================================
      !       Inverse Compton (isotropic)
      !===========================================================

      cte_ic = 8.d0*pi*re**2*c

      call vector_log(max(1.d3*E_sy(1),1.d1*eV), E_e(mE_e), E_ic_int, E_ic)

      L_ic_tot = 0.d0
      do i = 1, mE_ph    ! E_f(emitidos)
         q_gamma = 0.d0
         do j = 1, mE_e ! E_e
            dE_e = E_e(j)*(Eint_e-1.d0)
            P_g = 0.d0
            do i2 = 1, mE_sy
               dE_sy = E_sy(i2)*(Eint_sy-1.d0)  
               P_g = P_g + cte_ic * f(E_ic(i), E_e(j), E_sy(i2)) * n_sy(i2)*dE_sy
            enddo
            q_gamma = q_gamma + N_e(j)*P_g*dE_e
         enddo
         L_ic(i) = E_ic(i)*q_gamma 

         dE_ic = E_ic(i)*(E_ic_int-1.d0)
         L_ic_tot = L_ic_tot + L_ic(i)*dE_ic
      enddo

      if (write_files) then
         open(17,file='radssc.dat')
         do i = 1, mE_ph
            write(17,*) E_ic(i)/eV, L_ic(i)
         enddo
         close(17)
      endif


   END subroutine radssc_run


   !=============================================================================
   ! radssc_data_run
   !=============================================================================
   ! This routine calculates the intrinsic SSC emission in the fluid frame
   ! at the energies E_data, and returns the flux at each of those energies
   !-----------------------------------------------------------------------------
   subroutine radssc_data_run(B_new, E_data, S_data)
      implicit none
      real(dp), intent(in) :: B_new
      real(dp), intent(in) :: E_data(:)
      real(dp), intent(out) :: S_data(m_data_max)
      real(dp) :: f, q_gamma, P_g, a
      real(dp) :: cte_ic, dotNsy(mE_sy), n_sy(mE_sy), E_sy(mE_sy), dE_sy, Eint_sy
      integer :: i2, m_data

      interface
      function f(E_g, E_e, E_f) 
         real(kind(1.d0)), intent(in) :: E_g, E_e, E_f
      end function f
      end interface

      !------------------------------------------------------------------------
      !  Reading information here  

      if (.not. dist_calc ) then
         ! N_e(i) = number of e- with energy E_e(i) per unit energy
         open (14,file='dist_e.dat')
         do i = 1, mE_e
            read(14,*) E_e(i), N_e(i)
            E_e(i) = E_e(i)*eV    !   [erg]
         enddo
         close(14) 
      endif
      Eint_e = E_e(2)/E_e(1)

      !------------------------------------------------------------------------
      ! Obtain the synchrotron SED 

      call radsyn_run(B_new, E_sy, L_sy) 

      ! L_sy(i) = luminosity emitted with energy E_sy(i) per unit energy
      do i = 1, mE_sy
         ! dotNsy(i) = photon rate emitted with energy E_sy(i) per unit energy
         dotNsy(i) = L_sy(i)/E_sy(i)
      enddo
      close(14) 
      Eint_sy = E_sy(2)/E_sy(1)
      n_sy = dotNsy/(4.d0*pi*R**2*c) * 2.24d0  ! Same as naima

      !------------------------------------------------------------------------
      !  Radiative calculations start here  

      m_data = SIZE(E_data)
      if (m_data > m_data_max) then
         print*, 'Too many points to fit, make m_data_max larger in global.f90!'
         stop
      endif
      allocate( L_ic_data(m_data) )

      cte_ic = 8.d0*pi*re**2*c

      do i = 1, m_data    ! E_f(emitidos)
         q_gamma = 0.d0
         do j = 1, mE_e ! E_e
            dE_e = E_e(j)*(Eint_e-1.d0)
            P_g = 0.d0
            do i2 = 1, mE_sy
               dE_sy = E_sy(i2)*(Eint_sy-1.d0)  
               P_g = P_g + cte_ic * f(E_data(i), E_e(j), E_sy(i2)) * n_sy(i2)*dE_sy
            enddo
            q_gamma = q_gamma + N_e(j)*P_g*dE_e
         enddo
         L_ic_data(i) = E_data(i)*q_gamma 
         !write (*,*) E_ic(i)/eV, L_ic(i)
      enddo

      ! Define the output vector for the first m_data positions
      S_data = 0.d0
      a = 4.d0*pi*d**2
      do i = 1, m_data
         S_data(i) = L_ic_data(i)/a
      enddo

      deallocate( L_ic_data )

   END subroutine radssc_data_run


   !==============================

END module radssc

!====================================================================================
real(dp) function f(E_ic, Ee, Eph)     !f from Eq. 5.86 
use constants
real(dp) :: E_ic, Ee, Eph, z, b_, t_  !E_ic = gamma photon, Ee = electron, Eph = target
z = E_ic/Ee                      
b_ = (4.d0*Eph*Ee)/(mec2**2)    
t_ = z/(1.-z)
if( ((Eph/Ee) <= z) .and. (z <= (b_/(1.d0+b_))) ) then    ! Eq. 5.87
   f = ( 1.d0 + 0.5*z*t_ + t_/b_ - 2.d0*(t_/b_)**2 - 0.5*(z/b_)*t_**2 - 2.d0*(t_/b_)*log(b_/t_) )/(b_*Ee)   !1/erg
else
   f = 0.d0
endif

end function f

