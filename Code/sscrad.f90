module sscrad
   use constants
   use globals
   use numerical_utils, only: vector_log
   implicit none
   public :: SEDssc

   contains

   !============================================================================
   ! SEDssc(START)
   !============================================================================
   ! This routine calculates the intrinsic SSC emission in the fluid frame
   ! The input is the '' (jet) or 'CJ' (counter-jet). 
   !============================================================================
   subroutine SEDssc(jet)  
      implicit none
      character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet 
      real(dp) :: Ee(nIn), Ne(nIn,Ns), Ne_cut(nIn), Eint_e, dEe(nIn)
      real(dp) :: EphFF(NSED_sy,Ns), nph(NSED_sy,Ns), nph_cut(NSED_sy)
      real(dp) :: LicFF(NSED_IC,Ns), Lic(NSED_IC,Ns), EgFF(NSED_IC), Eg(NSED_IC)
      real(dp) :: Eg_int, EgFF_int, dEgFF, q_gamma, P_g, Lic_tot
      real(dp) :: Lic_jet(NSED_IC), Lic2_jet(NSED_IC), Lic3_jet(NSED_IC)
      real(dp) :: cte_ic, dEph
      integer :: i, l, i2, j
      ! Read from a file: Ee(i), Ne(i,l)
      ! N_e(i) = number of e- with energy Ee(i) per unit energy
      open(70,file='Ne_cuts.dat')
      
      do i = 1, nIn
         do l = 1, Ns
            read(70,*) Ee(i), Ne(i,l), dummy
            Ne(i,l) = Ne(i,l) * deltaZ(l)
         enddo
      enddo
      close(70)
      Eint_e = Ee(2)/Ee(1)
      do i = 1, nIn
         dEe(i) = Ee(i)*(Eint_e-1.d0)
      enddo
      

      ! Read from a file: Esy(i,l), E*Lsy_unabs(i,l), E*Lsy_abs(i,l), l
      ! Lsy(i) = luminosity emitted with energy Eph(i) per unit energy
      open(45,file='n_ph_syn_cuts_FF_'//jet//'.dat')
      
      do i = 1, NSED_sy
         do l = 1, Ns
            read(45,*) EphFF(i,l), nph(i,l)
         enddo
      enddo

      close(45)

      !===========================================================
      !       Inverse Compton (isotropic)
      !===========================================================

      open(26,file='SEDssc_cuts'//jet//'.dat')

      cte_ic = 8.d0*pi*re**2*c

 
      ! Define the energy vector in the OF in which to calculate the SED
      call vector_log(1e-12*MINVAL(delta_obs), Ee(nIn)*MAXVAL(delta_obs), Eg_int, Eg)
   
      EgFF_int = Eg_int ! This is a ratio, the delta_obs factor simplifies

      Lic_tot = 0.d0
      LicFF = 0d0
      Lic = 0d0
      ! All calculations are performed in the FF
      do l = 2, Ns
         EgFF = Eg /delta_obs(l)
         Ne_cut = Ne(:,l)  ! N(Ee) for the segment l
         nph_cut = nph(:,l)! nph(Eph) for the segment l
         do i = 1, NSED_IC    ! Energy of the emitted SSC photons
            q_gamma = 0.d0
            do j = 1, nIn ! Ee
               P_g = 0.d0
               i2 = 1
               do while (EphFF(i2, l) < EgFF(i) .and. i2 <= NSED_sy)
                  dEph = EphFF(i2,l)*(EphFF(2,l)/EphFF(1,l)-1.d0)  
                  P_g = P_g + f(EgFF(i), Ee(j), EphFF(i2,l)) * nph_cut(i2)*dEph
                  i2 = i2 + 1
               enddo
               q_gamma = q_gamma + Ne_cut(j)*P_g*dEe(j)
            enddo
            ! Convert from the FF to the observer frame before saving
            LicFF(i,l) = cte_ic*EgFF(i)*q_gamma
            Lic(i,l) = LicFF(i,l) * boost(l)
            write (26,*) Eg(i)/eV, Eg(i)*Lic(i,l), l

            dEgFF = EgFF(i)*(Eg_int-1.d0)
            Lic_tot = Lic_tot + Lic(i,l)*dEgFF
         enddo
      enddo
      close(26)

      write(*,'(3X,A,ES10.3E2,1X,A)') 'Jet SSC luminosity (FF):',Lic_tot,'erg/s'

      !Now save a total SED for the jet
      open(27,file='SEDssc_jet'//jet//'.dat')
      Lic_jet = 0.d0
      Lic2_jet = 0.d0
      Lic3_jet = 0.d0
      do i = 1, NSED_IC
         do l = 1, Ns
            Lic_jet(i) = Lic_jet(i) + Lic(i,l)            
         enddo
         write(27,*) Eg(i)/eV, Eg(i)*Lic_jet(i), Eg(i)*Lic2_jet(i), Eg(i)*Lic3_jet(i)
      enddo
      close(27)

      write(*,*) "FINISHED SSC"
   END subroutine SEDssc
   !============================================================================
   ! SEDssc (END)
   !============================================================================  

   !============================================================================
   ! SEDssc_data(START)
   !============================================================================
   ! This routine calculates the intrinsic SSC emission in the fluid frame
   ! at the energies E_data [cgs] (vector), and returns the specific flux at each 
   ! of those energies [cgs] (also a vector). The number of points in E_data has
   ! to be smaller than m_data_max defined in globals.f90.
   ! The input is the '' (jet) or 'CJ' (counter-jet). 
   !============================================================================
   subroutine SEDssc_data(E_data, S_data, jet)  
      implicit none
      real(dp), intent(in) :: E_data(:)
      real(dp), intent(out) :: S_data(:)  ! size m_data_max
      character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet 
      real(dp) :: Ee(nIn), Ne(nIn,Ns), Ne_cut(nIn), Eint_e, dEe(nIn)
      real(dp) :: EphFF(NSED_sy,Ns), nph(NSED_sy,Ns), nph_cut(NSED_sy)
      real(dp) :: LicFF(NSED_IC,Ns), Lic(NSED_IC,Ns), EgFF(NSED_IC), Eg(NSED_IC)
      real(dp) :: Eg_int, EgFF_int, dEgFF, q_gamma, P_g, Lic_tot
      real(dp) :: Lic_jet(NSED_IC), Lic2_jet(NSED_IC), Lic3_jet(NSED_IC)
      real(dp) :: cte_ic, dEph
      integer :: i, l, i2, j
      ! Read from a file: Ee(i), Ne(i,l)
      ! N_e(i) = number of e- with energy Ee(i) per unit energy
      open(70,file='Ne_cuts.dat')
      
      do i = 1, nIn
         do l = 1, Ns
            read(70,*) Ee(i), Ne(i,l), dummy
            Ne(i,l) = Ne(i,l) * deltaZ(l)
         enddo
      enddo
      close(70)
      Eint_e = Ee(2)/Ee(1)
      do i = 1, nIn
         dEe(i) = Ee(i)*(Eint_e-1.d0)
      enddo
      

      ! Read from a file: Esy(i,l), E*Lsy_unabs(i,l), E*Lsy_abs(i,l), l
      ! Lsy(i) = luminosity emitted with energy Eph(i) per unit energy
      open(45,file='n_ph_syn_cuts_FF_'//jet//'.dat')
      do i = 1, NSED_sy 
         do l = 1, Ns
            read(45,*) EphFF(i,l), nph(i,l)
         enddo
      enddo

      close(45)

      !===========================================================
      !       Inverse Compton (isotropic)
      !===========================================================

      open(26,file='SEDssc_cuts'//jet//'.dat')

      cte_ic = 8.d0*pi*re**2*c

 
      ! Define the energy vector in the OF in which to calculate the SED
      call vector_log(1e-12*MINVAL(delta_obs), Ee(nIn)*MAXVAL(delta_obs), Eg_int, Eg)
   
      EgFF_int = Eg_int ! This is a ratio, the delta_obs factor simplifies

      Lic_tot = 0.d0
      LicFF = 0d0
      Lic = 0d0
      ! All calculations are performed in the FF
      do l = 2, Ns
         EgFF = Eg /delta_obs(l)
         Ne_cut = Ne(:,l)  ! N(Ee) for the segment l
         nph_cut = nph(:,l)! nph(Eph) for the segment l
         do i = 1, NSED_IC    ! Energy of the emitted SSC photons
            q_gamma = 0.d0
            do j = 1, nIn ! Ee
               P_g = 0.d0
               i2 = 1
               do while (EphFF(i2, l) < EgFF(i) .and. i2 <= NSED_sy)
                  dEph = EphFF(i2,l)*(EphFF(2,l)/EphFF(1,l)-1.d0)  
                  P_g = P_g + f(EgFF(i), Ee(j), EphFF(i2,l)) * nph_cut(i2)*dEph
                  i2 = i2 + 1
               enddo
               q_gamma = q_gamma + Ne_cut(j)*P_g*dEe(j)
            enddo
            ! Convert from the FF to the observer frame before saving
            LicFF(i,l) = cte_ic*EgFF(i)*q_gamma
            Lic(i,l) = LicFF(i,l) * boost(l)
            write (26,*) Eg(i)/eV, Eg(i)*Lic(i,l), l

            dEgFF = EgFF(i)*(Eg_int-1.d0)
            Lic_tot = Lic_tot + Lic(i,l)*dEgFF
         enddo
      enddo
      close(26)

      write(*,'(3X,A,ES10.3E2,1X,A)') 'Jet SSC luminosity (FF):',Lic_tot,'erg/s'

      !Now save a total SED for the jet
      open(27,file='SEDssc_jet'//jet//'.dat')
      Lic_jet = 0.d0
      Lic2_jet = 0.d0
      Lic3_jet = 0.d0
      do i = 1, NSED_IC
         do l = 1, Ns
            Lic_jet(i) = Lic_jet(i) + Lic(i,l)            
         enddo
         write(27,*) Eg(i)/eV, Eg(i)*Lic_jet(i), Eg(i)*Lic2_jet(i), Eg(i)*Lic3_jet(i)
      enddo
      close(27)

      write(*,*) "FINISHED SSC"
   END subroutine SEDssc_data
   !============================================================================
   ! SEDssc_data (END)
   !============================================================================

   !********************************************************************************
   real(dp) function f(Eg, Ee, Eph)     !f from Eq. 5.86 
   real(dp) :: Eg, Ee, Eph, z_, b_, t_  !Eg = gamma photon, Ee = electron, Eph = target
   z_ = Eg/Ee                      
   b_ = (4.d0*Eph*Ee)/(mec2**2)    
   t_ = z_/(1.d0-z_)
   if( ((Eph/Ee) < z_) .and. (z_ < (b_/(1.d0+b_))) ) then    ! Eq. 5.87
      f = ( 1.d0 + 0.5d0*z_*t_ + t_/b_ - 2.d0*(t_/b_)**2 - 0.5d0*(z_/b_)*t_**2 - 2.d0*(t_/b_)*log(b_/t_) )/(b_*Ee)   !1/erg
   else
      f = 0.d0
   endif
   end function f
   !********************************************************************************

!=======================================================================================

END module sscrad
