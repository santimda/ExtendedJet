module sscrad
   use constants
   use globals
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
      real(dp) :: Eph(NSED_sy), Lsy(NSED_sy,Ns), nph(NSED_sy,Ns), dotNsy, nph_cut(NSED_sy)
      real(dp) :: LicFF(NSED_IC,Ns), Lic(NSED_IC,Ns), EgFF(NSED_IC), Eg(NSED_IC)
      real(dp) :: Eg_int, EgFF_int, dEgFF, q_gamma, P_g, Lic_tot
      real(dp) :: Lic_jet(NSED_IC), Lic2_jet(NSED_IC), Lic3_jet(NSED_IC)
      real(dp) :: boost(Ns), cte_ic, dEph, Eint_ph
      allocate ( delta_obs(Ns) )

      !------------------------------------------------------------------------
      !  Reading information here  

      ! Read Doppler factor
      open(20, file='deltas'//jet//'.dat')
      do l = 1, Ns
         read(20,*) delta_obs(l), boost(l)
      enddo
      close(20)

      ! Read from a file: Ee(i), Ne(i,l)
      ! N_e(i) = number of e- with energy Ee(i) per unit energy
      open(70,file='Ne_cuts.dat')
      do i = 1, nIn
         do l = 1, Ns
            read(70,*) Ee(i), Ne(i,l), dummy
         enddo
      enddo
      close(70)
      Eint_e = Ee(2)/Ee(1)
      do i = 1, nIn
         dEe(i) = Ee(i)*(Eint_e-1.d0)
      enddo

      ! Read from a file: Esy(i,l), E*Lsy_unabs(i,l), E*Lsy_abs(i,l), l
      ! Lsy(i) = luminosity emitted with energy Eph(i) per unit energy
      open(27, file='SEDsyn_cuts_ord'//jet//'.dat')
      do l = 1, Ns 
         do i = 1, NSED_sy
            read(27,*) Eph(i), dummy, Lsy(i,l), dummy
            Eph(i) = Eph(i)*eV ! Convert from eV to cgs
            ! We need L_sy_FF, so we unboost
            ! LA ENERGIA A LA QUE SE CALCULA ESTO TAMBIEN ES LA DEBOOSTEADA?
            Lsy(i,l) = Lsy(i,l)/(Eph(i)*boost(l)) ! [erg/s/erg]
            dotNsy = Lsy(i,l)/Eph(i)              ! [1/s/erg]
            nph(i,l) = dotNsy/VFF(l) * (rj(l)/c) !!! MAKE SURE IT IS CONSISTENT WITH U_sy
         enddo
      enddo
      close(27)
      Eint_ph = Eph(2)/Eph(1)

      !===========================================================
      !       Inverse Compton (isotropic)
      !===========================================================

      open(26,file='SEDssc_cuts'//jet//'.dat')

      cte_ic = 8.d0*pi*re**2*c

      ! Define the energy vector in the OF in which to calculate the SED
      call vector_log( 1.d3*Eph(1)*MINVAL(delta_obs), Ee(nIn)*MAXVAL(delta_obs), Eg_int, Eg)
      EgFF_int = Eg_int ! This is a ratio, the delta_obs factor simplifies

      Lic_tot = 0.d0
      ! All calculations are performed in the FF
      do l = 1, Ns
         EgFF = Eg/delta_obs(l)
         Ne_cut = Ne(:,l)  ! N(Ee) for the segment l
         nph_cut = nph(:,l)! nph(Eph) for the segment l
         do i = 1, NSED_IC    ! Energy of the emitted SSC photons
            q_gamma = 0.d0
            do j = 1, nIn ! Ee
               P_g = 0.d0
               do i2 = 1, NSED_sy
                  dEph = Eph(i2)*(Eint_ph-1.d0)  
                  P_g = P_g + cte_ic * f(EgFF(i), Ee(j), Eph(i2)) * nph_cut(i2)*dEph
               enddo
               q_gamma = q_gamma + Ne_cut(j)*P_g*dEe(j)
            enddo
            ! Convert from the FF to the observer frame before saving
            LicFF(i,l) = EgFF(i)*q_gamma 
            Lic(i,l) = LicFF(i,l) * boost(l)
            write (26,*) Eg(i)/eV, Eg(i)*Lic(i,l), l

            dEgFF = EgFF(i)*(Eg_int-1.d0)
            Lic_tot = Lic_tot + Lic(i,l)*dEgFF
         enddo
      enddo
      close(26)

      !write(*,'(3X,A,ES10.3E2,1X,A)') 'Jet unabsorbed synchrotron luminosity (FF):',LBolFF,'erg/s'
      write(*,'(3X,A,ES10.3E2,1X,A)') 'Jet SSC luminosity (FF):',Lic_tot,'erg/s'

      !Now save a total SED for the jet
      open(27,file='SEDssc_jet'//jet//'.dat')
      Lic_jet = 0.d0
      Lic2_jet = 0.d0
      Lic3_jet = 0.d0
      do i = 1, NSED_IC
         do l = 1, Ns
            Lic_jet(i) = Lic_jet(i) + Lic(i,l)            
            if (z(l) <= 1.d1*z0) Lic2_jet(i) = Lic_jet(i)
            if (z(l) <= 1.d2*z0) Lic3_jet(i) = Lic_jet(i)
         enddo
         write(27,*) Eg(i)/eV, Eg(i)*Lic_jet(i), Eg(i)*Lic2_jet(i), Eg(i)*Lic3_jet(i)
      enddo
      close(27)

   deallocate ( delta_obs )

   END subroutine SEDssc
   !============================================================================
   ! SEDssc (END)
   !============================================================================  

   !********************************************************************************
   real(dp) function f(Eg, Ee, Eph)     !f from Eq. 5.86 
   real(dp) :: Eg, Ee, Eph, z_, b_, t_  !Eg = gamma photon, Ee = electron, Eph = target
   z_ = Eg/Ee                      
   b_ = (4.d0*Eph*Ee)/(mec2**2)    
   t_ = z_/(1.-z_)
   if( ((Eph/Ee) <= z_) .and. (z_ <= (b_/(1.d0+b_))) ) then    ! Eq. 5.87
      f = ( 1.d0 + 0.5*z_*t_ + t_/b_ - 2.d0*(t_/b_)**2 - 0.5*(z_/b_)*t_**2 - 2.d0*(t_/b_)*log(b_/t_) )/(b_*Ee)   !1/erg
   else
      f = 0.d0
   endif
   end function f
   !********************************************************************************

!=======================================================================================

END module sscrad
