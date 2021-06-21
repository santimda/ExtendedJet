MODULE evolution
   use constants
   use globals
   implicit none

   contains

   !============================================================================
   ! Quantities (START)
   !============================================================================
   subroutine Quantities(jet)
   !============================================================================
   ! This subroutine calculates the evolution of different quantities 
   ! from cell to cell. 
   ! The input is the '' (jet) or 'CJ' (counter-jet). 
   !============================================================================
      character(LEN=*), intent(in) :: jet ! 'CJ' for counterjet or '' for jet
      real(dp) :: dir, cjo, boost, bj_, gj_
      allocate( joFF(Ns) )
      allocate( delta_obs(Ns) )

      open(20,file='deltas'//jet//'.dat')
      open(30,file='quant_FF'//jet//'.dat')
      open(31,file='quant_FF_norm'//jet//'.dat')

      ! The jet goes upwards, the CJ goes downwards (negative v_z)
      if (jet == '') then
         dir = 1.d0
      elseif (jet == 'CJ') then
         dir = -1.d0
      else
         write(*,*) 'Invalid jet type'
         stop
      endif

      cjo = cos(jo)
      do i = 1, Ns

         if (abs(z(i)) < zmax_inj) Ninj = i

         gj_ = lorentzFactor(z(i))
         bj_ = sqrt(1d0-1d0/gj_**2)
         deltaZ(i) = abs(z(i)) * (sqrt(dlogZ)-1d0/sqrt(dlogZ)) ! Added abs value for CJ
         deltatFF(i) = deltaZ(i) / (c*bj_*gj_)
         !print*, (deltaZ(i)/gj)/rj(i), z(i)/rj(i)/gj
         !Doppler factor between the jet and the observer
         delta_obs(i) = 1d0/(gj_*(1d0-dir*bj_*cjo))

         !Magnetic field
         BFF(i) = BFF_fun(z(i))

         !Stellar photon density and energy 
         nfFF(i) = nf0/delta_obs(i)  !!! CHECK!!!!!
         EphMFF(i) = EphM/delta_obs(i)

         !Dilution factor in the FF (Khangulyan 2014)
         Kappa(i) = delta_obs(i)**2 ! CHECK, FOR CMB

         !Velocity perpendicular to the jet propagation
         vnorm(i) = v_lat(z(i))

         !Jet observer angle in the FF
         joFF(i) = acos( (cjo-bj_)/(1d0-bj_*cjo) )
         sjoFF(i) = sin(joFF(i))        !Sinus of joFF

         !Segment volume in the FF
         VFF(i) = pi*rj(i)**2*(deltaZ(i)/gj_)

         !boost = delta_obs(i)**2/gj  ! CHECK THIS!
         boost = delta_obs(i)**2  ! CHECK THIS! (delta z defined in LF)

         write(20,*) delta_obs(i), boost
         write(30,*) z(i), BFF(i), rj(i), vnorm(i), Kappa(i), sjoFF(i), VFF(i), EphMFF(i), nfFF(i), deltaZ(i), deltatFF(i)
         write(31,*) z(i)/R_g, (BFF(i)/BFF(1))**2, (vnorm(i)/rj(i))/(vnorm(1)/rj(1)), VFF(i)/VFF(1)

      end do

      close(20)
      close(30)
      close(31)

      deallocate( delta_obs )
      deallocate( joFF )

   end subroutine Quantities
   !============================================================================
   ! Quantities (END)
   !============================================================================

END MODULE evolution
