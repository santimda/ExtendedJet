MODULE evolution
   use globals
   implicit none

   contains

   !*************************************************************************************
   !   Evolution of different quantities from cell to cell (jet)
   !*************************************************************************************
   subroutine Quantities(jet)

      character(LEN=*), intent(in) :: jet ! 'cj' or 'jet'
      real(dp) :: B, nf, joFF
      allocate( z(Ns), rj(Ns) )
      allocate( delta_obs(Ns) )
      allocate( BFF(Ns), vnorm(Ns), Kappa(Ns), sjoFF(Ns), VFF(Ns) )

      open(10,file='trajectory'//jet//'.dat')
      open(20,file='deltas'//jet//'.dat')
      open(30,file='quant_FF'//jet//'.dat')
      open(31,file='quant_FF_norm'//jet//'.dat')

      do i = 1, Ns

         read(10,*) z(i),rj(i)

         !Doppler factor between the jet and the observer
         delta_obs(i) = 1d0/(gj*(1d0-Bj*cos(jo)))

         !Magnetic field
         B = B0*(rj(1)/rj(i)) ! AGREGAR OTRO INDICE ACA!
         BFF(i) = B/gj

         ! LO QUE SIGUE HACERLO PARA EL CMB ?
         !Stellar photon density and energy
         !nf = nf0*(d(1)/d(i))**2
         !nfFF(i) = nf/delta_star(i)
         !EphMFF(i) = EphM/delta_star(i)

         !Dilution factor in the FF (Khangulyan 2014)
         Kappa(i) = delta_obs(i)**2 ! CHECK, FOR CMB

         !Velocity perpendicular to the jet propagation
         if (i == 1) then
            vnorm(i) = vj*tan(thetaj)   !Conical expansion up to zrec (at least)
         else
            vnorm(i) = (rj(i)-rj(i-1))/dt
         endif

         !Jet observer angle in the FF
         joFF = acos( (cos(jo)-Bj)/(1d0-Bj*cos(jo)) )
         sjoFF(i) = sin(joFF)        !Sinus of joFF

         !Segment volume in the FF
         VFF(i) = pi*rj(i)**2*dtFF*vj

         write(20,*) delta_obs(i), delta_obs(i)**3/gj
         write(30,*) z(i), BFF(i), rj(i), vnorm(i), Kappa(i), sjoFF(i), VFF(i)
         write(31,*) z(i)/R_g, (BFF(i)/BFF(1))**2, (vnorm(i)/rj(i))/(vnorm(1)/rj(1))

      end do

      close(10)
      close(20)
      close(30)
      close(31)

      deallocate( z, rj )
      deallocate( delta_obs, BFF, vnorm, Kappa, sjoFF, VFF )
      return
   end subroutine Quantities


END MODULE evolution
