MODULE evolution
	use globals
	implicit none

contains

!*************************************************************************************************************
!	Evolution of different quantities from cell to cell (jet)
!*************************************************************************************************************
	subroutine Quantities ()
	real*8 :: B,nf,joFF
	allocate (x(Ns),y(Ns),z(Ns),rj(Ns),lenj(Ns),d(Ns),jw(Ns),wo(Ns),jo(Ns))
	allocate (delta_star(Ns),delta_obs(Ns))
	allocate (BFF(Ns),nfFF(Ns),vnorm(Ns),EphMFF(Ns),cwoFF(Ns),Kappa(Ns),sjoFF(Ns),VFF(Ns))
	
	open(10,file='Trajectory.dat')
	open(20,file='Deltas.dat')
	open(30,file='Quant_FF.dat')
	open(31,file='Quant_FF_norm.dat')
	
	do i = 1,Ns
		
		read(10,*) x(i),y(i),z(i),rj(i),lenj(i),d(i),jw(i),wo(i),jo(i)
	
		!Doppler factor between the jet and the stellar wind
		delta_star(i) = 1d0/(gj*(1d0-Bj*cos(jw(i))))
		
		!Doppler factor between the jet and the observer
		delta_obs(i) = 1d0/(gj*(1d0-Bj*cos(jo(i))))

		!Magnetic field
		B = B0*(rj(1)/rj(i))
		BFF(i) = B/gj
	
		!Stellar photon density and energy
		nf = nf0*(d(1)/d(i))**2
		nfFF(i) = nf/delta_star(i)
		EphMFF(i) = EphM/delta_star(i)
		
		!Dilution factor in the FF (Khangulyan 2014)
		Kappa(i) = (Rstar/(2d0*d(i)))**2 * delta_star(i)**2
		
		!Velocity perpendicular to the jet propagation
		if (i == 1) then
			vnorm(i) = vj*tan(thetaj)	!Conical expansion up to zrec (at least)
		else
			vnorm(i) = (rj(i)-rj(i-1))/dt
		endif
		
		!Jet observer angle in the FF
		joFF = acos((cos(jo(i))-Bj)/(1d0-Bj*cos(jo(i))))
		sjoFF(i) = sin(joFF)		!Sinus of joFF

		!Segment volume in the FF
		VFF(i) = pi*rj(i)**2*dtFF*vj
		
		!Cosine of of the angle between the star-jet and jet-observer vectors (C.11 in delaCita 2016)
		cwoFF(i) = 1d0 - (1d0 - cos(wo(i)))*delta_star(i)*delta_obs(i)
		
		write(20,*) delta_star(i),delta_obs(i),lenj(i)/a,delta_obs(i)**3/gj
		write(30,*) lenj(i),BFF(i),nfFF(i),rj(i),vnorm(i),EphMFF(i),cwoFF(i),Kappa(i),wo(i),sjoFF(i),VFF(i)
		write(31,*) lenj(i)/a,(BFF(i)/BFF(1))**2,(EphMFF(i)*nfFF(i))/(EphMFF(1)*nfFF(1)),	&
		&			(vnorm(i)/rj(i))/(vnorm(1)/rj(1)),cos(wo(i)),cwoFF(i)
	
	end do
	
	close(10)
	close(20)
	close(30)
	close(31)
	
	deallocate (x,y,z,rj,lenj,d,jw,wo,jo)
	deallocate (delta_star,delta_obs,BFF,nfFF,vnorm,EphMFF,cwoFF,Kappa,sjoFF,VFF)
	return
	end subroutine Quantities
	
!*************************************************************************************************************
!	Evolution of different quantities from cell to cell (counter-jet)
!*************************************************************************************************************
	subroutine QuantitiesCJ ()
	real*8 :: B,nf,joFF
	allocate (x(Ns),y(Ns),z(Ns),rj(Ns),lenj(Ns),d(Ns),jw(Ns),wo(Ns),jo(Ns))
	allocate (delta_star(Ns),delta_obs(Ns))
	allocate (BFF(Ns),nfFF(Ns),vnorm(Ns),EphMFF(Ns),cwoFF(Ns),Kappa(Ns),sjoFF(Ns),VFF(Ns))
	
	open(10,file='TrajectoryCJ.dat')
	open(20,file='DeltasCJ.dat')
	open(30,file='Quant_FFCJ.dat')
	open(31,file='Quant_FF_normCJ.dat')
	
	do i = 1,Ns
		
		read(10,*) x(i),y(i),z(i),rj(i),lenj(i),d(i),jw(i),wo(i),jo(i)
	
		!Doppler factor between the jet and the stellar wind
		delta_star(i) = 1d0/(gj*(1d0-Bj*cos(jw(i))))
		
		!Doppler factor between the jet and the observer
		delta_obs(i) = 1d0/(gj*(1d0-Bj*cos(jo(i))))

		!Magnetic field
		B = B0*(rj(1)/rj(i))
		BFF(i) = B/gj
	
		!Stellar photon density and energy
		nf = nf0*(d(1)/d(i))**2
		nfFF(i) = nf/delta_star(i)
		EphMFF(i) = EphM/delta_star(i)
		
		!Dilution factor in the FF (Khangulyan 2014)
		Kappa(i) = (Rstar/(2d0*d(i)))**2 * delta_star(i)**2
		
		!Velocity perpendicular to the jet propagation
		if (i == 1) then
			vnorm(i) = vj*tan(thetaj)	!Conical expansion up to zrec (at least)
		else
			vnorm(i) = (rj(i)-rj(i-1))/dt
		endif
		
		!Jet observer angle in the FF
		joFF = acos((cos(jo(i))-Bj)/(1d0-Bj*cos(jo(i))))
		sjoFF(i) = sin(joFF)		!Sinus of joFF

		!Segment volume in the FF
		VFF(i) = pi*rj(i)**2*dtFF*vj
		
		!Cosine of of the angle between the star-jet and jet-observer vectors (C.11 in delaCita 2016)
		cwoFF(i) = 1d0 - (1d0 - cos(wo(i)))*delta_star(i)*delta_obs(i)
		
		write(20,*) delta_star(i),delta_obs(i),lenj(i)/a,delta_obs(i)**3/gj
		write(30,*) lenj(i),BFF(i),nfFF(i),rj(i),vnorm(i),EphMFF(i),cwoFF(i),Kappa(i),wo(i),sjoFF(i),VFF(i)
		write(31,*) lenj(i)/a,(BFF(i)/BFF(1))**2,(EphMFF(i)*nfFF(i))/(EphMFF(1)*nfFF(1)),	&
		&			(vnorm(i)/rj(i))/(vnorm(1)/rj(1)),cos(wo(i)),cwoFF(i)
	
	enddo
	
	close(10)
	close(20)
	close(30)
	close(31)
	
	deallocate (x,y,z,rj,lenj,d,jw,wo,jo)
	deallocate (delta_star,delta_obs,BFF,nfFF,vnorm,EphMFF,cwoFF,Kappa,sjoFF,VFF)
	return
	end subroutine QuantitiesCJ
	
END MODULE evolution
