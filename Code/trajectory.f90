MODULE trajectory
	use globals
	implicit none
  
contains
  
!*************************************************************************************************************
!	Computation of the jet trajectory
!*************************************************************************************************************
	subroutine JetPath ()
    real*8 :: SO(3),r(Ns),phi(Ns),Pr,Pphi,Pz,Px,Py,theta,Ptot,jorb,vW,PdotW,rhoW,S,Sphi,	&
    &		  Fwr,Fphi,Fwz,vjr,vjphi,vjz,y0,z0
    integer :: nOcculted
	allocate (x(Ns),y(Ns),z(Ns),rj(Ns),lenj(Ns),d(Ns),jw(Ns),wo(Ns),jo(Ns),isOccult(Ns))
	
	write(*,'(X,A,X,I4)') 'Number of segments:',Ns
    
    !System-observer vector (unitary)
    SO = (/sin(inc*pi/180d0)*cos(azi*pi/180d0),sin(inc*pi/180)*sin(azi*pi/180d0),cos(inc*pi/180d0)/)
    
    !Initial conditions (cylindrical coordinates)
    r(1) = a
    phi(1) = orbPhase*pi/180d0
    z(1) = zrec
    Pr = 0d0		!Blob momentum projection in the r direction
    Pphi = 0		!Blob momentum projection in the azimuth direction
    Pz = Pdotj*dt	!Blob momentum projection in the z direction
    lenj(1) = zrec
    
    nOcculted = 0	!Number of segments occulted by the star
    
    do i = 1,Ns
    
    	!Cartesian coordinates
    	x(i) = r(i)*cos(phi(i))
    	y(i) = r(i)*sin(phi(i))
    	Px = Pr*cos(phi(i)) + Pphi*sin(phi(i))
    	Py = Pr*sin(phi(i)) - Pphi*cos(phi(i))
    	
    	Ptot = sqrt(Pr**2+Pphi**2+Pz**2)	!Total momentum of the jet
       	d(i) = sqrt(r(i)**2+z(i)**2)		!Distance to the star
       	theta = atan(z(i)/r(i))   			!Angle between the wind and the orbital plane
       	vjr = vj*Pr/Ptot					!Jet speed in the x direction
       	vjphi = vj*Pphi/Ptot				!Jet speed in the azimuth direction
       	vjz = vj*Pz/Ptot					!Jet speed in the z direction
	   	rj(i) = lenj(i)*tan(thetaj)			!Jet radius
	   	
       	jw(i) = acos((Px*x(i) + Py*y(i) + Pz*z(i))/(Ptot*d(i)))		!Jet-wind angle
       	wo(i) = acos((x(i)*SO(1) + y(i)*SO(2) + z(i)*SO(3))/d(i))	!Wind-observer angle
	   	jo(i) = acos((Px*SO(1) + Py*SO(2) + Pz*SO(3))/Ptot)			!Jet-observer angle
	   	jorb = acos(vjphi/vj)										!Jet-orbit angle	   	

    	vW = vW_inf*(1d0-Rstar/d(i))**0.8d0			!Wind speed at each point
    	PdotW = Mdot*vW								!Wind thrust
    	rhoW = Mdot/(4d0*pi*d(i)**2*vW)				!Wind density
       	S = 2d0*rj(i)*dt*vj*sin(jw(i))				!Blob surface normal to the star direction
       	Sphi = 2d0*rj(i)*dt*vj*sin(jorb)    		!Blob surface normal to the orbital direction
       	
       	Fwr = S*PdotW/(4d0*pi*d(i)**2)*cos(theta)	!Force made by the wind in the r direction
       	Fwz = S*PdotW/(4d0*pi*d(i)**2)*sin(theta)	!Force made by the wind in the z direction
       	Fphi = rhoW*Sphi*min(2d0*pi*r(i)/T,4d0*pi*(r(i)-a)/T)**2	!Coriolis + centrifugal forces in the phi
       																!direction (v_perp = min[w*r , 2*w*vW*t])
		!Increase in the momentum
       	Pr = Pr + Fwr*dt
       	Pphi = Pphi + Fphi*dt
       	Pz = Pz + Fwz*dt
       
       	!Jet trajectory
       	if (i /= Ns) then
          	r(i+1) = r(i) + vjr*dt
          	phi(i+1) = phi(i) - vjphi/r(i)*dt
          	z(i+1) = z(i) + vjz*dt
          	lenj(i+1) = lenj(i) + vj*dt
       	endif
       	
	   	!Coordinates at x=0 of a photon emitted at (x(i),y(i),z(i)) and directed towards the observer (valid
    	!only for SO(1) /= 0)
    	y0 = y(i)-SO(2)/SO(1)*x(i)
    	z0 = z(i)-SO(3)/SO(1)*x(i)
    	
    	!Check if the segment is occulted by the star (only valid when SO(1) /= 0 => azi /= +-90º)
    	isOccult(i) = .false.
    	if (y0**2+z0**2 <= Rstar**2 .and. rj(i) < 1.5d0*Rstar .and. x(i)*SO(1) < 0d0) then
    		isOccult(i) = .true.
    		nOcculted = nOcculted + 1
    	endif
       	
    enddo
    
    open(10,file='Trajectory.dat')
    do i = 1,Ns
       	write(10,*) x(i),y(i),z(i),rj(i),lenj(i),d(i),jw(i),wo(i),jo(i),isOccult(i)
    enddo
    close(10)
    
    write(*,'(X,A,X,F4.1,A)') 'Jet incliation:',180/pi*atan(Pr/Pz),'º'
    !write(*,'(X,A,X,F4.1,A)') 'Jet incliation (BR):',180/pi*Xj*pi**2/(2d0*pi*Xj+4d0*sqrt(Xj)+pi**2),'º'
    if (nOcculted > 0) write(*,'(X,A,X,F4.1,X,A)') '% of occulted jet =',dble(nOcculted)/dble(Ns)*100d0,'%'
    
	deallocate(x,y,z,rj,lenj,d,jw,wo,jo,isOccult)
    return
	end subroutine JetPath
	
	
!*************************************************************************************************************
!	Computation of the counter-jet trajectory
!*************************************************************************************************************
	subroutine JetPathCJ ()
    real*8 :: SO(3),r(Ns),phi(Ns),Pr,Pphi,Pz,Px,Py,theta,Ptot,jorb,vW,PdotW,rhoW,S,Sphi,	&
    &		  Fwr,Fphi,Fwz,vjr,vjphi,vjz,y0,z0
    integer :: nOcculted
	allocate (x(Ns),y(Ns),z(Ns),rj(Ns),lenj(Ns),d(Ns),jw(Ns),wo(Ns),jo(Ns),isOccult(Ns))
	
    !System-observer vector (unitary)
    SO = (/sin(inc*pi/180d0)*cos(azi*pi/180d0),sin(inc*pi/180)*sin(azi*pi/180d0),cos(inc*pi/180d0)/)
    
    !Initial conditions (cylindrical coordinates)
    r(1) = a
    phi(1) = orbPhase*pi/180d0
    z(1) = -zrec
    Pr = 0d0		!Blob momentum projection in the r direction
    Pphi = 0d0		!Blob momentum projection in the azimuth direction
    Pz = -Pdotj*dt	!Blob momentum projection in the z direction
    lenj(1) = zrec
    
    nOcculted = 0	!Number of segments occulted by the star
    
     do i = 1,Ns
    
    	!Cartesian coordinates
    	x(i) = r(i)*cos(phi(i))
    	y(i) = r(i)*sin(phi(i))
    	Px = Pr*cos(phi(i)) + Pphi*sin(phi(i))
    	Py = Pr*sin(phi(i)) - Pphi*cos(phi(i))
    	
    	Ptot = sqrt(Pr**2+Pphi**2+Pz**2)	!Total momentum of the jet
       	d(i) = sqrt(r(i)**2+z(i)**2)		!Distance to the star
       	theta = atan(z(i)/r(i))   			!Angle between the wind and the orbital plane
       	vjr = vj*Pr/Ptot					!Jet speed in the x direction
       	vjphi = vj*Pphi/Ptot				!Jet speed in the azimuth direction
       	vjz = vj*Pz/Ptot					!Jet speed in the z direction
	   	rj(i) = lenj(i)*tan(thetaj)			!Jet radius
	   	
       	jw(i) = acos((Px*x(i)+Py*y(i)+Pz*z(i))/(Ptot*d(i)))		!Jet-wind angle
       	wo(i) = acos((x(i)*SO(1)+y(i)*SO(2)+z(i)*SO(3))/d(i))	!Wind-observer angle
	   	jo(i) = acos((Px*SO(1)+Py*SO(2)+Pz*SO(3))/Ptot)			!Jet-observer angle
	   	jorb = acos(vjphi/vj)									!Jet-orbit angle	   	

    	vW = vW_inf*(1d0-Rstar/d(i))**0.8d0			!Wind speed at each point
    	PdotW = Mdot*vW								!Wind thrust
    	rhoW = Mdot/(4d0*pi*d(i)**2*vW)				!Wind density
       	S = 2d0*rj(i)*dt*vj*sin(jw(i))				!Blob surface normal to the star direction   
       	Sphi = 2d0*rj(i)*dt*vj*sin(jorb)    		!Blob surface normal to the orbital direction
       	
       	Fwr = S*PdotW/(4d0*pi*d(i)**2)*cos(theta)	!Force made by the wind in the r direction
       	Fwz = S*PdotW/(4d0*pi*d(i)**2)*sin(theta)	!Force made by the wind in the z direction
       	Fphi = rhoW*Sphi*min(2d0*pi*r(i)/T,4d0*pi*(r(i)-a)/T)**2	!Coriolis + centrifugal forces in the phi
       																!direction (v_perp = min[w*r , 2*w*vW*t])       
       	!Increase in the momentum
       	Pr = Pr + Fwr*dt
       	Pphi = Pphi + Fphi*dt
       	Pz = Pz + Fwz*dt
       
       	!Jet trajectory
       	if (i /= Ns) then
          	r(i+1) = r(i) + vjr*dt
          	phi(i+1) = phi(i) - vjphi/r(i)*dt
          	z(i+1) = z(i) + vjz*dt
          	lenj(i+1) = lenj(i) + vj*dt
       	endif
       	    	
    	!Coordinates at x=0 of a photon emitted at (x(i),y(i),z(i)) and directed towards the observer (valid
    	!only for SO(1) /= 0)
    	y0 = y(i)-SO(2)/SO(1)*x(i)
    	z0 = z(i)-SO(3)/SO(1)*x(i)
    	
    	!Check if the segment is occulted by the star (only valid when SO(1) /= 0 => azi /= +-90º)
    	isOccult(i) = .false.
    	if (y0**2+z0**2 <= Rstar**2 .and. rj(i) < 1.5d0*Rstar .and. x(i)*SO(1) < 0d0) then
    		isOccult(i) = .true.
    		nOcculted = nOcculted + 1
    	endif
       	
    enddo
    
    open(10,file='TrajectoryCJ.dat')
    do i = 1,Ns
       	write(10,*) x(i),y(i),z(i),rj(i),lenj(i),d(i),jw(i),wo(i),jo(i),isOccult(i)
    enddo
    close(10)
    
    if (nOcculted > 0) write(*,'(X,A,X,F4.1,X,A)') '% of occulted c-jet =',dble(nOcculted)/dble(Ns)*100d0,'%'
    
	deallocate(x,y,z,rj,lenj,d,jw,wo,jo,isOccult)
    return
	end subroutine JetPathCJ
	
!*************************************************************************************************************
!	Computation of a vertical jet trajectory
!*************************************************************************************************************
	subroutine JetPath_nocurv ()
    real*8 :: SO(3),Px,Py,Pz,Ptot,y0,z0
    integer :: nOcculted
	allocate (x(Ns),y(Ns),z(Ns),rj(Ns),lenj(Ns),d(Ns),jw(Ns),wo(Ns),jo(Ns),isOccult(Ns))
	
	write(*,'(X,A,X,I4)') 'Number of segments:',Ns
    
    !System-observer vector (unitary)
    SO = (/sin(inc*pi/180d0)*cos(azi*pi/180d0),sin(inc*pi/180)*sin(azi*pi/180d0),cos(inc*pi/180d0)/)
    
    !Initial conditions
    x(1) = a*cos(orbPhase*pi/180d0)
    y(1) = a*sin(orbPhase*pi/180d0)
    z(1) = zrec
    Px = 0d0		!Blob momentum projection in the x direction
    Py = 0d0		!Blob momentum projection in the y direction
    Pz = Pdotj*dt	!Blob momentum projection in the z direction
    lenj(1) = zrec
    
    nOcculted = 0	!Number of segments occulted by the star
    
    do i = 1,Ns
    	
    	Ptot = sqrt(Px**2+Py**2+Pz**2)			!Total momentum of the jet      
       	d(i) = sqrt(x(i)**2+y(i)**2+z(i)**2)	!Distance to the star
	   	rj(i) = lenj(i)*tan(thetaj)				!Jet radius	   	

       	jw(i) = acos((Px*x(i)+Py*y(i)+Pz*z(i))/(Ptot*d(i)))		!Jet-wind angle
       	wo(i) = acos((x(i)*SO(1)+y(i)*SO(2)+z(i)*SO(3))/d(i))	!Wind-observer angle
	   	jo(i) = acos((Px*SO(1)+Py*SO(2)+Pz*SO(3))/Ptot)			!Jet-observer angle
       
       	!Jet trajectory
       	if (i /= Ns) then
          	x(i+1) = x(i)
          	y(i+1) = y(i)
          	z(i+1) = z(i) + vj*dt
          	lenj(i+1) = lenj(i) + vj*dt
       	endif
    	
    	!Coordinates at x=0 of a photon emitted at (x(i),y(i),z(i)) and directed towards the observer (valid
    	!only for SO(1) /= 0)
    	y0 = y(i)-SO(2)/SO(1)*x(i)
    	z0 = z(i)-SO(3)/SO(1)*x(i)
    	
    	!Check if the segment is occulted by the star (only valid when SO(1) /= 0 => azi /= +-90º)
    	isOccult(i) = .false.
    	if (y0**2+z0**2 <= Rstar**2 .and. rj(i) < 1.5d0*Rstar .and. x(i)*SO(1) < 0d0) then
    		isOccult(i) = .true.
    		nOcculted = nOcculted + 1
    	endif
    		
    enddo
    
    open(10,file='Trajectory.dat')
    do i = 1,Ns
       	write(10,*) x(i),y(i),z(i),rj(i),lenj(i),d(i),jw(i),wo(i),jo(i),isOccult(i)
    enddo
    close(10)
    
    write(*,'(X,A,X,F4.1,A)') 'phi_inf =',180/pi*atan(Px/Pz),'º'	!Jet inclination at infinity
    if (nOcculted > 0) write(*,'(X,A,X,F4.1,X,A)') '% of occulted jet =',dble(nOcculted)/dble(Ns)*100d0,'%'
    
	deallocate(x,y,z,rj,lenj,d,jw,wo,jo,isOccult)
    return
	end subroutine JetPath_nocurv
	
	
!*************************************************************************************************************
!	Computation of a vertical counter-jet trajectory
!*************************************************************************************************************
	subroutine JetPathCJ_nocurv ()
    real*8 :: SO(3),Px,Py,Pz,Ptot,y0,z0
    integer :: nOcculted
	allocate (x(Ns),y(Ns),z(Ns),rj(Ns),lenj(Ns),d(Ns),jw(Ns),wo(Ns),jo(Ns),isOccult(Ns))
	
    !System-observer vector (unitary)
    SO = (/sin(inc*pi/180d0)*cos(azi*pi/180d0),sin(inc*pi/180)*sin(azi*pi/180d0),cos(inc*pi/180d0)/)
    
    !Initial conditions
    x(1) = a*cos(orbPhase*pi/180d0)
    y(1) = a*sin(orbPhase*pi/180d0)
    z(1) = -zrec
    Px = 0d0		!Blob momentum projection in the x direction
    Py = 0d0		!Blob momentum projection in the y direction
    Pz = -Pdotj*dt	!Blob momentum projection in the z direction
    lenj(1) = zrec
    
    nOcculted = 0	!Number of segments occulted by the star
    
    do i = 1,Ns
    	
    	Ptot = sqrt(Px**2+Py**2+Pz**2)			!Total momentum of the jet      
       	d(i) = sqrt(x(i)**2+y(i)**2+z(i)**2)	!Distance to the star
	   	rj(i) = lenj(i)*tan(thetaj)				!Jet radius	   	

       	jw(i) = acos((Px*x(i)+Py*y(i)+Pz*z(i))/(Ptot*d(i)))		!Jet-wind angle
       	wo(i) = acos((x(i)*SO(1)+y(i)*SO(2)+z(i)*SO(3))/d(i))	!Wind-observer angle
	   	jo(i) = acos((Px*SO(1)+Py*SO(2)+Pz*SO(3))/Ptot)			!Jet-observer angle
       
       	!Jet trajectory
       	if (i /= Ns) then
          	x(i+1) = x(i)
          	y(i+1) = y(i)
          	z(i+1) = z(i) - vj*dt
          	lenj(i+1) = lenj(i) + vj*dt
       	endif
       	
    	!Coordinates at x=0 of a photon emitted at (x(i),y(i),z(i)) and directed towards the observer (valid
    	!only for SO(1) /= 0)
    	y0 = y(i)-SO(2)/SO(1)*x(i)
    	z0 = z(i)-SO(3)/SO(1)*x(i)
    	
    	!Check if the segment is occulted by the star (only valid when SO(1) /= 0 => azi /= +-90º)
    	isOccult(i) = .false.
    	if (y0**2+z0**2 <= Rstar**2 .and. rj(i) < 1.5d0*Rstar .and. x(i)*SO(1) < 0d0) then
    		isOccult(i) = .true.
    		nOcculted = nOcculted + 1
    	endif
    		
    enddo
    
    open(10,file='TrajectoryCJ.dat')
    do i = 1,Ns
       	write(10,*) x(i),y(i),z(i),rj(i),lenj(i),d(i),jw(i),wo(i),jo(i),isOccult(i)
    enddo
    close(10)
    
    if (nOcculted > 0) write(*,'(X,A,X,F4.1,X,A)') '% of occulted c-jet =',dble(nOcculted)/dble(Ns)*100d0,'%'
    
	deallocate(x,y,z,rj,lenj,d,jw,wo,jo,isOccult)
    return
	end subroutine JetPathCJ_nocurv
  
END MODULE trajectory
