MODULE losses
	use globals
	implicit none

contains

!*************************************************************************************************************
!	Total losses
!*************************************************************************************************************
	real*8 function TLosses (BFF,rj,vnorm,EphMFF,Kappa,EFF)
	real*8, intent(in) :: BFF,rj,vnorm,EphMFF,Kappa,EFF

	TLosses = SyncLoss(BFF,EFF) + ICLoss(EphMFF,Kappa,EFF) + AdLoss(rj,vnorm,EFF)

	return
	end function TLosses

!*************************************************************************************************************
!	Synchrotron losses
!*************************************************************************************************************
	real*8 function SyncLoss (BFF,EFF)
	real*8, intent(in) :: BFF,EFF
	real*8 :: Lgamma

	!Lorentz factor of the particles in the fluid frame
	Lgamma = EFF/(me*c**2)

	!Synchrotron losses in erg/s
	SyncLoss = - 1.05804d-15 * BFF**2 * Lgamma**2	!4/9*c*re^2 * B^2 * gamma^2

	return
	end function SyncLoss
	
!*************************************************************************************************************
!	Inverse Compton losses (Eq. (36) in Khangulyan 2014)					
!*************************************************************************************************************
	real*8 function ICLoss (EphMFF,Kappa,EFF)
	real*8, intent(in) :: EphMFF,Kappa,EFF
	real*8 :: TFF,t_iso,Eiso_ct
	
	TFF = EphMFF/(2.7d0*mec2)	!Dimensionless stellar temperature in the FF (EphMFF = 2.7*kB*TFF)
	
	t_iso = 4d0*TFF*(EFF/mec2)
	Eiso_ct = cte_IC_K*Kappa*TFF**2
	
	!IC losses (in erg/s due to the mec2 multiplying)
	ICLoss = - Eiso_ct*F_iso(t_iso)*mec2

	return
	end function ICLoss

!*************************************************************************************************************
!	Adiabatic losses
!*************************************************************************************************************
	real*8 function AdLoss (rj,vnorm,EFF)
	real*8, intent(in) :: rj,vnorm,EFF

	!Adiabatic losses
	AdLoss = - 2d0/3d0 * vnorm/rj * EFF * gj	!E/t = E'/t' = ct

	end function AdLoss
	
!*************************************************************************************************************
!	Write the different cooling times for a given segment
!*************************************************************************************************************
	subroutine Print_tloss(segN)
	integer, intent(in) :: segN
	real*8 :: EFF,Eint,tSyn,tIC,tAd,tTot,tacc,tdiff
	real*8 :: dummy
	allocate (BFF(Ns),rj(Ns),vnorm(Ns),EphMFF(Ns),Kappa(Ns))
	
	!Quantities at each segment
	open(10,file='Quant_FF.dat')
   	do i = 1,Ns
		read(10,*) dummy,BFF(i),dummy,rj(i),vnorm(i),EphMFF(i),dummy,Kappa(i)
	enddo
	close(10)
	
	!Initialization
	EFF = E1
	Eint = (E2/E1)**(1d0/(nNE-1))
	
	open(20,file='Losses.dat')	
	do i = 1,nNE
	
		!Cooling time for each process for segment segN
		tSyn = EFF/abs(SyncLoss(BFF(segN),EFF))
		tIC = EFF/abs(ICLoss(EphMFF(segN),Kappa(segN),EFF))
		tAd = EFF/abs(AdLoss(rj(segN),vnorm(segN),EFF))
		tTot = EFF/abs(TLosses(BFF(segN),rj(segN),vnorm(segN),EphMFF(segN),Kappa(segN),EFF))
		tacc = EFF/(facc*BFF(segN)*qe*c)					!Acceleration timescale
		tdiff = 3d0*qe*BFF(segN)*rj(segN)**2/(2d0*c*EFF)	!Diffusion timescale
		
		write(20,*) EFF/eV,tSyn,tIC,tAd,tTot,dtFF,tacc,tdiff
					
		EFF = EFF*Eint	!Energy increase
			
	enddo
	close(20)

	deallocate (BFF,rj,vnorm,EphMFF,Kappa)
	return
	end subroutine 

END MODULE losses
