MODULE fluxes
  use globals
  implicit none

contains

!*************************************************************************************************************
!	Computation of the fluxes at different bands
!*************************************************************************************************************
  subroutine fluxes_calc (Eph_sy,Lsy_jet,Lsy_CJ,Eph_IC,LIC_jet,LIC_CJ,Flux)
    real*8, intent(in) :: Eph_sy(NSED_sy),Lsy_jet(NSED_sy),Lsy_CJ(NSED_sy)
    real*8, intent(in) :: Eph_IC(NSED_IC),LIC_jet(NSED_IC),LIC_CJ(NSED_IC)
    real*8, intent(out) :: Flux(N_bands)
    real*8 :: Lsy(NSED_sy),LIC(NSED_IC),Emins(N_bands),Emaxs(N_bands),F_sy(N_bands),F_IC(N_bands)
    real*8 :: dilution,Ein_sy,Ein_IC,dEph_sy(NSED_sy),dEph_IC(NSED_IC)

	!Total luminositites
    Lsy = Lsy_jet + Lsy_CJ
    LIC = LIC_jet + LIC_CJ
    
    dilution = 4d0*pi*dist**2
       
    Ein_sy = (Emax_sy/Emin_sy)**(1d0/(NSED_sy-1))	!Energy increase factor
    do i = 1,NSED_sy
       dEph_sy(i) = Eph_sy(i)*(Ein_sy-1d0)			!Energy interval
    enddo
   
    Ein_IC = (Emax_IC/Emin_IC)**(1d0/(NSED_IC-1))	!Energy increase factor
    do i = 1,NSED_IC
       dEph_IC(i) = Eph_IC(i)*(Ein_IC-1d0)			!Energy interval
    enddo

    !Bands:  1-100 keV, 0.1-100 MeV, 0.1-10 GeV, 10-100 GeV, 0.1-30 TeV
    Emins = (/ 1d3*eV,1d5*eV,1d8*eV,1d10*eV,1d11*eV /)
    Emaxs = (/ 1d5*eV,1d8*eV,1d10*eV,1d11*eV,3d13*eV /)

    do j = 1,N_bands
       F_sy(j) = flux_int(NSED_sy,Eph_sy,dEph_sy,Lsy,Emins(j),Emaxs(j))    
       F_IC(j) = flux_int(NSED_IC,Eph_IC,dEph_IC,LIC,Emins(j),Emaxs(j))
    enddo
    Flux = (F_sy + F_IC)/dilution
    
    return
  end subroutine fluxes_calc
  
!*************************************************************************************************************
!	Integrate the fluxes at a given band
!*************************************************************************************************************  
  real*8 function flux_int (N,Eph,dEph,L,Eph_min,Eph_max)
    integer, intent(in) :: N
    real*8, intent(in) :: Eph(N),dEph(N),L(N),Eph_min,Eph_max
    
    flux_int = 0d0
    do i = 1,N
       if (Eph(i) >= Eph_min .and. Eph(i) < Eph_max) then
          flux_int = flux_int + L(i)*dEph(i)/Eph(i)		!L = E*L_E
       endif
    enddo
    
    return
  end function flux_int
  
END MODULE fluxes
