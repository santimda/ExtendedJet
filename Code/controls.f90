MODULE controls
	implicit none
	
	logical, parameter :: no = .FALSE.
	logical, parameter :: yes = .TRUE.

!===================================================
! Programs to run

	logical, parameter :: curvedJet = yes
	
	logical, parameter :: NE_calc = no
	logical, parameter :: NEint_calc = no!yes
	
	logical, parameter :: radsyn_calc = no!yes
	logical, parameter :: radic_calc = no
	
	logical, parameter :: remake_table_ff = no!yes	!Recompute tau for ff absorption?
	logical, parameter :: remake_table_gg = no	!Recompute tau for gg absorption?
	
!===================================================

END MODULE controls
