module controles

implicit none
public

!===================================================
! Qué programas voy a correr
logical, parameter :: no = .FALSE.
logical, parameter :: yes = .TRUE.

logical, parameter :: dist = yes
logical, parameter :: radsyn = no!yes
logical, parameter :: radsyn_data = yes

!===================================================

END module controles
