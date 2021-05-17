PROGRAM main

  use global
  use controles
  use dist, only: diste_run
  use radsyn, only: radsyn_run

  implicit none

  !================================================================


  write(*,*) 'ELECTRONES:'
        
  if (dist.eqv.yes) then
     write(*,*) 'Calculo la distribución de electrones NT'
     call diste_run 
     write(*,*) '' 
  endif

  if (radsyn.eqv.yes) then
     write(*,*) 'Calculo la emisión sincrotrón de los electrones'
     call radsyn_run 
     write(*,*) '' 
  endif

  !================================================================

  write(*,*) 'Terminado todo'
  
END PROGRAM main

   
    
