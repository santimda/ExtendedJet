DEBUG  = N

FC     = gfortran -fopenmp
FFLAGS = -ffree-line-length-none -fbackslash
LIBS   = -llapack -lblas

ifeq ($(DEBUG), Y)
  FFLAGS += -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
else
  FFLAGS += -Wno-tabs -fPIC -fmax-errors=1 -O2 -march=native -ffast-math -funroll-loops 
#  FFLAGS += -Wall -Wno-tabs -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops
#  FFLAGS += -Wall -Wextra -Wimplicit-interface -fPIC -Werror -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-none -fbackslash
endif

SOURCES = $(wildcard *.f90)
OBJECTS = $(subst .f90,.o,$(SOURCES)) 

.PHONY: clean test

run.x: $(OBJECTS)
	$(FC) $(OBJECTS) $(LIBS) -o run.x

.SUFFIXES:
.SUFFIXES: .f90 .f .o
.f90.o:
	$(FC) -c $(FFLAGS) $<

.f.o:
	$(FC) -c $(FFLAGS) $<

test:
	sh test_set.sh

clean:
	rm -f *.x *.o *.mod *.dat *.bin *.log *~
	rm -rf ./test

include Makefile.depend
