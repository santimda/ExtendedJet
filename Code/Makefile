DEBUG  = N

FC     = gfortran -fopenmp
FFLAGS = -ffree-line-length-none -fbackslash
LIBS   = -llapack -lblas
MODDIR = mod

ifeq ($(DEBUG), Y)
  FFLAGS += -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
else
  FFLAGS += -Wno-tabs -fPIC -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops
endif

SOURCES = $(wildcard *.f90)
OBJECTS = $(SOURCES:.f90=.o)

.PHONY: clean test depend

run.x: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) $(LIBS) -o $@

%.o: %.f90 | $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

%.o: %.f | $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -c $< -o $@

$(MODDIR):
	mkdir -p $(MODDIR)

run: run.x
	@echo "Running program..."
	@start_time=$$(date +%s); \
	./run.x; \
	end_time=$$(date +%s); \
	elapsed=$$((end_time - start_time)); \
	echo "Execution time: $${elapsed} s"

test:
	sh test_set.sh

depend:
	$(FC) $(FFLAGS) -cpp -MM $(SOURCES) > Makefile.depend

clean:
	$(RM) *.x *.o *.mod *.dat *.bin *.log *~
	$(RM) -r $(MODDIR) ./test

include Makefile.depend
