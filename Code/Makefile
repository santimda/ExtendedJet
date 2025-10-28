# ===============================
# Makefile for JetCode (Fortran + f2py)
# ===============================

DEBUG   = N
FC      = gfortran -fopenmp
LIBS    = -llapack -lblas

# Compiler flags
FFLAGS  = -ffree-line-length-none -fbackslash $(INCDIR)
ifeq ($(DEBUG), Y)
  FFLAGS += -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
else
  FFLAGS += -Wno-tabs -fPIC -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops
endif

# Sources
SOURCES = $(wildcard *.f90)
OBJECTS = $(SOURCES:.f90=.o)

.PHONY: all clean run f2py depend

# Default target: compile executable
all: run.x

# Compile executable
run.x: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) $(LIBS) -o $@

# Compile each object file and generate .mod files in $(MODDIR)
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Run the executable with timing
run: run.x
	@echo "Running program..."
	@start_time=$$(date +%s); \
	./run.x; \
	end_time=$$(date +%s); \
	elapsed=$$((end_time - start_time)); \
	echo "Execution time: $${elapsed} s"

# --------------------------
# Python module via f2py
# --------------------------
f2py: clean
	@echo "Building Python module JetCode..."
	f2py -c \
	  controls.f90 \
	  constants.f90 \
	  system_parameters.f90 \
	  globals.f90 \
	  initialize.f90 \
	  trajectory.f90 \
	  evolution.f90 \
	  losses.f90 \
	  mathematics.f90 \
	  edistribution.f90 \
	  synrad.f90 \
	  icrad.f90 \
	  sscrad.f90 \
	  main.f90 \
	  -m JetCode \
	  --f90flags="$(FFLAGS)"

# f2py -c controls.f90 constants.f90 system_parameters.f90 globals.f90 initialize.f90 trajectory.f90 evolution.f90 losses.f90 mathematics.f90 edistribution.f90 synrad.f90 icrad.f90 sscrad.f90 main.f90 -m JetCode

# Clean everything
clean:
	$(RM) *.x *.o *.mod *.dat *.bin *.log *~
	$(RM) -rf $(MODDIR)

# Generate dependencies
depend:
	$(FC) $(FFLAGS) -cpp -MM $(SOURCES) > Makefile.depend

# Include dependencies
-include Makefile.depend

