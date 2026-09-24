# ===============================
# Makefile for JetCode (Fortran + f2py)
# ===============================

DEBUG   = N
FC      = gfortran -fopenmp

# Compiler flags
FFLAGS  = -ffree-line-length-none -fbackslash
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
	$(FC) $(FFLAGS) $(OBJECTS) -o $@

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
# Python module via f2py (meson backend, numpy >= 1.26 / py3.12)
# --------------------------
F2PY_SOURCES = constants.f90 system_parameters.f90 globals.f90 controls.f90 \
	mathematics.f90 losses.f90 initialize.f90 trajectory.f90 evolution.f90 \
	edistribution.f90 synrad.f90 icrad.f90 sscrad.f90 main.f90

f2py:
	@echo "Building Python module JetCode..."
	f2py -c \
	  $(F2PY_SOURCES) \
	  -m JetCode \
	  --backend meson \
	  --f90flags="$(FFLAGS) -fopenmp" \
	  -lgomp

# Clean everything (build products only)
clean:
	$(RM) *.x *.o *.mod *.bin *.log *.so *~
	$(RM) -rf build

# Generate dependencies
depend:
	python3 gen_depend.py

# Include dependencies
-include Makefile.depend
