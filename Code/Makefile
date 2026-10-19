# ===============================
# Makefile for JetCode (Fortran + f2py)
# ===============================

DEBUG    ?= N
PROFILE  ?= N
FC        = gfortran -fopenmp
BUILD_DIR = build
CONFIG   ?= ../configs/test.nml

FFLAGS = -ffree-line-length-none -fbackslash
ifeq ($(DEBUG), Y)
  FFLAGS += -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
else
  FFLAGS += -Wno-tabs -fPIC -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops
  ifeq ($(PROFILE), Y)
    FFLAGS += -g -fno-omit-frame-pointer
  endif
endif

# Production sources, ordered from low-level modules to the program entry point.
SOURCES = constants.f90 \
          controls.f90 \
          system_parameters.f90 \
          numerical_utils.f90 \
          globals.f90 \
          losses.f90 \
          initialize.f90 \
          trajectory.f90 \
          evolution.f90 \
          edistribution.f90 \
          synrad.f90 \
          icrad.f90 \
          sscrad.f90 \
          main.f90
OBJECTS = $(SOURCES:%.f90=$(BUILD_DIR)/%.o)

F2PY_SOURCES = constants.f90 controls.f90 system_parameters.f90 \
               numerical_utils.f90 globals.f90 losses.f90 initialize.f90 \
               trajectory.f90 evolution.f90 edistribution.f90 synrad.f90 \
               icrad.f90 sscrad.f90 main.f90

.PHONY: all clean run f2py depend test test-unit test-fortran test-f2py validate

all: run.x

run.x: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

.SUFFIXES:
$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

run: run.x
	@echo "Running program with $(CONFIG)..."
	@start_time=$$(date +%s); \
	./run.x $(CONFIG); \
	status=$$?; \
	end_time=$$(date +%s); \
	elapsed=$$((end_time - start_time)); \
	echo "Execution time: $${elapsed} s"; \
	exit $$status

test: test-unit test-fortran test-f2py

test-unit: $(BUILD_DIR)/test_numerical_utilities
	./$(BUILD_DIR)/test_numerical_utilities

test-fortran: run.x
	python3 ../tests/fortran_regression.py

# Physics validation: N(E) solver vs analytic/characteristics solutions (needs the f2py module: make f2py)
validate: run.x
	python3 ../validation/ne_validation.py

# f2py syn_model wrapper vs run.x on the same namelist (needs the f2py module: make f2py)
test-f2py: run.x
	python3 ../tests/f2py_check.py

$(BUILD_DIR)/test_numerical_utilities: ../tests/test_numerical_utilities.f90 \
                                      $(BUILD_DIR)/constants.o \
                                      $(BUILD_DIR)/numerical_utils.o | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) $^ -o $@

f2py:
	@echo "Building Python module JetCode..."
	f2py -c \
	  $(F2PY_SOURCES) \
	  -m JetCode \
	  --backend meson \
	  --f90flags="$(FFLAGS) -fopenmp" \
	  -lgomp

clean:
	$(RM) *.x *.o *.mod *.bin *.log *.so *~
	$(RM) -rf $(BUILD_DIR)

depend:
	python3 gen_depend.py

-include Makefile.depend
