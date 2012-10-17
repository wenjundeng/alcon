# Copyright 2011, 2012 Wenjun Deng <wdeng@wdeng.info>
#
# This file is part of ALCON
#
# ALCON is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ALCON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ALCON.  If not, see <http://www.gnu.org/licenses/>.


# Makefile to compile and run ALCON
# use with GNU make or a compatible alternative

# adjust the following parameters according to your system environment before
#   compiling and running ALCON

# MPI Fortran 95 compiler
MPIF90 := mpif90

# compiling options
FFLAGS := -O3

# MPI executor
MPIEXEC := mpiexec

# number of MPI processes to run ALCON
NPE := 2


.PHONY : build run clean cleanoutput

build : alcon

# SLEPc variables
-include $(SLEPC_DIR)/conf/slepc_variables

alcon : alcon.o alcon_input.o alcon_eqdata.o alcon_solver.o alcon_output.o \
 cspline_petsc.o
	$(MPIF90) $(FFLAGS) -o $@ $^ $(SLEPC_LIB)

testcsp : testcsp.o cspline_petsc.o
	$(MPIF90) $(FFLAGS) -o $@ $^ $(SLEPC_LIB)

testcsp.o : testcsp.F90 cspline_petsc.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

alcon.o : alcon.F90 alcon_input.o alcon_eqdata.o alcon_solver.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

cspline_petsc.o : cspline_petsc.F90
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

alcon_input.o : alcon_input.F90
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

alcon_eqdata.o : alcon_eqdata.F90 alcon_input.o cspline_petsc.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

alcon_solver.o : alcon_solver.F90 alcon_input.o alcon_output.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

alcon_output.o : alcon_output.F90 alcon_input.o
	$(MPIF90) $(FFLAGS) -c -o $@ $< $(FCPPFLAGS)

run : alcon cleanoutput
	mkdir -p output
	$(MPIEXEC) -n $(NPE) ./$< 2>&1 | tee $<.log

runtestcsp : testcsp
	$(MPIEXEC) -n $(NPE) ./$< 2>&1 | tee $<.log

clean :
	rm -f alcon testcsp *.o *.mod *.log

cleanoutput :
	rm -f output/*

