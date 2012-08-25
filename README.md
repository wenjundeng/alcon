ALCON
=====
version 2012-08-24 23:04:13-04:00
---------------------------------

ALCON is a code for solving ideal MHD Alfven continua in tokamaks,
i.e., Eq. (10) in
[Physics of Fluids 29, 3695 (1986)](http://dx.doi.org/10.1063/1.865801),
using a poloidal-spectral method described in Appendix A in
[Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117).
When referencing ALCON in a publication, please cite:

	@ARTICLE{Deng2012b,
		author = {W. Deng and Z. Lin and I. Holod and Z. Wang and Y. Xiao and H. Zhang},
		title = {Linear properties of reversed shear {A}lfv\'en eigenmodes in the {DIII-D} tokamak},
		journal = {Nuclear Fusion},
		year = {2012},
		volume = {52},
		pages = {043006},
		number = {4},
		doi = {10.1088/0029-5515/52/4/043006}
	}


Copying
-------

ALCON is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ALCON is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ALCON.  If not, see <http://www.gnu.org/licenses/>.


System requirements
-------------------

ALCON requires MPI-2, [PETSc](http://www.mcs.anl.gov/petsc/) and
[SLEPc](http://www.grycap.upv.es/slepc/) with complex scalars.  ALCON has been
tested with [MPICH2](http://www.mcs.anl.gov/research/projects/mpich2/) 1.4,
[OpenMPI](http://www.open-mpi.org/) 1.5, PETSc and SLEPc 3.2.  Although with a
few tweaks ALCON can probably work with PETSc and SLEPc 3.0 and 3.1, it is
recommended to use PETSc and SLEPc 3.2+.

ALCON also requires a Fortran 90 compiler for compilation.  ALCON has been
tested with [GNU Fortran](http://gcc.gnu.org/fortran/) 4.6.  [GNU core
utilities](http://www.gnu.org/software/coreutils/) is required to use GNU make
for automatic compilation.  Typical installation of any major GNU/Linux
distribution has GNU core utilities pre-installed.

Note that BSD-like systems, e.g., FreeBSD and Apple Mac OS X, have utilities
similar to GNU core utilities, but they differ from the GNU version.  You need
to install the GNU version of them if you want to compile ALCON using GNU make
on BSD-like systems.  On Apple Mac OS X, GNU core utilities can be installed
through [Fink](http://www.finkproject.org/) or
[MacPorts](http://www.macports.org/).  You can also compile ALCON manually
without using GNU make.

Unless you keep excessive (100+) poloidal m-harmonics and want excessive high
radial resolution (solving 5000+ radial grid points), ALCON should be able to
finish a run within half an hour on any desktop that has two or more physical
cores.  In most practical cases, ALCON can finish a run within a few minutes on
a desktop within a few years old.


Usage
-----

### Know about the files

+ `README.md` -- This instruction.
+ `Makefile` -- For using GNU make to compile and run ALCON.
+ `alcon.F90` -- The main program file.
+ `alcon_input.F90` -- Module for managing input parameters. Change input
parameters in this file.
+ `alcon_eqdata.F90` -- Module for managing equilibrium data.
+ `alcon_solver.F90` -- Module for solving Alfven continua.
+ `alcon_output.F90` -- Module for managing output data and writing to output
files.
+ `cspline_petsc.F90` -- Module for cubic spline interpolation using PETSc.
+ `testcsp.F90` -- A small program to test `cspline_petsc.F90`, not an
essential part of ALCON, but if you are interested in playing with the cubic
spline module, help yourself.
+ `alcon.dat` -- A sample equilibrium data input file.  The equilibrium is
taken from [Physics of Plasmas 17, 112504 (2010)](http://wdeng.info/?p=37).
This file is for testing when you use ALCON for the first time.  For practical
use, you want to replace this file by your own `alcon.dat` file generated from
your equilibrium.  See `acdgen_sample.F90` for a sample subroutine to generate
`alcon.dat`.
+ `profile_sa.dat` -- A sample equilibrium profiles for solving a simple
analytic equilibrium.  The profiles are taken from
[Physics of Plasmas 17, 112504 (2010)](http://wdeng.info/?p=37).  This file is
for testing when you use ALCON for the first time.  For practical use, you want
to replace this file by your own `profile_sa.dat` file.  See next section for
preparing your own `profile_sa.dat`.
+ `acdgen_sample.F90` -- A sample subroutine to show how to generate
`alcon.dat`.


### Input preparation

ALCON reads input parameters from `alcon_input.F90`, so set your input
parameters in this file.  Detailed description of each parameter is included in
the comments in `alcon_input.F90`.

Depending on the eqtype parameter set in `alcon_input.F90`, ALCON reads
equilibrium data from additional input files.  For currently implemented
options for eqtype, ALCON may read from `alcon.dat` or `profile_sa.dat` for
equilibrium data.  ALCON comes with a sample of each of these two files.

To prepare your own `alcon.dat`, see `acdgen_sample.F90` for a sample
subroutine to generate `alcon.dat`.

To prepare your own `profile_sa.dat`, follow this format: the first line puts
one real number, which is the inverse aspect ratio a/R0; the second line puts
two integers `nrad` and `nprofile`, `nrad` being the number of radial grid
points for the profiles, `nprofile` being the number of profiles (for current
version always put `nprofile` to be 4); starting from the third line is a
`nrad * nprofile` matrix giving the profile data.  The first column is the radial
coordinate r/a.  The second column is the q-profile.  The third column is the
pressure profile normalized as: (4 pi gamma P / B^2), where gamma is the
specific heat ratio, P is the plasma pressure, B is the on-axis magnetic field,
and P and B are in CGS units.  The fourth column is the mass density normalized
by proton mass times on-axis electron density (Eq. (A.35) in
[Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117)).


### Compilation

It is recommended to use GNU make (or a compatible alternative) to compile the
code with the included Makefile.  Before compiling, open the Makefile and go to
the place after the license part, you will see a few parameters, such as
compiler name, compiling options, MPI executor, etc., defined there.  Adjust
them according to your system environment.  After saving your adjustments, you
can compile ALCON by executing in the terminal:

	make

To remove the compiled binary and all other temporary files during compilation,
execute:

	make clean


### Running

After successfully compiling ALCON, and making sure your MPI environment is
correctly configured, you can then run ALCON by executing:

	make run


### Output data and plotting

All output data files are put in the `output` sub-directory.  If you choose to
generate a plotting script (`outformat` option in `alcon_input.F90`), it is
also put in the `output` sub-directory.  Then you can use the corresponding
plotting software to produce an Alfven continuum figure.

**When you execute `make run`, the content of the output sub-directory will be
erased before ALCON runs.  If you think the output data files are useful, move
them immediately to a safe place, or simply rename the output sub-directory.**

You can also manually erase the output sub-directory by executing:

	make cleanoutput


Contact
-------

Get latest updates on the project website:
<http://wdeng.info/codes/alcon>.

Report bugs, request new features at:
<https://github.com/wenjundeng/alcon/issues>.

Contribute patches, new features by forking the project on GitHub:
<https://github.com/wenjundeng/alcon>, applying your contributions to the
forked project, and then submitting a pull request at:
<https://github.com/wenjundeng/alcon/pulls>.


Acknowledgments
---------------

ALCON was originally developed when I was a graduate student at University of
California, Irvine, supported by the US Department of Energy (DOE) SciDAC
Center for Gyrokinetic Simulation of Energetic Particle Turbulence and
Transport (GSEP).  ALCON is being further developed when I am a research
physicist at Princeton Plasma Physics Laboratory, supported by the US DOE
SciDAC Center for Nonlinear Simulation of Energetic Particles in Burning
Plasmas (CSEP).  I'd also like to acknowledge useful discussions with Eric
Bass, Guoyong Fu, Zhihong Lin, Don Spong, Xin Wang, Zhixuan Wang, and
Huasen Zhang.

Wenjun Deng

