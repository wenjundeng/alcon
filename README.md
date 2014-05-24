ALCON
=====
version 2014-05-24 13:22:42-04:00
---------------------------------

ALCON has been rewritten in Python.  All users including those of the former
Fortran version are encouraged to read through this README to get familiar with
the new Python version.

ALCON is a code for solving ideal MHD Alfven continua in tokamaks, i.e.,
Eq. (10) in [Physics of Fluids 29, 3695 (1986)][], using a poloidal-spectral
method described in Appendix A in [Nuclear Fusion 52, 043006 (2012)][].  When
referencing ALCON in a publication, please cite:

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

[Physics of Fluids 29, 3695 (1986)]: http://dx.doi.org/10.1063/1.865801
[Nuclear Fusion 52, 043006 (2012)]: http://wdeng.info/?p=117


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

ALCON requires:

+ [Python](https://www.python.org/) 2 (2.7+) or 3 (3.4+)
+ [NumPy](http://www.numpy.org/) 1.8+
+ [SciPy](http://www.scipy.org/scipylib/index.html) 0.13+

Additionally, [matplotlib](http://matplotlib.org/) 1.3+ is recommended for
quickly plotting solved continua.

Unless you keep excessive (100+) poloidal m-harmonics and want excessive high
radial resolution (solving 5000+ radial grid points), ALCON should be able to
finish a run within half an hour on any desktop that has two or more physical
cores.  Running the default case on an Apple Mac Mini 2011 takes about 14
minutes with a single process and about 10 minutes with two processes.


Usage
-----

### Know about the files

+ `README.md` -- This instruction.
+ `COPYING` -- The license (GNU General Public License) for ALCON.
+ `alcon.sh` -- A bash script for conveniently launching alcon.py without the
need to type all the input parameters in the command line.  This script can be
used as an input file.
+ `alcon.py` -- The main program file.
+ `alcon_input.py` -- Module for managing input parameters.
+ `alcon_eqdata.py` -- Module for managing equilibrium data.
+ `alcon_solver.py` -- Module for solving Alfven continua.
+ `alcon_output.py` -- Module for writing solutions to output files.
+ `wctimer.py` -- Module for a wall-clock timer.
+ `alcon_plot.py` -- Plotting utility for plotting output data from ALCON.
+ `alcon_plot_input.py` -- Module for managing input parameters for the
plotting utility.
+ `alcon_plot_core.py` -- Module for core functionalities of the plotting
utility.
+ `alcon.dat` -- A sample equilibrium data input file.  The equilibrium is
taken from [Physics of Plasmas 17, 112504 (2010)][].  This file is for testing
when you use ALCON for the first time.  For practical use, you want to replace
this file by your own `alcon.dat` file generated from your equilibrium.  See
`acdgen_sample.F90` for a sample subroutine to generate `alcon.dat`.
+ `profile_sa.dat` -- A sample equilibrium profiles for solving a simple
analytic equilibrium.  The profiles are taken from
[Physics of Plasmas 17, 112504 (2010)][].  This file is for testing when you
use ALCON for the first time.  For practical use, you want to replace this file
by your own `profile_sa.dat` file.  See next section for preparing your own
`profile_sa.dat`.
+ `acdgen_sample.F90` -- A sample subroutine to show how to generate
`alcon.dat`.

[Physics of Plasmas 17, 112504 (2010)]: http://wdeng.info/?p=37


### Running

The main program file is `alcon.py`.  To launch it, you can execute
`python alcon.py [options]`, or simply `./alcon.py [options]` if your system
supports direct script execution (with the script parser specified by the
[shebang][] in the script; supported by most UNIX and Linux systems including
Mac OS X) and the execution permission is granted for `alcon.py`.  To see all
the available `[options]`, use the `-h` option, i.e., execute
`python alcon.py -h`.  First-time user can simply run `python alcon.py` to try
out ALCON with default parameters and a sample tokamak equilibrium.  To get
meaningful results from the execution of `alcon.py`, careful input preparation
is needed.

Note that the above paragraph assumed the Python interpreter to be `python`.
If this is not the case, you need to replace the interpreter specifier
accordingly.  For example, if you use Python 3 and its interpreter is named
`python3` instead of `python`, you need to execute
`python3 alcon.py [options]` instead of `python alcon.py [options]`, or replace
the shebang line `#!/usr/bin/env python` by `#!/usr/bin/env python3` in
`alcon.py` for direct script execution.  This also applies to the plotting
utility `alcon_plot.py`.

[shebang]: http://en.wikipedia.org/wiki/Shebang_(Unix)


### Input preparation

ALCON (`alcon.py`) reads input parameters from the command line.  Execute
`python alcon.py -h` to see all the command line parameters.

A bash script `alcon.sh` is provided for easily setting input parameters.  You
can directly edit the input parameters in this file, and then execute this file
(`bash alcon.sh` or `./alcon.sh`) in the command line to launch ALCON with your
specified input parameters.

Depending on the parameter passed to the `--eqtype` option, ALCON reads
equilibrium data from additional input data file.  For currently implemented
parameters for `--eqtype`, ALCON may read from `alcon.dat` or `profile_sa.dat`
for equilibrium data.  ALCON comes with a sample of each of these two files.

To prepare your own `alcon.dat`, see `acdgen_sample.F90` for a sample
subroutine to generate `alcon.dat`.

To prepare your own `profile_sa.dat`, follow this format: the first line puts
one real number, which is the inverse aspect ratio a/R0; the second line puts
two integers `nrad` and `nprofile`, `nrad` being the number of radial grid
points for the profiles, `nprofile` being the number of profiles (for current
version always put `nprofile` to be 4); starting from the third line is a
`nrad * nprofile` matrix giving the profile data.  The first column is the
radial coordinate r/a.  The second column is the q-profile.  The third column
is the pressure profile normalized as: (4 pi gamma P / B^2), where gamma is the
specific heat ratio, P is the plasma pressure, B is the on-axis magnetic field,
and P and B are in CGS units.  The fourth column is the mass density normalized
by proton mass times on-axis electron density (Eq. (A.35) in
[Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117)).


### Output data and plotting

All output data files are put in the sub-directory specified by the `--dirout`
option.  Note that if the specified directory preexists, ALCON will exit with a
"file exists" error, unless the `--erasedirout` option is specified, in which
case the preexisting directory will be erased.

A plotting utility (`alcon_plot.py`) is provided for plotting and generating
plotting scripts for other visualization softwares.  It supports these
operations:

+ plotting on screen or saving a .png or .pdf figure through matplotlib;
+ generating a .m script for MATLAB;
+ generating a .pro script for IDL;
+ generating a .gp script and a Makefile for gnuplot.

Run `python alcon_plot.py -h` for detailed usage information.  First-time user
can simply run `python alcon_plot.py --dirout <directory of ALCON output>` to
get the continuum plot on screen, then gradually try out other options.


Known issue
-----------

ALCON can run in parallel by specifying an integer larger than 1 to the `-n`
option.  The parallelization is implemented using the multiprocessing module in
the Python Standard Library.  However, on Mac OS X, the subprocesses of a
program that is parallelized in this way and has certain NumPy and SciPy
function calls would crash, if NumPy and SciPy are built with
[Apple's Accelerate Framework][].  More information about this issue can be
found at these links:

+ <https://github.com/numpy/numpy/issues/654>
+ <http://mail.scipy.org/pipermail/numpy-discussion/2012-August/063590.html>

This issue does not apply if NumPy and SciPy are built with [OpenBLAS][].  If
you use [Homebrew][] to install NumPy and SciPy, you can choose to build them
with OpenBLAS by passing the `--with-openblas` option:

	brew install numpy --with-openblas
	brew install scipy --with-openblas

The above commands may complain:

	Error: Operation already in progress for openblas
	Another active Homebrew process is already using openblas.
	Please wait for it to finish or terminate it to continue.

To work around this, run `brew edit numpy`, then the NumPy formula would be
opened in the default text editor.  In the editor, replace
`homebrew/science/openblas` by `openblas`, save and quit.  Then run
`brew edit scipy` and do the same replacement.  Afterwards, run the above
commands to install NumPy and SciPy.

[Apple's Accelerate Framework]: https://developer.apple.com/library/mac/documentation/Accelerate/Reference/AccelerateFWRef/_index.html
[OpenBLAS]: http://www.openblas.net/
[Homebrew]: http://brew.sh/


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
California, Irvine, where my research was supported by the U.S. Department of
Energy (DOE) SciDAC Center for Gyrokinetic Simulation of Energetic Particle
Turbulence and Transport (GSEP).  The development continued after I became a
postdoc at Princeton Plasma Physics Laboratory, where my research was supported
by the U.S. DOE SciDAC Center for Nonlinear Simulation of Energetic Particles
in Burning Plasmas (CSEP).  I would like to acknowledge useful discussions with
Eric Bass, Guo-Yong Fu, Zhihong Lin, Don Spong, Xin Wang, Zhixuan Wang, and
Huasen Zhang.

Wenjun Deng
