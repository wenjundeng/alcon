! Copyright 2011, 2012 Wenjun Deng <wdeng@wdeng.info>
!
! This file is part of ALCON
!
! ALCON is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! ALCON is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with ALCON.  If not, see <http://www.gnu.org/licenses/>.

module alcon_input
implicit none
#include "finclude/petscdef.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the kind # of PetscReal !
! no need to change           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a real constant to test the kind # of PetscReal
PetscReal, parameter :: testkindPetscReal = 0d0
! kind # of PetscReal
PetscInt, parameter :: kpr = kind(testkindPetscReal)


!!!!!!!!!!!!!!!!!!!!
! equilibrium type !
!!!!!!!!!!!!!!!!!!!!
! eqtype = 20: load numerical equilibrium data from alcon.dat
! eqtype = 40: use a simple analytic equilibrium model described in Appendix E
!   in [Nuclear Fusion 52, 023005 (2012)](http://wdeng.info/?p=76) (extended to
!   second order), load profiles from profile_sa.dat
! eqtype = 60: load numerical equilibrium from spdata.dat and profile.dat (not
!   implemented yet)
PetscInt, parameter :: eqtype = 20


!!!!!!!!!!!!!!!!!!!!!
! solver parameters !
!!!!!!!!!!!!!!!!!!!!!
! finitebeta = 0: pure Alfven continua, no acoustic coupling
! finitebeta = 1: Alfven continua with acoustic coupling with slow sound
!   approximation [Physics of Fluids B: Plasma Physics 4, 3713 (1992)]
! finitebeta = 2: pure sound continua
! finitebeta = 3: accurate Alfven continua with acoustic coupling
PetscInt, parameter :: finitebeta = 3

! solving radial range (usually in the normalized radial coordinate rho)
! rho is square root of normalized toroidal flux:
!   rho = sqrt(psi_tor / psi_tor_wall)
! if eqtype = 20, the radial coordinate is determined by alcon.dat
! if eqtype = 40, the radial coordinate is rho, i.e., r/a in this simple
!   equilibrium
! if eqtype = 60, the radial coordinate is determined by profile.dat
!   note that if rad1 == 0, on the first radial grid, all matrix elements are
!   0, so there would be no output for the first radial grid
PetscReal, parameter :: rad1 = 0.01_kpr, rad2 = 1.0_kpr

! number of radial grid points
PetscInt, parameter :: nrad = 2000

! toroidal mode number n
PetscInt, parameter :: n = 4

! poloidal mode numbers from m1 to m2
PetscInt, parameter :: m1 = -20, m2 = 47

! # of off-diagonal stripes, representing cut-off of m-harmonic coupling
! noffdiag = 0: no coupling, trivial solution: \omega = v_A k_\parallel
! noffdiag = 1: m-harmonic only couples with (m - 1) and (m + 1), TAE gap shows
!   up
! noffdiag = 2: m-harmonic only couples with (m - 1), (m + 1), (m - 2), and
!   (m + 2)
PetscInt, parameter :: noffdiag = 20

! # of MPI processes for each radial grid point
! currently LAPACK is used for SLEPc solver, because when an eigenvalue is near
!   0, LAPACK gives good solution. However, when LAPACK is used, increasing #
!   of processes would not speed up.
PetscInt, parameter :: npe_solver = 1

! # of lowest Alfven continua are solved
PetscInt, parameter :: nev = 10

! for a solved omega, if the ratio of the absolute value of the imaginary part
!   to the real part is larger than this value, drop this solution due to too
!   large error.  because the Alfven continuum equation is Hermitian, the
!   imaginary part of the eigenvalue should be 0.
PetscReal, parameter :: imagrealratiocutoff = 0.1_kpr


!!!!!!!!!!!!!!!!!!!!!
! output parameters !
!!!!!!!!!!!!!!!!!!!!!
! scale omega before output
! 1.0 for no scaling, in which case omega is normalized by v_{Ap}/R_0, i.e.,
!   Eq. (A.36) in [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117)
PetscReal, parameter :: omegascale = 1.0_kpr
!PetscReal, parameter :: omegascale = 7.8384e8_kpr / &
!(173.86_kpr * 2.0_kpr * 3.14159265358979323846264_kpr * 1e3_kpr) ! DIII-D #142111 745ms, convert to kHz

! if scaled solution of omega is larger than this cut-off frequency, then drop
!   this solution.
! set to negative to disable cutting-off, in which case it is recommended to
!   set a vertical plotting range in your plotting script to avoid making a
!   figure dominated by continua of extremely high frequencies.
PetscReal, parameter :: omegacutoff = 1.0_kpr

! output format
! outformat = 0: generate a set of .dat files, each has n and m numbers in its
!   file name.  each file has two columns, first one for radial coordinate
!   (usually rho), second one for omega (continuum angular frequency).
! outformat = 10: in addition to .dat files, a MATLAB script file (alcon.m) is
!   generated for plotting in MATLAB.  tested with MATLAB 7.8.
! outformat = 20: in addition to .dat files, an IDL script file (alcon.pro) is
!   generated for plotting in IDL.  in IDL, simply executing "alcon" would plot
!   on screen; executing "alcon, /eps" would generate an alcon.eps file.  tested
!   with IDL 7.0.  note that IDL has no intrinsic function for making a legend
!   for a plot before 8.0, this outformat will put labels at their data average
!   positions.  as a result, the label positions will probably be not very
!   good.  you may want to open the alcon.pro file to adjust the label
!   positions.  they are near the end of the alcon.pro file.
! outformat = 30: in addition to .dat files, a gnuplot script file (alcon.gp)
!   and a Makefile (for GNU make) are generated for plotting using gnuplot.
!   execute "make" or "make alcon.eps" in the output directory will generate a
!   alcon.eps figure.  tested with gnuplot 4.4 and 4.6.
! outformat = 31: the same as outformat = 30, except that the generated
!   Makefile will make use of Omega Tools
!   <http://wdeng.info/codes/omegatools/>.
! outformat = 32: the same as outformat = 30, except that the gnuplot script
!   will use metapost terminal for output.  a modern LaTeX distribution is
!   required for this outformat.  the default figure from "make" is alcon.mps.
!   you can also use "make alcon.eps" to generate .eps file (requires mps2eps
!   <http://www.ctan.org/tex-archive/support/mps2eps>), or "make alcon.pdf" to
!   generate .pdf file, or "make alcon.png" to generate .png file (requires
!   "convert" in ImageMagick <http://www.imagemagick.org>).
! outformat = 33: the same as outformat = 32, except that the generated
!   Makefile will make use of Omega Tools
!   <http://wdeng.info/codes/omegatools/>.
PetscInt, parameter :: outformat = 10

! labels to put in figure
character(len=*), parameter :: &
  out_xlabel = "square root of normalized toroidal flux", &
  out_ylabel = "normalized angular frequency"
! if outformat = 33, all labels are in LaTeX math mode, so use something like
!   these instead:
!character(len=*), parameter :: &
!  out_xlabel = "\rho \textrm{ (square root of normalized toroidal flux)}", &
!  out_ylabel = "\omega / (v_{\mathrm{Ap}} / R_0)"

! starting m number that has labeling in output figure
PetscInt, parameter :: out_label_m1 = 5

! plotting radial range
PetscReal, parameter :: out_rad1 = 0.0_kpr, out_rad2 = 1.0_kpr

! verbose mode, controlling how much information is shown to stdout
! verbosity = 0: not showing any information except warnings and errors
! verbosity = 1: besides warnings and errors, show general information about
!   the program
! verbosity = 2: show additional diagnositc information (only for debugging, a
!   lot of variables printed out)
! verbosity = 3: show more additional diagnositc information (again, only for
!   debugging)
PetscInt, parameter :: verbosity = 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! debugging parameters               !
! only for development and debugging !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! radial grid to monitor
PetscInt, parameter :: iraddiag = 10

end module alcon_input

