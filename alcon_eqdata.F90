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

module alcon_eqdata
use alcon_input
implicit none
#include "finclude/petscdef.h"

! # of profiles and FFT data sets
PetscInt, parameter :: nprofiledata = 5, nfftdata = 5

! these variables store essential equilibrium data for alcon_solver
! to assign values to them, call one of the loadeq_* subroutines
! then they are ready to be used by alcon_solver

! index for nprofiledata:
! 1: radial coordinate (usually rho); 2: q; 3: g q + I; 4: beta; 5: rho_M
PetscReal, dimension(nprofiledata, nrad) :: acdprofile = 0.0_kpr
! correspondence between index for nfftdata and matrix:
! 1: H; 2: J; 3: K; 4: L; 5: N
! for details about these matrixes, see Eqs. (A.23)--(A.26) and (A.42) in [Nulcear Fusion 52, 043006 (2012)]
PetscScalar, dimension(0 : noffdiag, nfftdata, nrad) :: acdfft = (0.0_kpr, 0.0_kpr)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!
! generate radial grids !
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alcon_eqdata_radgrids
implicit none

PetscInt :: irad

if (nrad == 1) then ! only one grid point, then take it on rad1
  acdprofile(1, 1) = rad1
else ! more than one grid points, then uniform grids ranging [rad1, rad2]
  acdprofile(1, :) = (/ (rad1 + (rad2 - rad1) * (irad - 1.0_kpr) / (max(nrad, 2) - 1.0_kpr), irad = 1, nrad) /)
end if
end subroutine alcon_eqdata_radgrids


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load simple analytic equilibrium !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load simple analytic equilibrium model described in Appendix E in [Nuclear Fusion 52, 023005 (2012)] (extended to second order), load profiles from profile_sa.dat
! this second order equilibrium is:
! B = 1 - epsilon * cos(theta0) + epsilon**2 / 2 * (1/q**2 + cos(2 theta0)) + O(epsilon**3)
!   = 1 - epsilon * cos(theta) + epsilon**2 / 2 * (1/q**2 + 1) + O(epsilon**3)
! g = 1 - epsilon**2 / 2 + O(epsilon**3)
! I = epsilon**2 / q + epsilon**4 / (2 q) + O(epsilon**5)
! delta = epsilon * sin(theta0) - epsilon**2 / 2 * sin(2 theta0) + O(epsilon**3)
! theta = theta0 - epsilon * sin(theta0) + epsilon**2 / 4 * sin(2 theta0) + O(epsilon**3)
! theta0 = theta + epsilon * sin(theta) + epsilon**2 / 4 * sin(2 theta) + O(epsilon**3)
! zeta = zeta0 + O(epsilon**5)
subroutine alcon_eqdata_load_sa(iacerr)
use cspline_petsc
implicit none
#include "finclude/petsc.h90"

PetscInt, intent(out) :: iacerr ! error code

PetscInt, parameter :: nprofiledata_sa = 4 ! required # of profiles in profile_sa.dat
! file handle number
PetscInt, parameter :: fprofile_sa = 221

! MPIU_REAL for MPI subroutines manipulating PetscReal type
#if defined (PETSC_USE_REAL_SINGLE)
integer, parameter :: MPIU_REAL = MPI_REAL
#elif defined (PETSC_USE_REAL_LONG_DOUBLE)
integer, parameter :: MPIU_REAL = MPI_2DOUBLE_PRECISION
#else
integer, parameter :: MPIU_REAL = MPI_DOUBLE_PRECISION
#endif

PetscInt :: mype ! rank of current MPI process
PetscReal :: aminor ! minor radius normalized to major radius (a/R0)
PetscInt :: nrad_sa ! # of radial grid points read from profile_sa.dat
PetscInt :: nprofiledata1_sa ! # of profiles in profile_sa.dat
PetscInt :: iprofile, ifft ! for indexing
PetscReal, dimension(:), allocatable :: profile_sa ! profile data read from profile_sa.dat
PetscReal, dimension(:), allocatable :: d2profile ! 2nd derivative of profile data, for cspline

PetscErrorCode :: ierr
character(400) :: msg
PetscInt, dimension(2) :: intbuf

iacerr = 0

call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
CHKERRQ(ierr)
if (mype == 0) then
  open (fprofile_sa, file = 'profile_sa.dat', status = 'old')
  read (fprofile_sa, *) aminor
  read (fprofile_sa, *) nrad_sa, nprofiledata1_sa
  if (verbosity >= 2) then
    ! simply use write since mype == 0
    write (msg, *) "[alcon_eqdata_load_sa] Info: aminor = ", aminor
    write (*, "(a)") trim(adjustl(msg))
    write (msg, *) "[alcon_eqdata_load_sa] Info: nrad_sa = ", nrad_sa
    write (*, "(a)") trim(adjustl(msg))
    write (msg, *) "[alcon_eqdata_load_sa] Info: nprofiledata1_sa = ", nprofiledata1_sa
    write (*, "(a)") trim(adjustl(msg))
  end if
  if (nprofiledata1_sa /= nprofiledata_sa) then
    ! # of profiles in profile_sa.dat not correct
    iacerr = 1
  else
    allocate (profile_sa(nprofiledata_sa * nrad_sa))
    read (fprofile_sa, *) profile_sa
  end if
  close (fprofile_sa)
end if ! if (mype == 0)
if (verbosity >= 2) then
  write (msg, *) "[alcon_eqdata_load_sa] Info: iacerr = ", iacerr, "\n"
  call PetscPrintf(MPI_COMM_WORLD, adjustl(msg), ierr)
  CHKERRQ(ierr)
end if
intbuf(1) = iacerr
intbuf(2) = nrad_sa
! broadcast error code
call MPI_Bcast(intbuf, 2, MPIU_INTEGER, 0, MPI_COMM_WORLD, ierr)
CHKERRQ(ierr)
iacerr = intbuf(1)
nrad_sa = intbuf(2)
if (verbosity >= 3) then
  call PetscPrintf(MPI_COMM_WORLD, "[alcon_eqdata_load_sa] Info: successfully broadcasted iacerr and nrad_sa.\n", ierr)
  CHKERRQ(ierr)
end if
if (iacerr > 0) then
  return
end if
! broadcast aminor
call MPI_Bcast(aminor, 1, MPIU_REAL, 0, MPI_COMM_WORLD, ierr)
CHKERRQ(ierr)

if (mype /= 0) then
  allocate (profile_sa(nprofiledata_sa * nrad_sa))
end if
! broadcast profile_sa
call MPI_Bcast(profile_sa, nprofiledata_sa * nrad_sa, MPIU_REAL, 0, MPI_COMM_WORLD, ierr)
CHKERRQ(ierr)
if (verbosity >= 3) then
  call PetscPrintf(MPI_COMM_WORLD, "[alcon_eqdata_load_sa] Info: successfully broadcasted profile_sa.\n", ierr)
  CHKERRQ(ierr)
end if

! calculate acdprofile
call alcon_eqdata_radgrids ! generate radial grids (acdprofile(1, :))
if (verbosity >= 2) then
  call PetscPrintf(MPI_COMM_WORLD, "[alcon_eqdata_load_sa] Info: successfully gene", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "rated radial grids (alcon_eqdata_radgrids).\n", ierr)
  CHKERRQ(ierr)
end if
! initialize cspline_petsc
call cspline_init(nrad_sa)
cspline_warning_extrapolation = 1 ! let cspline_petsc to show warning of extrapolation only once
if (verbosity >= 2) then
  call PetscPrintf(MPI_COMM_WORLD, "[alcon_eqdata_load_sa] Info: successfully initialized cspline_petsc.\n", ierr)
  CHKERRQ(ierr)
end if
allocate (d2profile(nrad_sa))
if (verbosity >= 2) then
  call PetscPrintf(MPI_COMM_WORLD, "[alcon_eqdata_load_sa] Info: successfully allocated d2profile.\n", ierr)
  CHKERRQ(ierr)
end if
! rho = profile_sa(1 : nprofiledata_sa * nrad_sa : nprofiledata_sa)
! q = profile_sa(2 : nprofiledata_sa * nrad_sa : nprofiledata_sa)
! beta# = profile_sa(3 : nprofiledata_sa * nrad_sa : nprofiledata_sa)
! beta# = 4 pi gamma P_0 / B_a**2
! rho_M = profile_sa(4 : nprofiledata_sa * nrad_sa : nprofiledata_sa)
do iprofile = 2, nprofiledata_sa
  call cspline_construct_real( &
    profile_sa(1 : nprofiledata_sa * nrad_sa : nprofiledata_sa), &
    profile_sa(iprofile : nprofiledata_sa * nrad_sa : nprofiledata_sa), &
    d2profile(:), 0, 0.0_kpr, 0.0_kpr &
  )
  call cspline_interpolate_inc_real( &
    profile_sa(1 : nprofiledata_sa * nrad_sa : nprofiledata_sa), &
    profile_sa(iprofile : nprofiledata_sa * nrad_sa : nprofiledata_sa), &
    d2profile(:), acdprofile(1, :), acdprofile(iprofile, :) &
  )
end do
! shift acdprofile(3 : 4, :) to acdprofile(4 : 5, :)
acdprofile(5, :) = acdprofile(4, :)
acdprofile(4, :) = acdprofile(3, :)
! acdprofile(3, :) = g q + I
! g = 1 - epsilon**2/2 + O(epsilon**3)
! epsilon = r/R0 = acdprofile(1, :) * aminor
! I = (epsilon**2 + epsilon**4/2)/q + O(epsilon**5)
acdprofile(3, :) = (1.0_kpr - (acdprofile(1, :) * aminor)**2 / 2.0_kpr) &
  * acdprofile(2, :) + ((acdprofile(1, :) * aminor)**2 + (acdprofile(1, :) * aminor)**4 / 2.0_kpr) / acdprofile(2, :)

call cspline_final
deallocate (d2profile)

! calculate acdfft, i.e., matrix elements in Eqs. (A.23)--(A.26) and (A.43) in [Nuclear Fusion 52, 043006 (2012)]
! acdfft should be initialized to be 0 in the main program
! for matrix H, Eq. (A.23)
! H = (epsilon/q)**2 / (g q + I)
acdfft(0, 1, :) = cmplx((acdprofile(1, :) * aminor / acdprofile(2, :))**2 / acdprofile(3, :), 0.0_kpr, kpr)

! for matrix J, Eq. (A.24)
! J = (epsilon/q)**2 * (g q + I) / B^4
! 1/B^4 = 1 + epsilon 4 cos(theta) + epsilon**2 * 2 * (5 cos(theta)**2 - 1 - 1/q**2) + O(epsilon**3)
!       = (1 + epsilon**2 * (3 - 2 / q**2))
!       + 2 epsilon (exp(i theta) + exp(-i theta))
!       + 5/2 epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
do ifft = 0, min(2, noffdiag)
  acdfft(ifft, 2, :) = cmplx((acdprofile(1, :) * aminor / acdprofile(2, :))**2 * acdprofile(3, :), 0.0_kpr, kpr) ! (epsilon/q)**2 * (g q + I)
end do
! multiply by 1/B^4
acdfft(0, 2, :) = acdfft(0, 2, :) * (1 + (acdprofile(1, :) * aminor)**2 * (3.0_kpr - 2.0_kpr / acdprofile(2, :)**2))
if (noffdiag >= 1) acdfft(min(noffdiag, 1), 2, :) = acdfft(min(noffdiag, 1), 2, :) * 2.0_kpr * acdprofile(1, :) * aminor
if (noffdiag >= 2) acdfft(min(noffdiag, 2), 2, :) = acdfft(min(noffdiag, 2) ,2, :) * 2.5_kpr * (acdprofile(1, :) * aminor)**2
! the min functions above are for suppressing the "out of bound" warning
! so are the min functions in the rest of this subroutine

! for matrix K, Eq. (A.25)
! K = 2 g / B^3 * (partial B/partial theta)
!   = 2 epsilon * sin(theta) + 3 epsilon**2 * sin(2 theta) + O(epsilon**3)
!   = -i epsilon (exp(i theta) - exp(-i theta)) - i 3/2 epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
acdfft(0, 3, :) = (0.0_kpr, 0.0_kpr)
if (noffdiag >= 1) acdfft(min(noffdiag, 1), 3, :) = cmplx(0.0_kpr, - acdprofile(1, :) * aminor, kpr)
if (noffdiag >= 2) acdfft(min(noffdiag, 2), 3, :) = cmplx(0.0_kpr, -1.5_kpr * (acdprofile(1, :) * aminor)**2, kpr)

! for matrix L, Eq. (A.26)
! 1 / B^2 = 1 + epsilon 2 cos(theta) + epsilon**2 (3 cos(theta)**2 - 1 - 1/q**2) + O(epsilon**3)
!         = 1 + epsilon**2 (1/2 - 1/q**2)
!         + epsilon (exp(i theta) + exp(-i theta))
!         + 3/4 epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
! 1 / B^4 = (1 + epsilon**2 * (3 - 2 / q**2))
!         + 2 epsilon (exp(i theta) + exp(-i theta))
!         + 5/2 epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
! L = (g q + I) (beta# / B^4 + 1 / B^2)
!   = (g q + I) * (
!       beta# (1 + epsilon**2 * (3 - 2 / q**2)) + 1 + epsilon**2 (1/2 - 1/q**2)
!     + (2 beta# + 1) epsilon (exp(i theta) + exp(-i theta))
!     + (5/2 beta# + 3/4) epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
!     )
do ifft = 0, min(2, noffdiag)
  acdfft(ifft, 4, :) = cmplx(acdprofile(3, :), 0.0_kpr, kpr)
end do
acdfft(0, 4, :) = acdfft(0, 4, :) * ( &
  acdprofile(4, :) * (1.0_kpr + (acdprofile(1, :) * aminor)**2 * (3.0_kpr - 2.0_kpr / acdprofile(2, :)**2)) &
  + 1.0_kpr + (acdprofile(1, :) * aminor)**2 * (0.5_kpr - 1.0_kpr / acdprofile(2, :)**2) &
)
if (noffdiag >= 1) acdfft(min(noffdiag, 1), 4, :) = acdfft(min(noffdiag, 1), 4, :) &
  * (2.0_kpr * acdprofile(4, :) + 1.0_kpr) * (acdprofile(1, :) * aminor)
if (noffdiag >= 2) acdfft(min(noffdiag, 2), 4, :) = acdfft(min(noffdiag, 2), 4, :) &
  * (2.5_kpr * acdprofile(4, :) + 0.75_kpr) * (acdprofile(1, :) * aminor)**2

! for matrix N, Eq. (A.43)
! N = 4 beta# g**2 / (g q + I) * (partial B/partial theta)**2 / ((beta# + B^2) B^2)
!   = 4 beta# g**2 / (g q + I) * epsilon**2 * (
!       sin(theta)**2 / (beta# + 1)
!     + epsilon 2 (beta# + 2) sin(theta)**2 * cos(theta) / (beta# + 1)^2 + O(epsilon**2)
!     )
!   = beta# g**2 / (g q + I) * epsilon**2 * (
!       2 / (beta# + 1)
!     + epsilon (beta# + 2) / (beta# + 1)**2 * (exp(i theta) + exp(-i theta))
!     - 1 / (beta# + 1) * (exp(i 2 theta) + exp(-i 2 theta))
!     - epsilon (beta# + 2) / (beta# + 1)**2 * (exp(i 3 theta) + exp(-i 3 theta))
!     )
do ifft = 0, min(3, noffdiag)
  acdfft(ifft, 5, :) = cmplx(acdprofile(4, :) * (1.0_kpr - (acdprofile(1, :) * aminor)**2) &
    / acdprofile(3, :) * (acdprofile(1, :) * aminor)**2, 0.0_kpr, kpr)
end do
acdfft(0, 5, :) = acdfft(0, 5, :) * 2.0_kpr / (acdprofile(4, :) + 1.0_kpr)
if (noffdiag >= 1) acdfft(min(noffdiag, 1), 5, :) = acdfft(min(noffdiag, 1), 5, :) &
  * (acdprofile(1, :) * aminor) * (acdprofile(4, :) + 2.0_kpr) / (acdprofile(4, :) + 1.0_kpr)**2
if (noffdiag >= 2) acdfft(min(noffdiag, 2), 5, :) = -acdfft(min(noffdiag, 2), 5, :) / (acdprofile(4, :) + 1.0_kpr)
if (noffdiag >= 3) acdfft(min(noffdiag, 3), 5, :) = -acdfft(min(noffdiag, 3), 5, :) &
  * (acdprofile(1, :) * aminor) * (acdprofile(4, :) + 2.0_kpr) / (acdprofile(4, :) + 1.0_kpr)**2

deallocate (profile_sa)

end subroutine alcon_eqdata_load_sa


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load equilibrium data from alcon.dat !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! alcon.dat can be generated from GTC electromagnetic run with fload=0
subroutine alcon_eqdata_load_alcondat(iacerr)
use cspline_petsc
implicit none
#include "finclude/petsc.h90"

PetscInt, intent(out) :: iacerr ! error code

! file handler number
PetscInt, parameter :: falcondat = 214

! MPIU_REAL for MPI subroutines manipulating PetscReal type
#if defined (PETSC_USE_REAL_SINGLE)
integer, parameter :: MPIU_REAL = MPI_REAL
#elif defined (PETSC_USE_REAL_LONG_DOUBLE)
integer, parameter :: MPIU_REAL = MPI_2DOUBLE_PRECISION
#else
integer, parameter :: MPIU_REAL = MPI_DOUBLE_PRECISION
#endif

PetscReal, dimension(:), allocatable :: acdprofile_alcondat, d2profile
PetscScalar, dimension(:), allocatable :: acdfft_alcondat, d2fft
PetscInt :: mype, ifft, ifftdata, iprofile
PetscInt :: nrad_alcondat, nfft_alcondat, nprofiledata_alcondat, nfftdata_alcondat

PetscErrorCode :: ierr
PetscInt, dimension(3) :: intbuf
character(400) :: msg

iacerr = 0

nrad_alcondat = 0
nfft_alcondat = 0
nprofiledata_alcondat = 0
nfftdata_alcondat = 0
call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
CHKERRQ(ierr)
if (mype == 0) then
  open (falcondat, file = 'alcon.dat', status = 'old')
  read (falcondat, *) nrad_alcondat, nfft_alcondat, nprofiledata_alcondat, nfftdata_alcondat
  if (verbosity >= 2) then
    write (msg, *) "[alcon_eqdata_load_alcondat] Info: nrad_alcondat = ", nrad_alcondat
    write (*, "(a)") trim(adjustl(msg))
    write (msg, *) "[alcon_eqdata_load_alcondat] Info: nfft_alcondat = ", nfft_alcondat
    write (*, "(a)") trim(adjustl(msg))
    write (msg, *) "[alcon_eqdata_load_alcondat] Info: noffdiag = ", noffdiag
    write (*, "(a)") trim(adjustl(msg))
    write (msg, *) "[alcon_eqdata_load_alcondat] Info: nprofiledata_alcondat = ", nprofiledata_alcondat
    write (*, "(a)") trim(adjustl(msg))
    write (msg, *) "[alcon_eqdata_load_alcondat] Info: nfftdata_alcondat = ", nfftdata_alcondat
    write (*, "(a)") trim(adjustl(msg))
  end if
  if (nprofiledata_alcondat /= nprofiledata .or. nfftdata_alcondat /= nfftdata) then
    iacerr = 1
  else
    allocate (acdprofile_alcondat(nprofiledata * nrad_alcondat))
    allocate (acdfft_alcondat((1 + nfft_alcondat) * nfftdata * nrad_alcondat))
    read (falcondat, *) acdprofile_alcondat
    read (falcondat, *) acdfft_alcondat
  end if
  close (falcondat)
end if
intbuf(1) = iacerr
intbuf(2) = nrad_alcondat
intbuf(3) = nfft_alcondat
! broadcast error code and other integer parameters
call MPI_Bcast(intbuf, 3, MPIU_INTEGER, 0, MPI_COMM_WORLD, ierr)
CHKERRQ(ierr)
iacerr = intbuf(1)
nrad_alcondat = intbuf(2)
nfft_alcondat = intbuf(3)
if (iacerr > 0) then
  if (mype == 0) deallocate(acdprofile_alcondat, acdfft_alcondat)
  return
end if
if (mype /= 0) then
  allocate (acdprofile_alcondat(nprofiledata * nrad_alcondat))
  allocate (acdfft_alcondat((1 + nfft_alcondat) * nfftdata * nrad_alcondat))
end if
! broadcast acdprofile_alcondat and acdfft_alcondat
call MPI_Bcast(acdprofile_alcondat, nprofiledata * nrad_alcondat, MPIU_REAL, 0, MPI_COMM_WORLD, ierr)
CHKERRQ(ierr)
call MPI_Bcast(acdfft_alcondat, (1 + nfft_alcondat) * nfftdata * nrad_alcondat, MPIU_SCALAR, 0, MPI_COMM_WORLD, ierr)
CHKERRQ(ierr)

! initialized cspline_petsc
call cspline_init(nrad_alcondat)
cspline_warning_extrapolation = 1 ! let cspline_petsc to show warning of extrapolation only once
allocate (d2profile(nrad_alcondat), d2fft(nrad_alcondat))

! calculate acdprofile 
! interpolate using acdprofile_alcondat
call alcon_eqdata_radgrids ! generate radial grids (acdprofile(1, :))
do iprofile = 2, nprofiledata
  call cspline_construct_real( &
    acdprofile_alcondat(1 : nprofiledata * nrad_alcondat : nprofiledata), &
    acdprofile_alcondat(iprofile : nprofiledata * nrad_alcondat : nprofiledata), &
    d2profile(:), 0, 0.0_kpr, 0.0_kpr &
  )
  call cspline_interpolate_inc_real( &
    acdprofile_alcondat(1 : nprofiledata * nrad_alcondat : nprofiledata), &
    acdprofile_alcondat(iprofile : nprofiledata * nrad_alcondat : nprofiledata), &
    d2profile(:), acdprofile(1, :), acdprofile(iprofile, :) &
  )
end do
! direct copy (for debug only; if use, make sure nrad == nrad_alcondat)
!do iprofile = 1, nprofiledata
!  acdprofile(iprofile, :) = acdprofile_alcondat(iprofile : nprofiledata * nrad_alcondat : nprofiledata)
!end do

! calculate acdfft 
! interpolate using acdfft_alcondat
do ifft = 0, min(nfft_alcondat, noffdiag)
  do ifftdata = 1, nfftdata
    call cspline_construct( &
      acdprofile_alcondat(1 : nprofiledata * nrad_alcondat : nprofiledata), &
      acdfft_alcondat( &
        (ifft + 1) + nfft_alcondat * (ifftdata - 1) &
        : (nfft_alcondat + 1) * nfftdata * nrad_alcondat : (nfft_alcondat + 1) * nfftdata &
      ), d2fft(:), 0, (0.0_kpr, 0.0_kpr), (0.0_kpr, 0.0_kpr) &
    )
    call cspline_interpolate_inc( &
      acdprofile_alcondat(1 : nprofiledata * nrad_alcondat : nprofiledata), &
      acdfft_alcondat( &
        (ifft + 1) + (nfft_alcondat + 1) * (ifftdata - 1) &
        : (nfft_alcondat + 1) * nfftdata * nrad_alcondat : (nfft_alcondat + 1) * nfftdata &
      ), d2fft(:), acdprofile(1, :), acdfft(ifft, ifftdata, :) &
    )
  end do
end do
! direct copy (for debug only; if use, make sure nrad == nrad_alcondat)
!do ifft = 0, min(nfft_alcondat, noffdiag)
!  do ifftdata = 1, nfftdata
!    acdfft(ifft, ifftdata, :) = acdfft_alcondat( &
!      (ifft + 1) + (nfft_alcondat + 1) * (ifftdata - 1) &
!      : (nfft_alcondat + 1) * nfftdata * nrad_alcondat : (nfft_alcondat + 1) * nfftdata &
!    )
!  end do
!end do
if (noffdiag > nfft_alcondat) then
  call PetscPrintf(MPI_COMM_WORLD, "[alcon_eqdata_load_alcondat] Warning: nfft_alcondat < noffdiag, ", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "not enough data to calculate specified ordering of coupling.\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "You may want to re-generate your alcon.dat file to contain ", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "higher m-harmonic amplitude data.\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "(increase # of poloidal grid points for FFT and keep more ", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "elements of the transformed array).\n", ierr)
  CHKERRQ(ierr)
  acdfft(nfft_alcondat + 1 : noffdiag, :, :) = (0.0_kpr, 0.0_kpr)
end if

call cspline_final
deallocate(acdprofile_alcondat, acdfft_alcondat, d2profile, d2fft)

end subroutine alcon_eqdata_load_alcondat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! load equilibrium from spdata.dat and profile.dat !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! not implemented yet
subroutine loadeq_spdata(iacerr)
! gamma for ions and electrons 
PetscReal, parameter :: gammai = 7.0_kpr / 4.0_kpr, gammae = 1.0_kpr
! number of poloidal (theta) grid points in the sampling for fft
PetscInt, parameter :: npol = 2560

PetscErrorCode :: iacerr

iacerr = 0

end subroutine loadeq_spdata

end module alcon_eqdata

