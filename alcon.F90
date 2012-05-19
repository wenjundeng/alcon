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

program alcon
use alcon_input
use alcon_eqdata
use alcon_solver
use alcon_output
implicit none

#include "finclude/petsc.h90"
#include "finclude/slepc.h"

character(25), parameter :: version = '2012-05-19 18:43:25-04:00'

! indexing variables
PetscInt :: i
PetscInt :: irad, iprofile, ifft, ifftdata, iprogress

PetscInt, dimension(10) :: progress ! for tracking running progress
character(10) :: datetime(4) ! times when program starts and ends
PetscErrorCode :: ierr ! error code for MPI, PETSc, and SLEPc
PetscInt :: iacerr ! error code for ALCON

MPI_Comm :: comm_solver ! MPI communicator within each radial domain
MPI_Comm :: comm_radhead ! MPI communicator for radial domain heads
MPI_Group :: MPI_GROUP_WORLD ! the MPI group of all processes
MPI_Group :: group_radhead ! MPI group of radial domain heads
PetscInt :: npe ! # of total MPI processes
PetscInt :: mype ! MPI rank of current process
PetscInt :: ndom_rad ! # of radial domains
PetscInt :: mydom_rad ! radial domain that current process belongs to
PetscInt :: mype_solver ! MPI rank or current process in solver communicator
PetscInt, dimension(:), allocatable :: radhead ! MPI ranks of radial domain heads
PetscInt :: nrad_dom ! # of radial grid points per domain
PetscInt :: ndom_rad_mod ! # of radial domains that have one more radial grid point in addition to nrad_dom
PetscInt :: irad_low, irad_high ! lower and higher limits of radial grid index for each domain

character(2048) :: msg ! messeage to be displayed


#ifndef PETSC_USE_COMPLEX
call MPI_Init(ierr)
CHKERRQ(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
CHKERRQ(ierr)
if (mype == 0) then
  write (*, "(a)") "Error: ALCON requires complex version of PETSc. Please re-compiled it with complex version of PETSc."
endif
call MPI_Finalize(ierr)
CHKERRQ(ierr)
stop 1
#endif

call date_and_time(datetime(1), datetime(2))

! initialization
call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
CHKERRQ(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
CHKERRQ(ierr)

if (verbosity >= 1) then
  call PetscPrintf(MPI_COMM_WORLD, "ALCON version " // version // "\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "Info: MPI, PETSc and SLEPc successfully initialized.\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "Info: loading equilibrium...\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "Info: if a warning about extrapolation shows up in the follow", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "ing, it means your requested solving radial range is larger th", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "an the radial range provided by equilibrium data.\n", ierr)
  CHKERRQ(ierr)
endif

select case (eqtype)
  case (20)
    iacerr = 0
    call alcon_eqdata_load_alcondat(iacerr)
    if (iacerr == 1) then
      call PetscPrintf(MPI_COMM_WORLD, "[alcon_eqdata_load_alcondat] Error: number o", ierr)
      CHKERRQ(ierr)
      call PetscPrintf(MPI_COMM_WORLD, "f profiles and/or FFT datasets is not corr", ierr)
      CHKERRQ(ierr)
      call PetscPrintf(MPI_COMM_WORLD, "ect in alcon.dat.\n", ierr)
      CHKERRQ(ierr)
      call SlepcFinalize(ierr)
      CHKERRQ(ierr)
      stop 2
    endif

  case (40)
    iacerr = 0
    call alcon_eqdata_load_sa(iacerr)
    if (iacerr == 1) then
      call PetscPrintf(MPI_COMM_WORLD, "[alcon_eqdata_load_sa] Error: number of profi", ierr)
      CHKERRQ(ierr)
      call PetscPrintf(MPI_COMM_WORLD, "les is not correct in profile_sa.dat.\n", ierr)
      CHKERRQ(ierr)
      call SlepcFinalize(ierr)
      CHKERRQ(ierr)
      stop 2
    endif

  case default
    call PetscPrintf(MPI_COMM_WORLD, "Error: undefined equilibrium type (eqtype).\n", ierr)
    CHKERRQ(ierr)
    call SlepcFinalize(ierr)
    CHKERRQ(ierr)
    stop 1
endselect
if (verbosity >= 2) then
  call PetscPrintf(MPI_COMM_WORLD, "Diagnositc info: acdprofile(iprofile, :):", ierr)
  CHKERRQ(ierr)
  do iprofile = 1, nprofiledata
    write (msg, *) "iprofile = ", iprofile, ": ", &
      (acdprofile(iprofile, irad), irad = 1, min(nrad, 3)), "...", &
      acdprofile(iprofile, min(iraddiag, nrad)), "...", &
      (acdprofile(iprofile, irad), irad = max(nrad - 2, 1), nrad), "\n"
    call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
    CHKERRQ(ierr)
  enddo
  call PetscPrintf(MPI_COMM_WORLD, "Diagnostic info: acdfft(ifft, ifftdata, :):", ierr)
  CHKERRQ(ierr)
  do ifftdata = 1, nfftdata
    do ifft = 0, 2
      write (msg, *) "ifftdata = ", ifftdata, ", ifft = ", ifft, ": ", &
        (acdfft(ifft, ifftdata, irad), irad = 1, min(nrad, 3)), "...", &
        acdfft(ifft, ifftdata, min(iraddiag, nrad)), "...", &
        (acdfft(ifft, ifftdata, irad), irad = max(nrad - 2, 1), nrad), "\n"
      call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
      CHKERRQ(ierr)
    enddo
  enddo
endif

if (verbosity >= 1) then
  call PetscPrintf(MPI_COMM_WORLD, "Info: initializing solver, creating PETSc and SLEPc objects...\n", ierr)
  CHKERRQ(ierr)
endif

! radial domain and solver decomposition
call MPI_Comm_size(MPI_COMM_WORLD, npe, ierr)
CHKERRQ(ierr)
ndom_rad = ceiling(real(npe) / npe_solver)
if (ndom_rad > nrad) then
  ndom_rad = nrad
  call PetscPrintf(MPI_COMM_WORLD, "Warning: requested # of radial domains is lar", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "ger than # of radial grid points. Shrinking # o", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "f radial domains to # of radial grid points.", ierr)
  CHKERRQ(ierr)
endif

mydom_rad = mod(mype, ndom_rad)
mype_solver = mype / ndom_rad

call MPI_Comm_split(MPI_COMM_WORLD, mydom_rad, mype_solver, comm_solver, ierr)
CHKERRQ(ierr)

allocate (radhead(0 : ndom_rad - 1))
radhead(:) = (/ (i, i = 0, ndom_rad - 1) /)
call MPI_Comm_group(MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr)
CHKERRQ(ierr)
call MPI_Group_incl(MPI_GROUP_WORLD, ndom_rad, radhead, group_radhead, ierr)
CHKERRQ(ierr)
call MPI_Comm_create(MPI_COMM_WORLD, group_radhead, comm_radhead, ierr)
CHKERRQ(ierr)
deallocate (radhead)

nrad_dom = nrad / ndom_rad
ndom_rad_mod = mod(nrad, ndom_rad)
irad_low = mydom_rad * nrad_dom + 1 + min(mydom_rad, ndom_rad_mod)
irad_high = (mydom_rad + 1) * nrad_dom + min(mydom_rad + 1, ndom_rad_mod)
if (irad_high > nrad) irad_high = nrad
if (verbosity >= 2) then
  write (msg, *) "Info: nrad_dom = ", nrad_dom, ", ndom_rad_mod = ", ndom_rad_mod, "\n"
  call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
  CHKERRQ(ierr)
  write (msg, *) "Info: mype = ", mype, "irad_low = ", irad_low, "irad_high = ", irad_high
  write (*, "(a)") trim(adjustl(msg))
endif

! solver initialization, create PETSc and SLEPc objects
call alcon_solver_init(comm_solver)

if (verbosity >= 1) then
  call PetscPrintf(MPI_COMM_WORLD, "Info: solving Alfven continua...\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(MPI_COMM_WORLD, "(only the progress of the first radial domain will be shown)\n", ierr)
  CHKERRQ(ierr)
endif

if (verbosity == 1) then
  write (msg, *) "Info: solving radial grids (first radial domain total", irad_high - irad_low + 1, ")..."
  call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
  CHKERRQ(ierr)
  do iprogress = 1, 10
    progress(iprogress)  = iprogress * nrad_dom / 10
  enddo
  iprogress = 1
endif

! main loop over radial grids
do irad = irad_low, irad_high
!do irad = 1, nrad
!do irad = iraddiag - 1, iraddiag
  if (verbosity == 1) then
    if (irad >= progress(iprogress)) then
      write (msg, "(i2, a)") iprogress, "0%%..."
      call PetscPrintf(MPI_COMM_WORLD, msg, ierr)
      CHKERRQ(ierr)
      iprogress = iprogress + 1
    endif
  elseif (verbosity >= 2) then
    write (msg, *) "Info: solving radial grid: ", irad, "/", nrad_dom, "\n"
    call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
    CHKERRQ(ierr)
  endif

  ! assign values to matrixes A and B as specified in Eqs. (A.20) and (A.21) in
  !   [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117)
  call alcon_solver_prepare(irad)

  ! this doesn't work anymore after radial domain decomposition
  ! so here set condition to "verbosity >= 10" to disable it
  if (verbosity >= 10 .and. mydom_rad == 0) then
    ! print out matA, matB on irad == 1 and irad == iraddiag
    if (irad == 1 .or. irad == iraddiag) then
      call PetscPrintf(MPI_COMM_WORLD, "Info: matA:\n", ierr)
      CHKERRQ(ierr)
      call MatView(acs_matA, PETSC_VIEWER_STDOUT_WORLD, ierr)
      CHKERRQ(ierr)
      call PetscPrintf(MPI_COMM_WORLD, "Info: matB:\n", ierr)
      CHKERRQ(ierr)
      call MatView(acs_matB, PETSC_VIEWER_STDOUT_WORLD, ierr)
      CHKERRQ(ierr)
    endif
  endif

  ! solve for eigenvalues (continua frequencies) and output
  call alcon_solver_solve(acdprofile(1, irad), acdprofile(5, irad))
enddo
if (verbosity == 1) then
  call PetscPrintf(MPI_COMM_WORLD, "\n", ierr)
  CHKERRQ(ierr)
endif


if (verbosity >= 1) then
  call PetscPrintf(MPI_COMM_WORLD, "Info: writing data output files...\n", ierr)
  CHKERRQ(ierr)
endif
if (mype < ndom_rad) then ! radial domain heads write data output files
  call alcon_output_data(comm_radhead, ndom_rad)
endif

if (verbosity >= 2) then
  call PetscPrintf(MPI_COMM_WORLD, "Info: writing extra output files...\n", ierr)
  CHKERRQ(ierr)
endif
if (mype == 0) then
  call alcon_output_extra
endif

if (verbosity >= 1) then
  call PetscPrintf(MPI_COMM_WORLD, "Info: destroying objects and finalizing...\n", ierr)
  CHKERRQ(ierr)
endif

! solver finalization, destroy PETSc and SLEPc objects
call alcon_solver_final

call SlepcFinalize(ierr)
CHKERRQ(ierr)

if (verbosity >= 1 .and. mype == 0) then
  call date_and_time(datetime(3), datetime(4))
  write (*, "(4a)") 'Info: ALCON started at date:', datetime(1), 'time:', datetime(2)
  write (*, "(4a)") 'Info: ALCON ended   at date:', datetime(3), 'time:', datetime(4)
endif

end program alcon

