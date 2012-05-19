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

module alcon_solver
use alcon_input
implicit none
#include "finclude/petscdef.h"

! these variables are private to this module and solely used for sharing data between the subroutines in this module
! DO NOT access them directly. use the subroutines provided by this module.
MPI_Comm :: acs_comm
Mat :: acs_matA, acs_matB ! matrixes described in Eqs. (A.20) and (A.21) in [Nuclear Fusion 52, 043006 (2012)]
Vec :: acs_vecEigenRe, acs_vecEigenIm ! eigenvectores. note that acs_vecEigenIm is only for conforming SLEPc syntax. both real and imaginary parts of the eigenvector are stored in acs_vecEigenRe
PetscInt :: acs_npe, acs_mype ! # of MPI processes; rank of current MPI process

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solver initialization, create PETSc and SLEPc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alcon_solver_init(comm)
implicit none
#include "finclude/petsc.h90"

MPI_Comm, intent(in) :: comm
PetscErrorCode :: ierr

acs_comm = comm

! initialize acs_matA
call MatCreate(acs_comm, acs_matA, ierr)
CHKERRQ(ierr)
call MatSetType(acs_matA, MATMPIAIJ, ierr)
CHKERRQ(ierr)
if (finitebeta > 2) then
  call MatSetSizes(acs_matA, PETSC_DECIDE, PETSC_DECIDE, (m2 - m1 + 1) * 2, (m2 - m1 + 1) * 2, ierr)
  CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(acs_matA, noffdiag * 2 + 1, PETSC_NULL_INTEGER, noffdiag * 3 + 1, PETSC_NULL_INTEGER, ierr)
  CHKERRQ(ierr)
else
  call MatSetSizes(acs_matA, PETSC_DECIDE, PETSC_DECIDE, m2 - m1 + 1, m2 - m1 + 1, ierr)
  CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(acs_matA, noffdiag * 2 + 1, PETSC_NULL_INTEGER, noffdiag, PETSC_NULL_INTEGER, ierr)
  CHKERRQ(ierr)
endif

! initialize acs_matB (duplicate from acs_matA)
! MatDuplicate does not work for unassembled matrix
!call MatDuplicate(acs_matA, MAT_DO_NOT_COPY_VALUES, acs_matB, ierr)
!CHKERRQ(ierr)

! initialize acs_matB
call MatCreate(acs_comm, acs_matB, ierr)
CHKERRQ(ierr)
call MatSetType(acs_matB, MATMPIAIJ, ierr)
CHKERRQ(ierr)
if (finitebeta > 2) then
  call MatSetSizes(acs_matB, PETSC_DECIDE, PETSC_DECIDE, (m2 - m1 + 1) * 2, (m2 - m1 + 1) * 2, ierr)
  CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(acs_matB, noffdiag * 2 + 1, PETSC_NULL_INTEGER, noffdiag * 3 + 1, PETSC_NULL_INTEGER, ierr)
  CHKERRQ(ierr)
else
  call MatSetSizes(acs_matB, PETSC_DECIDE, PETSC_DECIDE, m2 - m1 + 1, m2 - m1 + 1, ierr)
  CHKERRQ(ierr)
  call MatMPIAIJSetPreallocation(acs_matB, noffdiag * 2 + 1, PETSC_NULL_INTEGER, noffdiag, PETSC_NULL_INTEGER, ierr)
  CHKERRQ(ierr)
endif

! initialize acs_vecEigenRe and acs_vecEigenIm
call VecCreate(acs_comm, acs_vecEigenRe, ierr)
CHKERRQ(ierr)
if (finitebeta > 2) then
  call VecSetSizes(acs_vecEigenRe, PETSC_DECIDE, (m2 - m1 + 1) * 2, ierr)
  CHKERRQ(ierr)
else
  call VecSetSizes(acs_vecEigenRe, PETSC_DECIDE, m2 - m1 + 1, ierr)
  CHKERRQ(ierr)
endif
call VecSetFromOptions(acs_vecEigenRe, ierr)
CHKERRQ(ierr)
call VecDuplicate(acs_vecEigenRe, acs_vecEigenIm, ierr)
CHKERRQ(ierr)

! epsAlcon needs to be created and destroyed for every time of solving, so it is initialized in the main program

! get # of MPI processes
call MPI_Comm_size(acs_comm, acs_npe, ierr)
CHKERRQ(ierr)
! get rank of current MPI process
call MPI_Comm_rank(acs_comm, acs_mype, ierr)
CHKERRQ(ierr)

end subroutine alcon_solver_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solver finalization, destroy PETSc and SLEPc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alcon_solver_final
implicit none
#include "finclude/petsc.h90"

PetscInt :: idatout

PetscErrorCode :: ierr

! destroy PETSc and SLEPc objects
call VecDestroy(acs_vecEigenIm, ierr)
CHKERRQ(ierr)
call VecDestroy(acs_vecEigenRe, ierr)
CHKERRQ(ierr)
!call EPSDestroy(epsAlcon, ierr)
!CHKERRQ(ierr)
call MatDestroy(acs_matB, ierr)
CHKERRQ(ierr)
call MatDestroy(acs_matA, ierr)
CHKERRQ(ierr)

end subroutine alcon_solver_final


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare matrix values for the solver !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alcon_solver_prepare(irad)
use alcon_eqdata
implicit none
#include "finclude/petsc.h90"

PetscInt, intent(in) :: irad ! prepare solver for which radial grid

PetscInt :: ilow, ihigh, nindex, i, j
PetscInt, dimension(0 : (m2 - m1 + 1) * 2 - 1) :: indices
PetscScalar, dimension(0 : (m2 - m1 + 1) * 2 - 1) :: values

PetscErrorCode :: ierr
character(400) :: msg

! set values for acs_matA and acs_matB
call MatGetOwnershipRange(acs_matA, ilow, ihigh, ierr)
CHKERRQ(ierr)
if (verbosity >= 2 .and. irad == 1) then
  ! flush() function is not a standard Fortran function
  !flush(6)
  !call MPI_Barrier(acs_comm, ierr)
  !CHKERRQ(ierr)
  write (msg, *) "[alcon_solver_prepare] Info: acs_mype = ", acs_mype, ", ilow = ", ilow, ", ihigh = ", ihigh
  write (*, "(a)") trim(adjustl(msg))
  !flush(6)
  !call MPI_Barrier(acs_comm, ierr)
  !CHKERRQ(ierr)
endif
do i = ilow, ihigh - 1
  if (finitebeta == 2) then
    ! beta * G^dagger * M *G part in acs_matA
    nindex = 0
    indices(nindex) = i
    values(nindex) = acdprofile(4, irad) * (real(n, kpr) * acdprofile(2, irad) - real(m1 + i, kpr))**2 / acdprofile(3, irad)
    nindex = nindex + 1
    call MatSetValues(acs_matA, 1, i, nindex, indices, values, INSERT_VALUES, ierr)

    nindex = 0
    do j = max(0, i - noffdiag), min(m2 - m1, i + noffdiag)
      ! L part
      indices(nindex) = j
      if (j < i) then
        values(nindex) = conjg(acdfft(i - j, 4, irad))
      else
        values(nindex) = acdfft(j - i, 4, irad)
      endif
      nindex = nindex + 1
    enddo
    call MatSetValues(acs_matB, 1, i, nindex, indices, values, INSERT_VALUES, ierr)
  else
    nindex = 0
    if (i < m2 - m1 + 1) then
      do j = max(0, i - noffdiag), min(m2 - m1, i + noffdiag)
        ! G^dagger * H * G part in acs_matA
        indices(nindex) = j
        if (j < i) then
          values(nindex) = conjg(acdfft(i - j, 1, irad))
        else
          values(nindex) = acdfft(j - i, 1, irad)
        endif
        values(nindex) = values(nindex) &
          * (real(n, kpr) * acdprofile(2, irad) - real(m1 + i, kpr)) &
          * (real(n, kpr) * acdprofile(2, irad) - real(m1 + j, kpr))

        ! for slow sound approximation, A = G^dagger * H * G part + N
        if (finitebeta == 1) then
          if (j < i) then
            values(nindex) = values(nindex) + conjg(acdfft(i - j, 5, irad))
          else
            values(nindex) = values(nindex) + acdfft(j - i, 5, irad)
          endif
        endif
        nindex = nindex + 1

        ! -beta# * K part in acs_matA
        if (finitebeta > 2) then
          indices(nindex) = m2 - m1 + 1 + j
          if (j < i) then
            values(nindex) = -acdprofile(4, irad) * conjg(acdfft(i - j, 3, irad))
          else
            values(nindex) = -acdprofile(4, irad) * acdfft(j - i, 3, irad)
          endif
          nindex = nindex + 1
        endif
      enddo
      call MatSetValues(acs_matA, 1, i, nindex, indices, values, INSERT_VALUES, ierr)
      CHKERRQ(ierr)

      ! J part in acs_matB
      nindex = 0
      do j = max(0, i - noffdiag), min(m2 - m1, i + noffdiag)
        indices(nindex) = j
        if (j < i) then
          values(nindex) = conjg(acdfft(i - j, 2, irad))
        else
          values(nindex) = acdfft(j - i, 2, irad)
        endif
        nindex = nindex + 1
      enddo
      call MatSetValues(acs_matB, 1, i, nindex, indices, values, INSERT_VALUES, ierr)
      CHKERRQ(ierr)
    else ! (i < m2 - m1 + 1)
      ! this means (i >= m2 - m1 + 1), which implies finitebeta > 2
      ! beta# * G^dagger * M * G part
      nindex = 0
      indices(nindex) = i
      values(nindex) = acdprofile(4, irad) &
        * (real(n, kpr) * acdprofile(2, irad) - real(m1 + (i - (m2 - m1 + 1)), kpr))**2 / acdprofile(3, irad)
      nindex = nindex + 1
      call MatSetValues(acs_matA, 1, i, nindex, indices, values, INSERT_VALUES, ierr)
      CHKERRQ(ierr)

      nindex = 0
      do j = max(m2 - m1 + 1, i - noffdiag), min((m2 - m1 + 1) * 2 - 1, i + noffdiag)
        ! K part in acs_matB
        indices(nindex) = j - (m2 - m1 + 1)
        if (j < i) then
          values(nindex) = conjg(acdfft(i - j, 3, irad))
        else
          values(nindex) = acdfft(j - i, 3, irad)
        endif
        nindex = nindex + 1

        ! L part in acs_matB
        indices(nindex) = j
        if (j < i) then
          values(nindex) = conjg(acdfft(i - j, 4, irad))
        else
          values(nindex) = acdfft(j - i, 4, irad)
        endif
        nindex = nindex + 1
      enddo
      call MatSetValues(acs_matB, 1, i, nindex, indices, values, INSERT_VALUES, ierr)
      CHKERRQ(ierr)
    endif ! else of (i < m2 - m1 + 1)
  endif ! else of (finitebeta == 2)
enddo ! i = ilow, ihigh
call MatAssemblyBegin(acs_matA, MAT_FINAL_ASSEMBLY, ierr)
CHKERRQ(ierr)
call MatAssemblyBegin(acs_matB, MAT_FINAL_ASSEMBLY, ierr)
CHKERRQ(ierr)
call MatAssemblyEnd(acs_matA, MAT_FINAL_ASSEMBLY, ierr)
CHKERRQ(ierr)
call MatAssemblyEnd(acs_matB, MAT_FINAL_ASSEMBLY, ierr)
CHKERRQ(ierr)

end subroutine alcon_solver_prepare


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve eigenvalue and eigenvector !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alcon_solver_solve(rad, rho_M)
use alcon_output
implicit none
#include "finclude/petsc.h90"
#include "finclude/slepc.h"

PetscReal, intent(in) :: rad ! radial coordinate, usually the noralized radial coordinate rho = sqrt(psitor / psitor_wall)
PetscReal, intent(in) :: rho_M ! mass density

EPS :: epsAlcon
PetscInt :: nconv, maxit, nev_out, its
PetscScalar :: evalre, evalim
PetscReal :: eivecabsmax, omegarescaled
EPSType :: tname
PetscReal :: tol

PetscInt :: iconv, ieivecabsmax

PetscErrorCode :: ierr
character(400) :: msg

! initialize eigenvalue solver epsAlcon
call EPSCreate(acs_comm, epsAlcon, ierr)
CHKERRQ(ierr)
!call EPSSetProblemType(epsAlcon, EPS_PGNHEP, ierr)
!CHKERRQ(ierr)
call EPSSetWhichEigenpairs(epsAlcon, EPS_SMALLEST_MAGNITUDE, ierr)
CHKERRQ(ierr)
if (finitebeta > 2) then
  call EPSSetDimensions(epsAlcon, m2 - m1 + 1 + nev, PETSC_DECIDE, PETSC_DECIDE, ierr)
  CHKERRQ(ierr)
elseif (finitebeta == 2) then
  call EPSSetDimensions(epsAlcon, m2 - m1 + 1, PETSC_DECIDE, PETSC_DECIDE, ierr)
  CHKERRQ(ierr)
else
  call EPSSetDimensions(epsAlcon, nev, PETSC_DECIDE, PETSC_DECIDE, ierr)
  CHKERRQ(ierr)
endif
call EPSSetOperators(epsAlcon, acs_matA, acs_matB, ierr)
CHKERRQ(ierr)

! near mode rational surfaces, the lowest continuum has near-0 frequency,
!   the eigenvalue problem becomes near-singular.  EPSLAPACK is found to handle
!   this situation quite well.  however, it does not scale well.  therefore, if
!   you know an EPSType that can handle near-singular eigenvalue problem well
!   and meanwhile scales well, please let me know <wdeng@wdeng.info>.
!   thank you.
call EPSSetType(epsAlcon, EPSLAPACK, ierr)
CHKERRQ(ierr)
!call EPSSetFromOptions(epsAlcon, ierr)
!CHKERRQ(ierr)

! solve the eigenvalue problem
call EPSSolve(epsAlcon, ierr)
CHKERRQ(ierr)

if (verbosity >= 3) then
  ! write out information about the solver on irad == 1 or irad == iraddiag
  call EPSGetIterationNumber(epsAlcon, its, ierr)
  CHKERRQ(ierr)
  write (msg, *) "[alcon_solver_solve] Info: # of iterations: ", its, "\n"
  call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
  CHKERRQ(ierr)
  call EPSGetType(epsAlcon, tname, ierr)
  CHKERRQ(ierr)
  write (msg, *) "[alcon_solver_solve] Info: solution method: ", tname, "\n"
  call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
  CHKERRQ(ierr)
  call EPSGetDimensions(epsAlcon, nev_out, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
  CHKERRQ(ierr)
  write (msg, *) "[alcon_solver_solve] Info: # of requested eigenvalues: ", nev_out, "\n"
  call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
  CHKERRQ(ierr)
  call EPSGetTolerances(epsAlcon, tol, maxit, ierr)
  CHKERRQ(ierr)
  write (msg, *) "[alcon_solver_solve] Info: stopping condition: tol = ", tol, ", maxit = ", maxit, "\n"
  call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
  CHKERRQ(ierr)
endif

call EPSGetConverged(epsAlcon, nconv, ierr)
CHKERRQ(ierr)
if (verbosity >= 2) then
  write (msg, *) "[alcon_solver_solve] Info: # of converged solutions: nconv = ", nconv, "\n"
  call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
  CHKERRQ(ierr)
endif
! analyze converged solutions
do iconv = 0, nconv - 1
  call EPSGetEigenpair(epsAlcon, iconv, evalre, evalim, acs_vecEigenRe, acs_vecEigenIm, ierr)
  CHKERRQ(ierr)
  if (verbosity >= 3 .and. iconv == 0) then
    call VecView(acs_vecEigenRe, PETSC_VIEWER_STDOUT_WORLD, ierr)
    CHKERRQ(ierr)
    call VecView(acs_vecEigenIm, PETSC_VIEWER_STDOUT_WORLD, ierr)
    CHKERRQ(ierr)
  endif
  ! the absolute value of omega's imaginary part has to be sufficiently small
  if (abs(aimag(sqrt(evalre))) < real(sqrt(evalre), kpr) * imagrealratiocutoff) then
    omegarescaled = real(sqrt(evalre / rho_M), kpr) * omegascale
    if (omegacutoff < 0.0_kpr .or. omegarescaled < omegacutoff) then
      call VecAbs(acs_vecEigenRe, ierr)
      CHKERRQ(ierr)
      call VecMax(acs_vecEigenRe, ieivecabsmax, eivecabsmax, ierr)
      CHKERRQ(ierr)
      ! put data to aco_data of alcon_output module
      aco_ndata(ieivecabsmax) = aco_ndata(ieivecabsmax) + 1
      if (aco_ndata(ieivecabsmax) > aco_data_stacksize) then
        call PetscPrintf(MPI_COMM_WORLD, "[alcon_solver_solve] Error: aco_data is full. Incr", ierr)
        CHKERRQ(ierr)
        call PetscPrintf(MPI_COMM_WORLD, "ease aco_data_stacksize and re-", ierr)
        CHKERRQ(ierr)
        call PetscPrintf(MPI_COMM_WORLD, "compile and run again.\n", ierr)
        CHKERRQ(ierr)
        call alcon_solver_final
        ! too complicated to inform all other processes and terminate nicely, so just MPI_Abort
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        CHKERRQ(ierr)
      endif
      aco_data(ieivecabsmax, 1, aco_ndata(ieivecabsmax)) = rad
      aco_data(ieivecabsmax, 2, aco_ndata(ieivecabsmax)) = omegarescaled

      if (verbosity >= 2 .and. iconv == 0) then
        write (msg, *) "[alcon_solver] Info: lowest omega: ", omegarescaled, "\n"
        call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
        CHKERRQ(ierr)
        write (msg, *) "[alcon_solver] Info: corresponding imaginary part: ", aimag(sqrt(evalre / rho_M)) * omegascale, "\n"
        call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
        CHKERRQ(ierr)
        if (ieivecabsmax < m2 - m1 + 1) then
          write (msg, *) "[alcon_solver] Info: corresponding m: ", m1 + ieivecabsmax, "\n"
        else
          write (msg, *) "[alcon_solver] Info: corresponding m: ", m1 + ieivecabsmax - (m2 - m1 + 1), "\n"
        endif
        call PetscPrintf(MPI_COMM_WORLD, trim(adjustl(msg)), ierr)
        CHKERRQ(ierr)
      endif
    endif ! (omegacutoff < 0.0_kpr .or. omegarescaled < omegacutoff)
  endif ! (abs(aimag(sqrt(evalre))) < real(sqrt(evalre), kpr) * imagrealratiocutoff)
enddo ! do i = 0, nconv - 1
! destroy epsAlcon
call EPSDestroy(epsAlcon, ierr)
CHKERRQ(ierr)

end subroutine alcon_solver_solve

end module alcon_solver

