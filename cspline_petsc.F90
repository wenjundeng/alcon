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


! cubic spline interpolation using PETSc parallelization
module cspline_petsc
implicit none

#include "finclude/petscdef.h"

PetscInt, parameter :: cspline_ndata_default = 1
PetscReal, parameter :: cspline_testkindPetscReal = 0d0 ! a real constant to test the kind # of PetscReal
PetscInt, parameter :: cskpr = kind(cspline_testkindPetscReal) ! kind # of PetscReal

! when the requested x for interpolation is out of range of the tabulated x,
! extrapolation will be performed, this variable sets whether to show a warning
! this variable is public, set its value from the main program
! cspline_warning_extrapolation = 0: not show the warning
! cspline_warning_extrapolation = 1: only show the warning once, then this variable will be set to 0
! cspline_warning_extrapolation = 2: always show the warning
PetscInt :: cspline_warning_extrapolation = 2

! the following variables are private to this module and solely used for sharing data between the subroutines in this module.
! DO NOT access them directly. use the subroutines provided by this module.
PetscInt :: cspline_ndata = cspline_ndata_default
Mat :: cspline_matA
Vec :: cspline_vecc, cspline_vecc_seq, cspline_vecb
IS :: cspline_isfrom, cspline_isto
VecScatter :: cspline_gather
KSP :: cspline_ksp

save

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize PETSc objects for spline construction !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cspline_init(ndata)
implicit none
#include "finclude/petsc.h90"

PetscInt, intent(in) :: ndata ! # of data points

PetscErrorCode :: ierr

if (ndata < 3) then
  call PetscPrintf(PETSC_COMM_WORLD, "[cspline_init] Warning: # of data points (ndata) < 3.\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(PETSC_COMM_WORLD, "cspline is not initialized.\n", ierr)
  CHKERRQ(ierr)
  return
endif

cspline_ndata = ndata

! get MPI rank of current process
!call MPI_Comm_rank(PETSC_COMM_WORLD, cspline_mype, ierr)
!CHKERRQ(ierr)

! initialize cspline_matA
call MatCreate(PETSC_COMM_WORLD, cspline_matA, ierr)
CHKERRQ(ierr)
call MatSetType(cspline_matA, MATMPIAIJ, ierr)
CHKERRQ(ierr)
call MatSetSizes(cspline_matA, PETSC_DECIDE, PETSC_DECIDE, ndata, ndata, ierr)
CHKERRQ(ierr)
call MatMPIAIJSetPreallocation(cspline_matA, 3, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, ierr)
CHKERRQ(ierr)
  
! initialize cspline_vecc and cspline_vecb
call VecCreate(PETSC_COMM_WORLD, cspline_vecc, ierr)
CHKERRQ(ierr)
call VecSetSizes(cspline_vecc, PETSC_DECIDE, ndata, ierr)
CHKERRQ(ierr)
call VecSetFromOptions(cspline_vecc, ierr)
CHKERRQ(ierr)
call VecDuplicate(cspline_vecc, cspline_vecb, ierr)
CHKERRQ(ierr)

! initialize cspline_vecc_seq
call VecCreate(PETSC_COMM_SELF, cspline_vecc_seq, ierr)
CHKERRQ(ierr)
call VecSetSizes(cspline_vecc_seq, ndata, ndata, ierr)
CHKERRQ(ierr)
call VecSetType(cspline_vecc_seq, VECSEQ, ierr)
CHKERRQ(ierr)
call VecSetFromOptions(cspline_vecc_seq, ierr)
CHKERRQ(ierr)

call ISCreateStride(PETSC_COMM_WORLD, ndata, 0, 1, cspline_isfrom, ierr)
CHKERRQ(ierr)
call ISCreateStride(PETSC_COMM_SELF, ndata, 0, 1, cspline_isto, ierr)
CHKERRQ(ierr)

call VecScatterCreate(cspline_vecc, cspline_isfrom, cspline_vecc_seq, cspline_isto, cspline_gather, ierr)
CHKERRQ(ierr)

! initialize linear solver
call KSPCreate(PETSC_COMM_WORLD, cspline_ksp, ierr)
CHKERRQ(ierr)

end subroutine cspline_init


!!!!!!!!!!!!!!!!!!!!!!!!!!
! construct cubic spline !
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cspline_construct(x, y, d2y, bc, dyb1, dyb2)
implicit none
#include "finclude/petsc.h90"

! input tabulated x and y (y is a function of x)
PetscReal, dimension(0 : cspline_ndata - 1), intent(in) :: x
PetscScalar, dimension(0 : cspline_ndata - 1), intent(in) :: y
! output tabulated second derivative of y w.r.t. x
PetscScalar, dimension(0 : cspline_ndata - 1), intent(out) :: d2y
! boundary condition
! bc = 0: natural (free) boundary condition (second derivative d2y = 0)
! bc = 1: clamped boundary condition (first derivative dy1 and dyn specified)
PetscInt, intent(in) :: bc
! first derivative values for clamped boundary condition
PetscScalar, intent(in) :: dyb1, dyb2

PetscInt :: ilow, ihigh, k, nindex
PetscInt, dimension(0 : cspline_ndata - 1) :: indices
PetscScalar, dimension(0 : cspline_ndata - 1) :: values

PetscScalar, pointer :: vecc(:)

PetscErrorCode :: ierr

! check initialization
if (cspline_ndata == cspline_ndata_default) then
  call PetscPrintf(PETSC_COMM_WORLD, "[cspline_construct] Warning: not initialized.\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(PETSC_COMM_WORLD, "Call cspline_init first before calling cspline_construct.\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(PETSC_COMM_WORLD, "Returning 0 for d2y.\n", ierr)
  CHKERRQ(ierr)
  d2y = 0.0_cskpr
  return
endif

! check if x is strictly increasing
do k = 0, cspline_ndata - 2
  if (x(k) >= x(k + 1)) then
    call PetscPrintf(PETSC_COMM_WORLD, "[cspline_construct] Warning: tabulated x is not strictly increasing.\n", ierr)
    CHKERRQ(ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "Returning 0 for d2y.\n", ierr)
    CHKERRQ(ierr)
    d2y = 0.0_cskpr
    return
  endif
enddo

! assign values for cspline_vecb
call VecGetOwnershipRange(cspline_vecb, ilow, ihigh, ierr)
CHKERRQ(ierr)
nindex = 0
do k = ilow, ihigh - 1
  indices(nindex) = k
  if (k == 0) then
    if (bc == 0) then ! natural (free) boundary condition
      values(nindex) = 0.0_cskpr
    else ! clamped boundary condition
      values(nindex) = 3.0_cskpr * ((y(k + 1) - y(k)) / (x(k + 1) - x(k)) - dyb1)
    endif
  elseif (k == cspline_ndata - 1) then
    if (bc == 0) then ! natural (free) boundary condition
      values(nindex) = 0.0_cskpr
    else ! clamped boundary condition
      values(nindex) = 3.0_cskpr * (dyb2 - (y(k) - y(k - 1)) / (x(k) - x(k - 1)))
    endif
  else
    values(nindex) = 3.0_cskpr * ((y(k + 1) - y(k)) / (x(k + 1) - x(k)) - (y(k) - y(k - 1)) / (x(k) - x(k - 1)))
  endif
  nindex = nindex + 1
enddo
call VecSetValues(cspline_vecb, nindex, indices, values, INSERT_VALUES, ierr)
CHKERRQ(ierr)
call VecAssemblyBegin(cspline_vecb, ierr)
CHKERRQ(ierr)

call MatGetOwnershipRange(cspline_matA, ilow, ihigh, ierr)
CHKERRQ(ierr)

! assign values for cspline_matA
call VecGetOwnershipRange(cspline_vecb, ilow, ihigh, ierr)
CHKERRQ(ierr)
do k = ilow, ihigh - 1
  nindex = 0
  if (k == 0) then
    if (bc == 0) then ! natural (free) boundary condition
      indices(nindex) = k
      values(nindex) = 1.0_cskpr
      nindex = nindex + 1

      indices(nindex) = k + 1
      values(nindex) = 0.0_cskpr
      nindex = nindex + 1
    else ! clamped boundary condition
      indices(nindex) = k
      values(nindex) = 2.0_cskpr * (x(k + 1) - x(k))
      nindex = nindex + 1

      indices(nindex) = k + 1
      values(nindex) = x(k + 1) - x(k)
      nindex = nindex + 1
    endif
  elseif (k == cspline_ndata - 1) then
    if (bc == 0) then ! natural (free) boundary condition
      indices(nindex) = k
      values(nindex) = 1.0_cskpr
      nindex = nindex + 1

      indices(nindex) = k - 1
      values(nindex) = 0.0_cskpr
      nindex = nindex + 1
    else ! clamped boundary condition
      indices(nindex) = k
      values(nindex) = 2.0_cskpr * (x(k) - x(k - 1))
      nindex = nindex + 1

      indices(nindex) = k - 1
      values(nindex) = x(k) - x(k - 1)
      nindex = nindex + 1
    endif
  else ! (k /= 0 .and. k /= cspline_ndata - 1)
    indices(nindex) = k - 1
    values(nindex) = x(k) - x(k - 1)
    nindex = nindex + 1

    indices(nindex) = k
    values(nindex) = 2.0_cskpr * (x(k + 1) - x(k - 1))
    nindex = nindex + 1

    indices(nindex) = k + 1
    values(nindex) = x(k + 1) - x(k)
    nindex = nindex + 1
  endif
  call MatSetValues(cspline_matA, 1, k, nindex, indices, values, INSERT_VALUES, ierr)
  CHKERRQ(ierr)
enddo
call MatAssemblyBegin(cspline_matA, MAT_FINAL_ASSEMBLY, ierr)
CHKERRQ(ierr)

! finish assembling vector and matrix
call VecAssemblyEnd(cspline_vecb, ierr)
CHKERRQ(ierr)
call MatAssemblyEnd(cspline_matA, MAT_FINAL_ASSEMBLY, ierr)
CHKERRQ(ierr)

! prepare solver
call KSPSetOperators(cspline_ksp, cspline_matA, cspline_matA, DIFFERENT_NONZERO_PATTERN, ierr)
CHKERRQ(ierr)
! solve for cspline_vecc
call KSPSolve(cspline_ksp, cspline_vecb, cspline_vecc, ierr)
CHKERRQ(ierr)

call VecScatterBegin(cspline_gather, cspline_vecc, cspline_vecc_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)
CHKERRQ(ierr)
call VecScatterEnd(cspline_gather, cspline_vecc, cspline_vecc_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)
CHKERRQ(ierr)

call VecGetArrayF90(cspline_vecc_seq, vecc, ierr)
CHKERRQ(ierr)
d2y = 2.0_cskpr * vecc
call VecRestoreArrayF90(cspline_vecc_seq, vecc, ierr)
CHKERRQ(ierr)

end subroutine cspline_construct


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! construct cubic spline (wrapper for real spline)                !
! when PetscScalar is complex but real spline is wanted, use this !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cspline_construct_real(x, y, d2y, bc, dyb1, dyb2)
implicit none
#include "finclude/petsc.h90"

! input tabulated x and y (y is a function of x)
PetscReal, dimension(0 : cspline_ndata - 1), intent(in) :: x, y
! output tabulated second derivative of y w.r.t. x
PetscReal, dimension(0 : cspline_ndata - 1), intent(out) :: d2y
! boundary condition
! bc = 0: natural (free) boundary condition (second derivative d2y = 0)
! bc = 1: clamped boundary condition (first derivative dy1 and dyn specified)
PetscInt, intent(in) :: bc
! first derivative values for clamped boundary condition
PetscReal, intent(in) :: dyb1, dyb2

PetscScalar, dimension(0 : cspline_ndata - 1) :: d2y_cmplx

call cspline_construct(x, cmplx(y, 0.0_cskpr, cskpr), d2y_cmplx, bc, cmplx(dyb1, 0.0_cskpr, cskpr), cmplx(dyb2, 0.0_cskpr, cskpr))
d2y = real(d2y_cmplx, cskpr)

end subroutine cspline_construct_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolate for random argument !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cspline_interpolate(x, y, d2y, xip, yip)
implicit none
#include "finclude/petsc.h90"

! input tabulated x, y, and 2nd derivatives d2y
! d2y is calculated from cspline_construct
PetscReal, dimension(0 : cspline_ndata - 1), intent(in) :: x
PetscScalar, dimension(0 : cspline_ndata - 1), intent(in) :: y, d2y
! input x value to interpolate
PetscReal, intent(in) :: xip
! output interpolated y value
PetscScalar, intent(out) :: yip

PetscInt :: k1, k2, k3 ! indexes for bisection method

PetscErrorCode :: ierr

if (xip < x(0) .or. xip > x(cspline_ndata - 1)) then
  if (cspline_warning_extrapolation > 0) then
    call PetscPrintf(PETSC_COMM_WORLD, "[cspline_interpolate] Warning: xip is out of range of tabulated x.\n", ierr)
    CHKERRQ(ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "Performing extrapolation instead.\n", ierr)
    CHKERRQ(ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "If extrapolation is wanted, you can set cspline_warning_", ierr)
    CHKERRQ(ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "extrapolation = 0 in the main program.\n", ierr)
    CHKERRQ(ierr)
    if (cspline_warning_extrapolation == 1) cspline_warning_extrapolation = 0
  endif
  if (xip < x(0)) then
    k1 = 0
    k2 = 1
  else
    k1 = cspline_ndata - 2
    k2 = cspline_ndata - 1
  endif
else
  k1 = 0
  k2 = cspline_ndata - 1
  do while (k2 - k1 > 1)
    k3 = (k1 + k2) / 2
    if (xip < x(k3)) then
      k2 = k3
    else
      k1 = k3
    endif
  enddo
endif
yip = y(k1) &
  + (xip - x(k1)) * ((y(k2) - y(k1)) / (x(k2) - x(k1)) - (x(k2) - x(k1)) * (2.0_cskpr * d2y(k1) + d2y(k2)) / 6.0_cskpr) &
  + (xip - x(k1))**2 * d2y(k1) / 2.0_cskpr &
  + (xip - x(k1))**3 * (d2y(k2) - d2y(k1)) / (6.0_cskpr * (x(k2) - x(k1)))

end subroutine cspline_interpolate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolate for random argument (wrapper for real spline)       !
! when PetscScalar is complex but real spline is wanted, use this !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cspline_interpolate_real(x, y, d2y, xip, yip)
implicit none
#include "finclude/petsc.h90"

! input tabulated x, y, and 2nd derivatives d2y
! d2y is calculated from cspline_construct
PetscReal, dimension(0 : cspline_ndata - 1), intent(in) :: x, y, d2y
! input x value to interpolate
PetscReal, intent(in) :: xip
! output interpolated y value
PetscReal, intent(out) :: yip

PetscScalar :: yip_cmplx

call cspline_interpolate(x, cmplx(y, 0.0_cskpr, cskpr), cmplx(d2y, 0.0_cskpr, cskpr), xip, yip_cmplx)

yip = real(yip_cmplx, cskpr)

end subroutine cspline_interpolate_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolate for increasing arguments !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cspline_interpolate_inc(x, y, d2y, xip, yip)
implicit none
#include "finclude/petsc.h90"

! input tabulated x, y, and 2nd derivatives d2y
! d2y is calculated from cspline_construct
PetscReal, dimension(0 : cspline_ndata - 1), intent(in) :: x
PetscScalar, dimension(0 : cspline_ndata - 1), intent(in) :: y, d2y
! input x values to interpolate
PetscReal, dimension(:), intent(in) :: xip
! output interpolated y values
PetscScalar, dimension(:), intent(out) :: yip

PetscInt :: k1, k2, k3 ! indexes for bisection method
PetscInt :: ix1, ix2 ! index boundaries for xip
PetscInt :: i

PetscErrorCode :: ierr

ix1 = lbound(xip, 1)
ix2 = ubound(xip, 1)

if (ix1 /= lbound(yip, 1) .or. ix2 /= ubound(yip, 1)) then
  call PetscPrintf(PETSC_COMM_WORLD, "[cspline_interpolate_inc] Warning: dimensions of xip and yip do not match.\n", ierr)
  CHKERRQ(ierr)
  call PetscPrintf(PETSC_COMM_WORLD, "Returning 0 for yip.\n", ierr)
  CHKERRQ(ierr)
#if defined(PETSC_USE_COMPLEX)
  yip = (0.0_cskpr, 0.0_cskpr)
#else
  yip = 0.0_cskpr
#endif
  return
endif

! check if xip is strictly increasing
do i = ix1, ix2 - 1
  if (xip(i) >= xip(i + 1)) then
    call PetscPrintf(PETSC_COMM_WORLD, "[cspline_interpolate_inc] Warning: xip is not strictly increasing.\n", ierr)
    CHKERRQ(ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "Returning 0 for yip.\n", ierr)
    CHKERRQ(ierr)
#if defined(PETSC_USE_COMPLEX)
  yip = (0.0_cskpr, 0.0_cskpr)
#else
  yip = 0.0_cskpr
#endif
    return
  endif
enddo

if (xip(ix1) < x(0) .or. xip(ix2) > x(cspline_ndata - 1)) then
  if (cspline_warning_extrapolation > 0) then
    call PetscPrintf(PETSC_COMM_WORLD, "[cspline_interpolate] Warning: xip is out of range of tabulated x.\n", ierr)
    CHKERRQ(ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "Performing extrapolation instead.\n", ierr)
    CHKERRQ(ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "If extrapolation is wanted, you can set cspline_warning_", ierr)
    CHKERRQ(ierr)
    call PetscPrintf(PETSC_COMM_WORLD, "extrapolation = 0 in the main program\n", ierr)
    CHKERRQ(ierr)
    if (cspline_warning_extrapolation == 1) cspline_warning_extrapolation = 0
  endif
endif

k1 = 0
k2 = cspline_ndata - 1
do i = ix1, ix2
  if (i == ix1 .or. (xip(i) >= x(k2) .and. k2 < cspline_ndata - 1)) then
    if (i /= ix1 .and. xip(i) >= x(k2)) then
      k1 = k2
      k2 = cspline_ndata - 1
    endif
    do while (k2 - k1 > 1)
      k3 = (k1 + k2) / 2
      if (xip(i) < x(k3)) then
        k2 = k3
      else
        k1 = k3
      endif
    enddo
  endif
  yip(i) = y(k1) &
    + (xip(i) - x(k1)) * ((y(k2) - y(k1)) / (x(k2) - x(k1)) - (x(k2) - x(k1)) * (2.0_cskpr * d2y(k1) + d2y(k2)) / 6.0_cskpr) &
    + (xip(i) - x(k1))**2 * d2y(k1) / 2.0_cskpr &
    + (xip(i) - x(k1))**3 * (d2y(k2) - d2y(k1)) / (6.0_cskpr * (x(k2) - x(k1)))
enddo

end subroutine cspline_interpolate_inc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! interpolate for increasing arguments (wrapper for real spline)  !
! when PetscScalar is complex but real spline is wanted, use this !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cspline_interpolate_inc_real(x, y, d2y, xip, yip)
implicit none
#include "finclude/petsc.h90"

! input tabulated x, y, and 2nd derivatives d2y
! d2y is calculated from cspline_construct
PetscReal, dimension(0 : cspline_ndata - 1), intent(in) :: x, y, d2y
! input x values to interpolate
PetscReal, dimension(:), intent(in) :: xip
! output interpolated y values
PetscReal, dimension(:), intent(out) :: yip

PetscScalar, dimension(:), allocatable :: yip_cmplx

allocate (yip_cmplx(lbound(yip, 1) : ubound(yip, 1)))

call cspline_interpolate_inc(x, cmplx(y, 0.0_cskpr, cskpr), cmplx(d2y, 0.0_cskpr, cskpr), xip, yip_cmplx)

yip = real(yip_cmplx, cskpr)

deallocate (yip_cmplx)

end subroutine cspline_interpolate_inc_real


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! finalization, destroy PETSc objects !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cspline_final
implicit none
#include "finclude/petsc.h90"

PetscErrorCode :: ierr

call KSPDestroy(cspline_ksp, ierr)
CHKERRQ(ierr)
call ISDestroy(cspline_isfrom, ierr)
CHKERRQ(ierr)
call ISDestroy(cspline_isto, ierr)
CHKERRQ(ierr)
call VecScatterDestroy(cspline_gather, ierr)
CHKERRQ(ierr)
call VecDestroy(cspline_vecc, ierr)
CHKERRQ(ierr)
call VecDestroy(cspline_vecc_seq, ierr)
CHKERRQ(ierr)
call VecDestroy(cspline_vecb, ierr)
CHKERRQ(ierr)
call MatDestroy(cspline_matA, ierr)
CHKERRQ(ierr)

cspline_ndata = cspline_ndata_default

end subroutine cspline_final

end module cspline_petsc

