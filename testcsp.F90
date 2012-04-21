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


! test cspline_petsc
program testcsp
use cspline_petsc
implicit none
#include "finclude/petsc.h90"

PetscInt, parameter :: nx = 7, nxip = 101

PetscInt, parameter :: fexact = 231, ftab = 232, fip = 233

PetscReal, dimension(nx) :: x
!PetscScalar, dimension(nx) :: y, d2y
PetscReal, dimension(nx) :: y, d2y

PetscReal, dimension(nxip) :: xip
!PetscScalar, dimension(nxip) :: yip
PetscReal, dimension(nxip) :: yip

PetscInt :: i

PetscInt :: mype
PetscErrorCode :: ierr
character(400) :: msg

call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
CHKERRQ(ierr)
call MPI_Comm_rank(PETSC_COMM_WORLD, mype, ierr)
CHKERRQ(ierr)

x = (/ ((i - ((nx + 1) / 2.0_cskpr)) / ((nx - 1) / 2.0_cskpr), i = 1, nx) /)
write (msg, *) "x = ", (x(i), i = 1, nx), "\n"
call PetscPrintf(PETSC_COMM_WORLD, msg, ierr)
CHKERRQ(ierr)

y = x * sin(2.0_cskpr * PETSC_PI * x + 1.0_cskpr) ! this is the function to test

write (msg, *) "y = ", (y(i), i = 1, nx), "\n"
call PetscPrintf(PETSC_COMM_WORLD, msg, ierr)
CHKERRQ(ierr)
write (msg, *) "Re(y) = ", (real(y(i), cskpr), i = 1, nx), "\n"
call PetscPrintf(PETSC_COMM_WORLD, msg, ierr)
CHKERRQ(ierr)

call cspline_init(nx)
!call cspline_construct(x, y, d2y, 0, (0.0_cskpr, 0.0_cskpr), (0.0_cskpr, 0.0_cskpr))
!call cspline_construct(x, y, d2y, 1, (-2.533_cskpr, 0.0_cskpr), (4.2363_cskpr, 0.0_cskpr))
call cspline_construct_real(x, y, d2y, 0, 0.0_cskpr, 0.0_cskpr)
!call cspline_construct_real(x, y, d2y, 1, -2.533_cskpr, 4.2363_cskpr)

write (msg, *) "d2y = ", (d2y(i), i = 1, nx), "\n"
call PetscPrintf(PETSC_COMM_WORLD, msg, ierr)
CHKERRQ(ierr)
!write (msg, *) "Re(d2y) = ", (real(d2y(i), cskpr), i = 1, nx), "\n"
!call PetscPrintf(PETSC_COMM_WORLD, msg, ierr)
!CHKERRQ(ierr)

xip = (/ (2.0_cskpr * (i - ((nxip + 1) / 2.0_cskpr)) / ((nxip - 1) / 2.0_cskpr), i = 1, nxip) /)
write (msg, *) "xip = ", xip(1), xip(2), xip(3), "...", xip(nxip - 2), xip(nxip - 1), xip(nxip), "\n"
call PetscPrintf(PETSC_COMM_WORLD, msg, ierr)
CHKERRQ(ierr)

!call cspline_interpolate(x, y, d2y, 1.0_cskpr, yip(1))
!write (msg, *) "yip(1) = ", yip(1), "\n"
!call PetscPrintf(PETSC_COMM_WORLD, msg, ierr)
!CHKERRQ(ierr)

!call cspline_interpolate_inc(x, y, d2y, xip, yip)
call cspline_interpolate_inc_real(x, y, d2y, xip, yip)
write (msg, *) "yip = ", yip(1), yip(2), yip(3), "...", yip(nxip - 2), yip(nxip - 1), yip(nxip), "\n"
call PetscPrintf(PETSC_COMM_WORLD, msg, ierr)
CHKERRQ(ierr)

if (mype == 0) then
  open (fexact, file = 'fexact.dat', status = 'replace') ! exact function
  open (ftab, file = 'ftab.dat', status = 'replace') ! tabulated points
  open (fip, file = 'fip.dat', status = 'replace') ! interpolated function

  do i = 1, nxip
    write (fexact, *) xip(i), xip(i) * sin(2.0_cskpr * PETSC_PI * xip(i) + 1.0_cskpr)
    write (fip, *) xip(i), real(yip(i), cskpr)
  enddo
  do i = 1, nx
    write (ftab, *) x(i), real(y(i), cskpr)
  enddo

  close (fexact)
  close (ftab)
  close (fip)
endif

call cspline_final

call PetscFinalize(ierr)
CHKERRQ(ierr)

end program testcsp

