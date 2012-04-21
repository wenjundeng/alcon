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

module alcon_output
use alcon_input
implicit none
#include "finclude/petscdef.h"

PetscInt, parameter :: aco_data_stacksize = nrad * 10

! output data
PetscReal, dimension(0 : (m2 - m1 + 1) * 2 - 1, 2, aco_data_stacksize) :: aco_data
! # of output data for each m-harmonic
PetscInt, dimension(0 : (m2 - m1 + 1) * 2 - 1) :: aco_ndata = 0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! returns output data file name for given m-harmonic index !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(40) function alcon_output_filename(im, withpath, withext)
implicit none

! index in the 0 : (m2 - m1 + 1) * 2 - 1 m-harmonic space
PetscInt, intent(in) :: im
! whether to include path
PetscBool, intent(in) :: withpath
! whether to include extention (".dat")
PetscBool, intent(in) :: withext

character(10) :: path, ext

if (withpath) then
  path = "output/"
else
  path = ""
endif

if (withext) then
  ext = ".dat"
else
  ext = ""
endif

if (finitebeta > 2) then
  if (im < m2 - m1 + 1) then
    if (m1 + im < 0) then
      write (alcon_output_filename, '(a, "a_n", i3.3, "m_", i2.2, a)') trim(path), n, abs(m1 + im), trim(ext)
    else
      write (alcon_output_filename, '(a, "a_n", i3.3, "m", i3.3, a)') trim(path), n, m1 + im, trim(ext)
    endif
  else
    if (m1 + im - (m2 - m1 + 1) < 0) then
      write (alcon_output_filename, '(a, "s_n", i3.3, "m_", i2.2, a)') trim(path), n, abs(m1 + im - (m2 - m1 + 1)), trim(ext)
    else
      write (alcon_output_filename, '(a, "s_n", i3.3, "m", i3.3, a)') trim(path), n, m1 + im - (m2 - m1 + 1), trim(ext)
    endif
  endif
else
  if (m1 + im < 0) then
    write (alcon_output_filename, '(a, "n", i3.3, "m_", i2.2, a)') trim(path), n, abs(m1 + im), trim(ext)
  else
    write (alcon_output_filename, '(a, "n", i3.3, "m", i3.3, a)') trim(path), n, m1 + im, trim(ext)
  endif
endif
end function alcon_output_filename


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out Alfven continuum data                       !
! only comm_radhead members should call this subroutine !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alcon_output_data(comm_radhead, ndom_rad)
implicit none
#include "finclude/petsc.h90"

! communicator for radial domain heads, who are responsible for writing output data
MPI_Comm, intent(in) :: comm_radhead
! # of radial domains
PetscInt, intent(in) :: ndom_rad

PetscInt, parameter :: tag = 1
PetscInt, parameter :: nbyteperei = 53 ! # of bytes per eigenvalue output

PetscInt, dimension(0 : (m2 - m1 + 1) * 2 - 1) :: aco_ndata_total, aco_ndata_lower, intbuf
PetscInt :: im, mype, nm, fdat, idata
PetscInt, dimension(MPI_STATUS_SIZE) :: st
integer(kind = MPI_OFFSET_KIND) :: disp
PetscErrorCode :: ierr

character(nbyteperei) :: charbuf
character(400) :: msg

if (verbosity >= 2) then
  call PetscPrintf(comm_radhead, "[alcon_output_data] Info: preparing data about file structures...\n", ierr)
  CHKERRQ(ierr)
endif
call MPI_Comm_rank(comm_radhead, mype, ierr)
CHKERRQ(ierr)

if (finitebeta > 2) then
  nm = (m2 - m1 + 1) * 2
else
  nm = m2 - m1 + 1
endif

call MPI_Allreduce(aco_ndata, aco_ndata_total, nm, MPIU_INTEGER, MPI_SUM, comm_radhead, ierr)
CHKERRQ(ierr)

if (mype > 0) then
  call MPI_Recv(intbuf, nm, MPIU_INTEGER, mype - 1, tag, comm_radhead, st, ierr)
  CHKERRQ(ierr)
  aco_ndata_lower = intbuf
else
  aco_ndata_lower = 0
endif
if (mype < ndom_rad - 1) then
  intbuf = aco_ndata_lower + aco_ndata
  call MPI_Send(intbuf, nm, MPIU_INTEGER, mype + 1, tag, comm_radhead, ierr)
  CHKERRQ(ierr)
endif

if (verbosity >= 2) then
  call PetscPrintf(comm_radhead, "[alcon_output_data] Info: writing to files...\n", ierr)
  CHKERRQ(ierr)
endif

do im = 0, nm - 1
  if (aco_ndata_total(im) == 0) cycle
  if (verbosity >= 2) then
    write (msg, "(3a)") "[alcon_output_data] Info: opening file ", trim(alcon_output_filename(im, PETSC_TRUE, PETSC_TRUE)), " ...\n"
    call PetscPrintf(comm_radhead, trim(msg), ierr)
    CHKERRQ(ierr)
  endif
  call MPI_File_open(comm_radhead, trim(alcon_output_filename(im, PETSC_TRUE, PETSC_TRUE)), &
    MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fdat, ierr)
  CHKERRQ(ierr)
  disp = aco_ndata_lower(im) * nbyteperei
  if (verbosity >= 2) then
    write (msg, "(3a)") "[alcon_output_data] Info: setting view for file ", &
      trim(alcon_output_filename(im, PETSC_TRUE, PETSC_TRUE)), " ...\n"
    call PetscPrintf(comm_radhead, trim(msg), ierr)
    CHKERRQ(ierr)
    write (msg, *) "[alcon_output_data] Info: mype = ", mype, &
      ", aco_ndata_lower(im) = ", aco_ndata_lower(im), ", disp = ", disp
    write (*, "(a)") trim(adjustl(msg))
  endif
  call MPI_File_set_view(fdat, disp, MPI_CHARACTER, MPI_CHARACTER, 'native', MPI_INFO_NULL, ierr)
  CHKERRQ(ierr)
  if (verbosity >= 2) then
    write (msg, "(3a)") "[alcon_output_data] Info: writing data for file ", &
      trim(alcon_output_filename(im, PETSC_TRUE, PETSC_TRUE)), " ...\n"
    call PetscPrintf(comm_radhead, trim(msg), ierr)
    CHKERRQ(ierr)
  endif
  do idata = 1, aco_ndata(im)
    write (charbuf, "(2e26.17e3, a)") aco_data(im, :, idata), char(10)
    call MPI_File_write(fdat, charbuf, nbyteperei, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    CHKERRQ(ierr)
  enddo
  if (verbosity >= 2) then
    write (msg, "(3a)") "[alcon_output_data] Info: closing file ", &
      trim(alcon_output_filename(im, PETSC_TRUE, PETSC_TRUE)), " ...\n"
    call PetscPrintf(comm_radhead, trim(msg), ierr)
    CHKERRQ(ierr)
  endif
  call MPI_File_close(fdat, ierr)
  CHKERRQ(ierr)
enddo

if (verbosity >= 2) then
  call PetscPrintf(comm_radhead, "[alcon_output_data] Info: has successfully written to files...\n", ierr)
  CHKERRQ(ierr)
endif

end subroutine alcon_output_data


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate extra output files, such as plotting script and Makefile !
! should only be called by root MPI process                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alcon_output_extra
implicit none
#include "finclude/petsc.h90"

PetscInt, parameter :: fmatlab = 203, fidl = 204, fgp = 205, fmakefile = 206
PetscInt, parameter :: ialcon_output_filenamestart = 8

PetscInt :: ncolor ! number of colors available in a specified terminal of gnuplot
PetscInt :: nlabeled ! number of m harmonics that are labeled in gnuplot
PetscInt :: imlast1, imlast2 ! index of the output file that has the largest index in acs_fdatout, imlast1 for first half of acs_fdatout, imlast2 for overall
character(2) :: linestyle ! linestyle number in gnuplot or color number in MATLAB
character(100) :: linetrailing ! trailing string of a line
character(3) :: string_m
character(3) :: string_mindex ! m indexing for generating MATLAB legend

PetscInt :: im, i

do im = m2 - m1, 0, -1
  if (aco_ndata(im) > 0) then
    imlast1 = im
    exit
  endif
enddo
if (finitebeta > 2) then
  do im = (m2 - m1 + 1) * 2 - 1, m2 - m1 + 1, -1
    if (aco_ndata(im) > 0) then
      imlast2 = im
      exit
    endif
  enddo
else
  imlast2 = imlast1
endif

select case (outformat)
  case (10) ! MATLAB
    ncolor = 13
    nlabeled = 13

    open (fmatlab, file = 'output/alcon.m', status = 'replace')
    do im = 0, imlast2
      if (aco_ndata(im) > 0) then
        write (fmatlab, "(4a)") "load '", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), "'"
      endif
    enddo
    write (fmatlab, "(a)") ""
    write (fmatlab, "(a)") "colortable=zeros(15,3);"
    write (fmatlab, "(a)") "colortable( 1,:)=[0.75 0.75 0.0 ];"
    write (fmatlab, "(a)") "colortable( 2,:)=[0.8  0.8  0.8 ];"
    write (fmatlab, "(a)") "colortable( 3,:)=[0.0  0.0  0.0 ];"
    write (fmatlab, "(a)") "colortable( 4,:)=[1.0  0.0  0.0 ];"
    write (fmatlab, "(a)") "colortable( 5,:)=[0.2  0.2  1.0 ];"
    write (fmatlab, "(a)") "colortable( 6,:)=[1.0  0.66 0.0 ];"
    write (fmatlab, "(a)") "colortable( 7,:)=[0.0  0.45 0.0 ];"
    write (fmatlab, "(a)") "colortable( 8,:)=[0.9  0.0  0.9 ];"
    write (fmatlab, "(a)") "colortable( 9,:)=[0.0  0.85 0.0 ];"
    write (fmatlab, "(a)") "colortable(10,:)=[0.0  0.85 0.85];"
    write (fmatlab, "(a)") "colortable(11,:)=[0.6  0.0  0.0 ];"
    write (fmatlab, "(a)") "colortable(12,:)=[0.4  0.7  0.4 ];"
    write (fmatlab, "(a)") "colortable(13,:)=[0.0  0.0  0.5 ];"
    write (fmatlab, "(a)") "colortable(14,:)=[0.6  0.0  0.6 ];"
    write (fmatlab, "(a)") "colortable(15,:)=[0.0  0.5  1.0 ];"
    write (fmatlab, "(a)") ""
    write (fmatlab, "(a)") "figure"
    write (fmatlab, "(a)") "hold on"
    if (finitebeta > 2) then
      do im = m2 - m1 + 1, min((m2 - m1 + 1) + (out_label_m1 - m1) - 1, imlast2)
        if (aco_ndata(im) > 0) then
          write (fmatlab, "(5a)") "scatter(", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
            "(:,1),", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
            "(:,2),1,colortable( 2,:),'.')"
        endif
      enddo

      do im = (m2 - m1 + 1) + out_label_m1 - m1, min((m2 - m1 + 1) + out_label_m1 - m1 + nlabeled - 1, imlast2)
        if (aco_ndata(im) > 0) then
          write (linestyle, '(i2)') min(3 + im - ((m2 - m1 + 1) + out_label_m1 - m1), 3 + ncolor - 1)
          write (fmatlab, "(7a)") "scatter(", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
            "(:,1),", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
            "(:,2),1,colortable(", linestyle, ",:),'.')"
        endif
      enddo

      do im = (m2 - m1 + 1) + out_label_m1 - m1 + nlabeled, imlast2
        if (aco_ndata(im) > 0) then
          write (fmatlab, "(5a)") "scatter(", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
            "(:,1),", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
            "(:,2),1,colortable( 2,:),'.')"
        endif
      enddo
    endif ! (finitebeta > 2)

    do im = 0, (out_label_m1 - m1) - 1
      if (aco_ndata(im) > 0) then
        write (fmatlab, "(5a)") "scatter(", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
          "(:,1),", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
          "(:,2),20,colortable( 1,:),'filled')"
      endif
    enddo

    i = 1
    do im = out_label_m1 - m1, min(out_label_m1 - m1 + nlabeled - 1, imlast1)
      if (aco_ndata(im) > 0) then
        write (linestyle, '(i2)') min(3 + im - (out_label_m1 - m1), 3 + ncolor - 1)
        write (string_mindex, "(i3)") i
        write (fmatlab, "(9a)") "hm(", string_mindex, ")=scatter(", &
          trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
          "(:,1),", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
          "(:,2),20,colortable(", linestyle, ",:),'filled');"
        i = i + 1
      endif
    enddo

    do im = out_label_m1 - m1 + nlabeled, imlast1
      if (aco_ndata(im) > 0) then
        write (fmatlab, "(5a)") "scatter(", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
          "(:,1),", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_FALSE)), &
          "(:,2),20,colortable( 1,:),'filled')"
      endif
    enddo
    write (fmatlab, "(a)") ""

    ! generate legend
    write (fmatlab, "(a)") "legend(hm, ..."
    do im = out_label_m1 - m1, min(out_label_m1 - m1 + nlabeled - 1, imlast1)
      if (aco_ndata(im) > 0) then
        write (string_m, "(i3)") im + m1
        write (fmatlab, "(3a)") "'m=", trim(adjustl(string_m)), "', ..."
      endif
    enddo
    write (fmatlab, "(a)") "'Location', 'NorthEastOutside')"
    ! set labels
    write (fmatlab, "(3a)") "xlabel('", out_xlabel, "')"
    write (fmatlab, "(3a)") "ylabel('", out_ylabel, "')"
    ! set plotting ranges
    if (omegacutoff > 0.0_kpr) then
      write (fmatlab, "(a, 4f20.16, a)") "axis([", out_rad1, out_rad2, 0.0_kpr, omegacutoff, "])"
    else
      write (fmatlab, "(a, 4f20.16, a)") "axis([", out_rad1, out_rad2, 0.0_kpr, 1.0_kpr, "])"
      write (fmatlab, "(a)") "axis 'auto y'"
    endif

    write (fmatlab, "(a)") ""
    close (fmatlab)


  case (20) ! IDL
    ncolor = 13
    nlabeled = 13

    open (fidl, file = 'output/alcon.pro', status = 'replace')
    write (fidl, "(a)") "pro alcon, eps = eps"
    write (fidl, "(a)") ""
    write (fidl, "(a)") "colortable = lonarr(16)"
    write (fidl, "(a)") "colortable( 0) = long(0.75 * 255l) + long(0.75 * 255 * 256l) + long(0.0  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable( 1) = long(0.8  * 255l) + long(0.8  * 255 * 256l) + long(0.8  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable( 2) = long(1.0  * 255l) + long(1.0  * 255 * 256l) + long(1.0  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable( 3) = long(0.0  * 255l) + long(0.0  * 255 * 256l) + long(0.0  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable( 4) = long(1.0  * 255l) + long(0.0  * 255 * 256l) + long(0.0  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable( 5) = long(0.2  * 255l) + long(0.2  * 255 * 256l) + long(1.0  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable( 6) = long(1.0  * 255l) + long(0.66 * 255 * 256l) + long(0.0  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable( 7) = long(0.0  * 255l) + long(0.45 * 255 * 256l) + long(0.0  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable( 8) = long(0.9  * 255l) + long(0.0  * 255 * 256l) + long(0.9  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable( 9) = long(0.0  * 255l) + long(0.85 * 255 * 256l) + long(0.0  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable(10) = long(0.0  * 255l) + long(0.85 * 255 * 256l) + long(0.85 * 255 * 256l^2)"
    write (fidl, "(a)") "colortable(11) = long(0.6  * 255l) + long(0.0  * 255 * 256l) + long(0.0  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable(12) = long(0.4  * 255l) + long(0.7  * 255 * 256l) + long(0.4  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable(13) = long(0.0  * 255l) + long(0.0  * 255 * 256l) + long(0.5  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable(14) = long(0.6  * 255l) + long(0.0  * 255 * 256l) + long(0.6  * 255 * 256l^2)"
    write (fidl, "(a)") "colortable(15) = long(0.0  * 255l) + long(0.5  * 255 * 256l) + long(1.0  * 255 * 256l^2)"
    write (fidl, "(a)") ""
    write (fidl, "(a)") "if keyword_set(eps) then begin"
    write (fidl, "(2a)") char(9), "old_device = !d.name"
    write (fidl, "(2a)") char(9), "set_plot, 'PS'"
    write (fidl, "(2a)") char(9), "device, file = 'alcon.eps', /encapsulated, /color, bits=8"
    write (fidl, "(2a)") char(9), "tvlct, 0.75 * 255, 0.75 * 255, 0.0  * 255,  0"
    write (fidl, "(2a)") char(9), "tvlct, 0.8  * 255, 0.8  * 255, 0.8  * 255,  1"
    write (fidl, "(2a)") char(9), "tvlct, 1.0  * 255, 1.0  * 255, 1.0  * 255,  2"
    write (fidl, "(2a)") char(9), "tvlct, 0.0  * 255, 0.0  * 255, 0.0  * 255,  3"
    write (fidl, "(2a)") char(9), "tvlct, 1.0  * 255, 0.0  * 255, 0.0  * 255,  4"
    write (fidl, "(2a)") char(9), "tvlct, 0.2  * 255, 0.2  * 255, 1.0  * 255,  5"
    write (fidl, "(2a)") char(9), "tvlct, 1.0  * 255, 0.66 * 255, 0.0  * 255,  6"
    write (fidl, "(2a)") char(9), "tvlct, 0.0  * 255, 0.45 * 255, 0.0  * 255,  7"
    write (fidl, "(2a)") char(9), "tvlct, 0.9  * 255, 0.0  * 255, 0.9  * 255,  8"
    write (fidl, "(2a)") char(9), "tvlct, 0.0  * 255, 0.85 * 255, 0.0  * 255,  9"
    write (fidl, "(2a)") char(9), "tvlct, 0.0  * 255, 0.85 * 255, 0.85 * 255, 10"
    write (fidl, "(2a)") char(9), "tvlct, 0.6  * 255, 0.0  * 255, 0.0  * 255, 11"
    write (fidl, "(2a)") char(9), "tvlct, 0.4  * 255, 0.7  * 255, 0.4  * 255, 12"
    write (fidl, "(2a)") char(9), "tvlct, 0.0  * 255, 0.0  * 255, 0.5  * 255, 13"
    write (fidl, "(2a)") char(9), "tvlct, 0.6  * 255, 0.0  * 255, 0.6  * 255, 14"
    write (fidl, "(2a)") char(9), "tvlct, 0.0  * 255, 0.5  * 255, 1.0  * 255, 15"
    write (fidl, "(2a)") char(9), "colortable = lindgen(16)"
    write (fidl, "(a)") "endif"
    write (fidl, "(a)") ""
    write (fidl, "(a)") "lun = 124"
    if (omegacutoff <= 0.0_kpr) then
      write (fidl, "(a, i3, a)") "arrayomegamax = dblarr(", imlast2 + 1, ")"
    endif
    do im = 0, imlast2
      if (aco_ndata(im) > 0) then
        write (fidl, "(3a)") "nlines = file_lines('", &
          trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), "', /noexpand_path)"
        write (string_mindex, "(i3)") im
        write (fidl, "(3a)") "data", trim(adjustl(string_mindex)), " = dblarr(2, nlines)"
        write (fidl, "(3a)") "openr, lun, '", &
          trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), "'"
        write (fidl, "(2a)") "readf, lun, data", trim(adjustl(string_mindex))
        write (fidl, "(a)") "close, lun"
        if (omegacutoff <= 0.0_kpr) then
          write (fidl, "(5a)") "arrayomegamax(", trim(adjustl(string_mindex)),") = max(data", &
            trim(adjustl(string_mindex)), "(1, *))"
        endif
      endif
    enddo
    write (fidl, "(a)") ""
    write (fidl, "(a)") "t = findgen(17) * (2.0 * !pi / 16.0)"
    write (fidl, "(a)") "usersym, 0.3 * cos(t), 0.3 * sin(t), /fill"
    write (fidl, "(a)") ""
    write (fidl, "(a)") "!p.psym = 8"
    write (fidl, "(a)") "!p.symsize = 1.0"
    write (fidl, "(a)") "!p.background = colortable(2)"
    write (fidl, "(a)") ""
    if (omegacutoff > 0.0_kpr) then
      write (linetrailing, "(a, f20.16, a)") "yrange = [0.0, ", omegacutoff, "], "
    else
      write (fidl, "(a)") "omegamax = max(arrayomegamax)"
      write (linetrailing, "(a)") "yrange = [0.0, omegamax], "
    endif
    write (fidl, "(a, f20.16, a, f20.16, 7a)") "plot, [-0.1], [-0.1], xrange = [", &
      out_rad1, ", ", out_rad2, &
      "], ", trim(linetrailing), "xtitle = '", out_xlabel, "', ytitle = '", &
      out_ylabel, "', color = colortable(3)"
    if (finitebeta > 2) then
      do im = m2 - m1 + 1, min((m2 - m1 + 1) + (out_label_m1 - m1) - 1, imlast2)
        if (aco_ndata(im) > 0) then
          write (string_mindex, "(i3)") im
          write (fidl, "(5a)") "oplot, data", trim(adjustl(string_mindex)), "(0, *), data", &
            trim(adjustl(string_mindex)), "(1, *), color = colortable( 1)"
        endif
      enddo

      do im = (m2 - m1 + 1) + out_label_m1 - m1, min((m2 - m1 + 1) + out_label_m1 - m1 + nlabeled - 1, imlast2)
        if (aco_ndata(im) > 0) then
          write (string_mindex, "(i3)") im
          write (linestyle, "(i2)") 3 + im - ((m2 - m1 + 1) + out_label_m1 - m1)
          write (fidl, "(7a)") "oplot, data", trim(adjustl(string_mindex)), "(0, *), data", &
            trim(adjustl(string_mindex)), "(1, *), color = colortable(", trim(linestyle), ")"
        endif
      enddo

      do im = (m2 - m1 + 1) + out_label_m1 - m1 + nlabeled, imlast2
        if (aco_ndata(im) > 0) then
          write (string_mindex, "(i3)") im
          write (fidl, "(5a)") "oplot, data", trim(adjustl(string_mindex)), "(0, *), data", &
            trim(adjustl(string_mindex)), "(1, *), color = colortable( 1)"
        endif
      enddo
    endif

    do im = 0, min((out_label_m1 - m1) - 1, imlast1)
      if (aco_ndata(im) > 0) then
        write (string_mindex, "(i3)") im
        write (fidl, "(5a)") "oplot, data", trim(adjustl(string_mindex)), "(0, *), data", &
          trim(adjustl(string_mindex)), "(1, *), color = colortable( 0), symsize = 3.0"
      endif
    enddo

    do im = out_label_m1 - m1, min(out_label_m1 - m1 + nlabeled - 1, imlast1)
      if (aco_ndata(im) > 0) then
        write (string_mindex, "(i3)") im
        write (linestyle, "(i2)") 3 + im - (out_label_m1 - m1)
        write (fidl, "(7a)") "oplot, data", trim(adjustl(string_mindex)), "(0, *), data", &
          trim(adjustl(string_mindex)), "(1, *), color = colortable(", trim(linestyle), "), symsize = 3.0"
      endif
    enddo

    do im = out_label_m1 - m1 + nlabeled, imlast1
      if (aco_ndata(im) > 0) then
        write (string_mindex, "(i3)") im
        write (fidl, "(5a)") "oplot, data", trim(adjustl(string_mindex)), "(0, *), data", &
          trim(adjustl(string_mindex)), "(1, *), color = colortable( 0), symsize = 3.0"
      endif
    enddo
    write (fidl, "(a)") ""

    write (fidl, "(a)") "; put label at the data average position"
    do im = out_label_m1 - m1, min(out_label_m1 - m1 + nlabeled - 1, imlast1)
      if (aco_ndata(im) > 0) then
        write (string_mindex, "(i3)") im
        write (string_m, "(i3)") im + m1
        write (linestyle, "(i2)") 3 + im - (out_label_m1 - m1)
        write (fidl, "(9a)") "xyouts, mean(data", trim(adjustl(string_mindex)), &
          "(0, *)), mean(data", trim(adjustl(string_mindex)), "(1, *)), 'm=", &
          trim(adjustl(string_m)), "', color = colortable(", linestyle, ")"
      endif
    enddo
    write (fidl, "(a)") ""
    write (fidl, "(a)") "if keyword_set(eps) then begin"
    write (fidl, "(2a)") char(9), "device, /close"
    write (fidl, "(2a)") char(9), "set_plot, old_device"
    write (fidl, "(a)") "endif"
    write (fidl, "(a)") ""

    write (fidl, "(a)") "end"
    write (fidl, "(a)") ""
    close (fidl)


  case (30 : 33) ! gnuplot
    ncolor = 13
    nlabeled = 13

    ! generate gp script
    open (fgp, file = 'output/alcon.gp', status = 'replace')
    if (outformat == 30 .or. outformat == 31) then
      write (fgp, "(a)") "set terminal postscript enhanced color"
      write (fgp, "(a)") "set output 'alcon.eps'"
      write (fgp, "(a)") "set style line 11 lc rgbcolor '#bfbf00' pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 12 lc 1 pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 13 lc rgbcolor '#007300' pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 14 lc 3 pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 15 lc 4 pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 16 lc 7 pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 17 lc 8 pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 18 lc rgbcolor '#66b266' pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 19 lc rgbcolor '#000080' pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 20 lc 6 pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 21 lc rgbcolor '#990099' pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 22 lc rgbcolor '#0080ff' pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 23 lc 2 pt 7 ps 0.6"
      write (fgp, "(a)") "set style line 24 lc 5 pt 7 ps 0.6"
      if (finitebeta > 2) then
        write (fgp, "(a)") "set style line 41 lc rgbcolor '#cccccc' pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 42 lc 1 pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 43 lc rgbcolor '#007300' pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 44 lc 3 pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 45 lc 4 pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 46 lc 7 pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 47 lc 8 pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 48 lc rgbcolor '#66b266' pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 49 lc rgbcolor '#000080' pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 50 lc 6 pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 51 lc rgbcolor '#990099' pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 52 lc rgbcolor '#0080ff' pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 53 lc 2 pt 7 ps 0.2"
        write (fgp, "(a)") "set style line 54 lc 5 pt 7 ps 0.2"
      endif
      write (fgp, "(3a)") "set xlabel '", out_xlabel, "' offset 0, 1"
      write (fgp, "(3a)") "set ylabel '", out_ylabel, "' offset 2, 0"
      write (fgp, "(a)") "set xtics offset 0, 0.4"
      write (fgp, "(a)") "set ytics offset 0.5, 0"
    else ! if (outformat == 30 .or. outformat == 31); implying (outformat == 32 or 33)
      write (fgp, "(a)") "set terminal mp color"
      write (fgp, "(a)") "set output 'alcon.mp'"
      write (fgp, "(a)") "set style line 11 lc rgbcolor '#bfbf00' pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 12 lc 1 pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 13 lc 2 pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 14 lc 3 pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 15 lc 4 pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 16 lc rgbcolor '#007300' pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 17 lc 6 pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 18 lc 5 pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 19 lc 7 pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 20 lc rgbcolor '#990000' pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 21 lc rgbcolor '#66b266' pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 22 lc rgbcolor '#000080' pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 23 lc rgbcolor '#990099' pt 9 ps 0.6"
      write (fgp, "(a)") "set style line 24 lc rgbcolor '#0080ff' pt 9 ps 0.6"
      if (finitebeta > 2) then
        write (fgp, "(a)") "set style line 41 lc rgbcolor '#cccccc' pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 42 lc 1 pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 43 lc 2 pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 44 lc 3 pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 45 lc 4 pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 46 lc rgbcolor '#007300' pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 47 lc 6 pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 48 lc 5 pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 49 lc 7 pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 50 lc rgbcolor '#990000' pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 51 lc rgbcolor '#66b266' pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 52 lc rgbcolor '#000080' pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 53 lc rgbcolor '#990099' pt 9 ps 0.001"
        write (fgp, "(a)") "set style line 54 lc rgbcolor '#0080ff' pt 9 ps 0.001"
      endif
      write (fgp, "(3a)") "set xlabel '", out_xlabel, "' offset 0, 1.1"
      write (fgp, "(3a)") "set ylabel '", out_ylabel, "' offset 2.5, 0"
      write (fgp, "(a)") "set xtics offset 0, 0.5"
      write (fgp, "(a)") "set ytics offset 0.5, 0"
    endif ! else of if (outformat == 30 .or. outformat == 31)
    write (fgp, "(a, f20.16, a, f20.16, a)") "set xrange [", out_rad1, ":", out_rad2, "]"
    if (omegacutoff > 0.0_kpr) then
      write (fgp, "(a, f20.16, a)") "set yrange [0:", omegacutoff, "]"
    endif
    write (fgp, "(a)") "set key outside"
    write (fgp, "(a)") "set key textcolor rgb variable"
    write (fgp, "(a)") "plot \"
    if (finitebeta > 2) then
      do im = m2 - m1 + 1, min((m2 - m1 + 1) + (out_label_m1 - m1) - 1, imlast2)
        if (aco_ndata(im) > 0) then
          write (fgp, "(3a)") "'", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), &
            "' notitle with points ls 41, \"
        endif
      enddo

      do im = (m2 - m1 + 1) + out_label_m1 - m1, min((m2 - m1 + 1) + out_label_m1 - m1 + nlabeled - 1, imlast2)
        if (aco_ndata(im) > 0) then
          write (linestyle, '(i2)') min(42 + im - ((m2 - m1 + 1) + out_label_m1 - m1), 42 + ncolor - 1)
          write (fgp, "(5a)") "'", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), &
            "' notitle with points ls ", linestyle, ", \"
        endif
      enddo

      do im = (m2 - m1 + 1) + out_label_m1 - m1 + nlabeled, imlast2
        if (aco_ndata(im) > 0) then
          write (fgp, "(3a)") "'", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), &
            "' notitle with points ls 41, \"
        endif
      enddo
    endif ! (finitebeta > 2)

    do im = 0, min((out_label_m1 - m1) - 1, imlast1)
      if (aco_ndata(im) > 0) then
        if (im < imlast1) then
          linetrailing = ", \"
        else
          linetrailing = ""
        endif
        write (fgp, "(4a)") "'", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), &
          "' notitle with points ls 11", trim(linetrailing)
      endif
    enddo

    do im = out_label_m1 - m1, min(out_label_m1 - m1 + nlabeled - 1, imlast1)
      if (aco_ndata(im) > 0) then
        if (im < imlast1) then
          linetrailing = ", \"
        else
          linetrailing = ""
        endif
        write (linestyle, '(i2)') min(12 + im - (out_label_m1 - m1), 12 + ncolor - 1)
        write (string_m, "(i3)") im + m1
        write (fgp, "(7a)") "'", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), &
          "' title 'm = ", trim(adjustl(string_m)), &
          "' with points ls ", linestyle, trim(linetrailing)
      endif
    enddo

    do im = out_label_m1 - m1 + nlabeled, imlast1
      if (aco_ndata(im) > 0) then
        if (im < imlast1) then
          linetrailing = ", \"
        else
          linetrailing = ""
        endif
        write (fgp, "(4a)") "'", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), &
          "' notitle with points ls 11", trim(linetrailing)
      endif
    enddo

    write (fgp, "(a)") ""

    close (fgp)

    ! generate Makefile
    open (fmakefile, file = 'output/Makefile', status = 'replace')

    if (outformat == 30 .or. outformat == 31) then
      write (fmakefile, "(a)") ".PHONY : clean"
    else
      write (fmakefile, "(a)") ".PHONY : clean cleanall"
      write (fmakefile, "(a)") ""
      write (fmakefile, "(a)") "alcon.mps : alcon.mp"
      write (fmakefile, "(2a)") char(9), "mpost $<"
      write (fmakefile, "(2a)") char(9), "mv alcon.0 $@"
      write (fmakefile, "(a)") ""
      write (fmakefile, "(a)") "alcon.eps : alcon.mps"
      write (fmakefile, "(2a)") char(9), "mps2eps $<"
      write (fmakefile, "(a)") ""
      write (fmakefile, "(a)") "alcon.pdf : alcon.mps"
      write (fmakefile, "(2a)") char(9), "mptopdf $<"
      write (fmakefile, "(2a)") char(9), "mv alcon-mps.pdf $@"
      write (fmakefile, "(a)") ""
      write (fmakefile, "(a)") "alcon.png : alcon.eps"
      write (fmakefile, "(2a)") char(9), "convert -density 400 -flatten $< $@"
    endif
    write (fmakefile, "(a)") ""
    if (outformat == 30 .or. outformat == 32) then
      if (outformat == 30) then
        write (fmakefile, "(a)") "alcon.eps : alcon.gp \"
      else
        write (fmakefile, "(a)") "alcon.mp : alcon.gp \"
      endif
      do im = 0, imlast2
        if (aco_ndata(im) > 0) then
          if (im < imlast2) then
            linetrailing = " \"
          else
            linetrailing = ""
          endif
          write (fmakefile, "(3a)") " ", trim(alcon_output_filename(im, PETSC_FALSE, PETSC_TRUE)), trim(linetrailing)
        endif
      enddo
    else ! implying outformat == 31 or 33
      if (outformat == 31) then
        write (fmakefile, "(a)") "alcon.eps : alcon.gp alcon.d"
      else ! implying outformat == 33
        write (fmakefile, "(a)") "alcon.mp : alcon.gp alcon.d"
      endif
    endif
    write (fmakefile, "(2a)") char(9), "gnuplot $<"
    if (outformat == 33) then
      write (fmakefile, "(2a)") char(9), "gpmp2latexmp $@"
    endif
    write (fmakefile, "(a)") ""
    if (outformat == 31 .or. outformat == 33) then
      write (fmakefile, "(a)") "alcon.d : alcon.gp"
      write (fmakefile, "(2a)") char(9), "makedependgp $<"
      write (fmakefile, "(a)") ""
      write (fmakefile, "(a)") "-include alcon.d"
      write (fmakefile, "(a)") ""
    endif
    write (fmakefile, "(a)") "clean :"
    if (outformat == 31 .or. outformat == 33) then
      linetrailing = " alcon.d"
    else
      linetrailing = ""
    endif
    if (outformat == 30 .or. outformat == 31) then
      write (fmakefile, "(3a)") char(9), "rm -f alcon.eps", trim(linetrailing)
    else
      write (fmakefile, "(3a)") char(9), "rm -f mpxerr.tex mpxerr.log alcon.0 alcon.mp~ alcon.mpx alcon.log", trim(linetrailing)
      write (fmakefile, "(2a)") char(9), &
        "rm -f alcon.mp"
      write (fmakefile, "(a)") ""
      write (fmakefile, "(a)") "cleanall : clean"
      write (fmakefile, "(2a)") char(9), "rm -f alcon.mps alcon.eps alcon.png"
    endif
    write (fmakefile, "(a)") ""

    close (fmakefile)
  case default
    return
endselect
end subroutine alcon_output_extra

end module alcon_output
