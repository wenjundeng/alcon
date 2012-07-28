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


! this is a sample subroutine to show how to generate an alcon.dat file.
subroutine acdgen
implicit none

! kind # to get double precision
integer, parameter :: kint = selected_int_kind(10), kreal = selected_real_kind(10)

real(kind = kreal), parameter :: pi = 3.14159265358979323846264d0

! gamma for ions and electrons
real(kind = kreal), parameter :: gammai = 7d0 / 4d0, gammae = 1d0

! mass ratio of ion to proton
real(kind = kreal), parameter :: mass_ratio_of_ion_to_proton = 1d0

! on-axis magnetic field (Gauss)
real(kind = kreal), parameter :: b0 = 1.91d4

! on-axis major radius (cm)
real(kind = kreal), parameter :: r0 = 1.4194d2

! poloidal flux at wall
real(kind = kreal), parameter :: psiw = 2.59808d-2

! number of radial grid points
! since ALCON will use cubic spline interpolation,
! the radial resolution in alcon.dat does not need to be very high.
! usually 200 is good enough
integer(kind = kint), parameter :: nrad = 200

! number of poloidal grid points in the sampling for FFT
integer(kind = kint), parameter :: npol = 2560

! for the Fourier components in Eqs. (A.23)--(A.26) and (A.43) in
!   [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117), up to nfft
!   will be written to alcon.dat
integer(kind = kint), parameter :: nfft = 100

! number of profiles and number of FFT data sets
integer(kind = kint), parameter :: nprofiledata = 5, nfftdata = 5

! file number
integer(kind = kint), parameter :: falcondat = 240

! two major arrays to write to alcon.dat
real(kind = kreal), dimension(nprofiledata, nrad) :: acdprofile ! profile data
complex(kind = kreal), dimension(0 : nfft, nfftdata, nrad) :: acdfft ! Fourier component data

! sampling over poloidal direction for matrixes H, J, K, L, and N.
! these matrixes are defined in Eqs. (A.23)--(A.26) and (A.42) in
!   [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117).
real(kind = kreal), dimension(0 : npol - 1) :: hsamp, jsamp, ksamp, lsamp, nsamp

! for temporary storage of the Fourier components obtained by performing FFT on
!   the above arrays
! since FFT is performed on real arrays, only the first half of the Fourier
!   components are meaningful.
complex(kind = kreal), dimension(0 : npol / 2) :: sampfft

integer(kind = kint) :: irad, ipol
real(kind = kreal) :: psi, theta

! external functions containing equilibrium information are needed:
! (here it is assumed that the poloidal flux psi is the basic radial coordinate
!   in your code)
! rho(psi): rho is the square root of the normalized (to wall) toroidal flux
! q(psi): safety factor q
! g(psi), ci(psi): current functions g and I defined in Eq. (A.8) in
!   [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117)
! ni(psi), ti(psi), ne(psi), te(psi): ion and electron densities and temperatures
! b(psi, theta): magnetic field amplitude, where theta is the magnetic poloidal
!   coordinate
! dbdt(psi, theta): partial b / partial theta
! gradpsi(psi, theta): magnitude of gradient of poloidal flux
! fftr1d(real_array_1d): performs FFT on a real 1D array
! these function should return values in CGS units
real(kind = kint), external :: rho, q, g, ci, ni, ti, ne, te, b, dbdt, gradpsi, fftr1d

! initialization
acdprofile = 0d0
acdfft = 0d0

do irad = 1, nrad
  ! this way the radial grids are uniform in psi, psi ranging in
  !   [1d-4 * psiw, psiw], which may have too sparse grids in real space near
  !   magnetic axis.
  psi = psiw * (1d-4 + (1d0 - 1d-4) * (irad - 1) / (nrad - 1))
  ! this alternative way generates radial grids more uniform in real space.
  psi = psiw * (1d-2 + (1d0 - 1d-2) * (irad - 1) / (nrad - 1))**2

  acdprofile(1, irad) = rho(psi)
  acdprofile(2, irad) = q(psi)
  ! for Jacobian calculation
  acdprofile(3, irad) = g(psi) * q(psi) + ci(psi)
  ! pressure profile
  acdprofile(4, irad) = 4 * pi * &
    (gammai * ni(psi) * ti(psi) + gammae * ne(psi) * te(psi)) / b0**2
  ! mass density
  acdprofile(5, irad) = (ni(psi) / ne(0d0)) * mass_ratio_of_ion_to_proton

  do ipol = 0, npol - 1
    theta = 2d0 * pi * ipol / npol
    ! operators H, J, K, L, and N defined in Eqs. (A.23)--(A.26) and (A.42) in
    !   [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117).
    ! note that the last factors in the following lines are the normalization
    !   factors in Eqs. (A.28)--(A.31) and (A.44)
    hsamp(ipol) = (gradpsi(psi, theta)**2 / acdprofile(3, irad)) / (b0 * r0)
    jsamp(ipol) = (gradpsi(psi, theta)**2 * acdprofile(3, irad) / (b(psi, theta)**4)) * (b0 / r0**3)
    ksamp(ipol) = (2d0 * g(psi) * dbdt(psi, theta) / (b(psi, theta)**3)) * (b0 / r0)
    lsamp(ipol) = ((acdprofile(4, irad) + b(psi, theta)**2) * acdprofile(3, irad) &
      / b(psi, theta)**4) * (b0 / r0)
    nsamp(ipol) = (acdprofile(4, irad) * ksamp(ipol)**2 / lsamp(ipol)) / (b0 * r0)
  end do
  ! the "/ npol" factors in the following lines are for normalization.
  ! there may be an additional factor needed to be divided by, such as 2 * pi.
  ! check your fftr1d function, make sure after normalization, performing it
  !   on cos(theta) would give 1 on the first Fourier component.
  sampfft = fftr1d(hsamp)
  acdfft(0 : nfft, 1, irad) = sampfft(0 : nfft) / npol
  sampfft = fftr1d(jsamp)
  acdfft(0 : nfft, 2, irad) = sampfft(0 : nfft) / npol
  sampfft = fftr1d(ksamp)
  acdfft(0 : nfft, 3, irad) = sampfft(0 : nfft) / npol
  sampfft = fftr1d(lsamp)
  acdfft(0 : nfft, 4, irad) = sampfft(0 : nfft) / npol
  sampfft = fftr1d(nsamp)
  acdfft(0 : nfft, 5, irad) = sampfft(0 : nfft) / npol
end do

! note that if in a parallel code, the writing part should be executed by only
!   one process.
! also note that here acdprofile and acdfft are written in column-major order
!   as Fortran deals with multi-dimension arrays in this manner.  if you use
!   another language to generate alcon.dat, make sure you still write the data
!   in column-major order.
open (falcondat, file = 'alcon.dat', status = 'replace')
write (falcondat, *) nrad, nfft, nprofiledata, nfftdata
write (falcondat, *) acdprofile
write (falcondat, *) acdfft
close (falcondat)

end subroutine acdgen

