! Copyright 2011-2014 Wenjun Deng <wdeng@wdeng.info>
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


! This is a sample subroutine to show how to generate an alcon.dat file.
subroutine acdgen
implicit none

! kind # to get double precision
integer, parameter :: kint = selected_int_kind(10), kreal = selected_real_kind(10)

real(kind = kreal), parameter :: pi = 3.14159265358979323846264d0

! gamma for ions and electrons
real(kind = kreal), parameter :: gammai = 7d0 / 4d0, gammae = 1d0

! Ion mass (normalized to proton mass)
real(kind = kreal), parameter :: massi = 1d0

! On-axis magnetic field (Gauss)
real(kind = kreal), parameter :: bb0 = 1.91d4

! On-axis major radius (cm)
real(kind = kreal), parameter :: rr0 = 1.4194d2

! Poloidal flux at wall
real(kind = kreal), parameter :: psiw = 2.59808d-2

! Number of radial (psi) grid points
! Since ALCON will use cubic spline interpolation, the radial resolution in
! alcon.dat does not need to be very high.  Usually 200 is good enough.
integer(kind = kint), parameter :: nrad = 200

! Number of poloidal (theta) grid points in the sampling for FFT
integer(kind = kint), parameter :: ntheta = 2048

! Number of FFT coefficents to be saved for the Fourier components in
! Eqs. (A.23)--(A.26) and (A.43) in
! [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117).
integer(kind = kint), parameter :: nfftcoef = 100

! Number of profiles and number of FFT data sets
integer(kind = kint), parameter :: nprofile = 5, nfft = 5

! File number
integer(kind = kint), parameter :: falcondat = 240

! Two major arrays to write to alcon.dat
real(kind = kreal), dimension(nprofile, nrad) :: profiles ! profile data
complex(kind = kreal), dimension(0:nfftcoef, nfft, nrad) :: fftcoefs ! FFT coefficent data

! Samplings over poloidal (theta) direction for matrixes H, J, K, L, and N.
! These matrixes are defined in Eqs. (A.23)--(A.26) and (A.42) in
! [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117).
real(kind = kreal), dimension(0:(ntheta - 1)) :: hhsamp, jjsamp, kksamp, llsamp, nnsamp

! For temporary storage of the Fourier components obtained by performing FFT on
! the above arrays.  Since FFT is performed on real arrays, only the first half
! of the coefficients are meaningful.
complex(kind = kreal), dimension(0:(ntheta / 2)) :: sampfft

integer(kind = kint) :: irad, itheta
real(kind = kreal) :: psi, theta

! External functions containing equilibrium information are needed:
! (here it is assumed that the poloidal flux psi is the basic radial coordinate
!   in your code)
! rho(psi): rho is the square root of the normalized (to wall) toroidal flux
! q(psi): safety factor q
! g(psi), ii(psi): current functions g and I defined in Eq. (A.8) in
!   [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117)
! deni(psi), tti(psi), dene(psi), tte(psi): ion and electron densities and
!   temperatures
! bb(psi, theta): magnetic field amplitude, where theta is the magnetic
!   poloidal coordinate
! dbbdtheta(psi, theta): partial b / partial theta
! gradpsi(psi, theta): magnitude of gradient of poloidal flux
! fftr1d(real_array_1d): performs FFT on a real 1D array
! In this sample subroutine, these functions should return values in CGS units
real(kind = kreal), external :: rho, q, g, ii, deni, tti, dene, tte, bb, &
  dbbdtheta, gradpsi, fftr1d

! Initialization
profiles = 0d0
fftcoefs = 0d0

do irad = 1, nrad
  ! This way the radial grids are uniform in psi, psi ranging in
  ! [1d-4 * psiw, psiw], which may have too sparse grids in real space near
  ! magnetic axis.
  psi = psiw * (1d-4 + (1d0 - 1d-4) * (irad - 1) / (nrad - 1))

  ! This alternative way generates radial grids more uniform in real space.
  psi = psiw * (1d-2 + (1d0 - 1d-2) * (irad - 1) / (nrad - 1))**2

  ! 1st profile is radial coordinate
  profiles(1, irad) = rho(psi)

  ! 2nd profile is q-profile
  profiles(2, irad) = q(psi)

  ! 3rd profile is "g q + I" for Jacobian calculation
  profiles(3, irad) = g(psi) * q(psi) + ii(psi)

  ! 4th profile is pressure
  profiles(4, irad) = 4 * pi * &
    (gammai * deni(psi) * tti(psi) + gammae * dene(psi) * tte(psi)) / bb0**2

  ! 5th profile is mass density
  profiles(5, irad) = (deni(psi) / dene(0d0)) * massi

  ! Sampling operators H, J, K, L, and N along poloidal (theta) direction
  ! defined in Eqs. (A.23)--(A.26) and (A.42) in
  ! [Nuclear Fusion 52, 043006 (2012)](http://wdeng.info/?p=117).
  do itheta = 0, ntheta - 1
    theta = 2d0 * pi * itheta / ntheta
    ! Note that the last factors in the following lines are for normalization
    ! factors in Eqs. (A.28)--(A.31) and (A.44)

    ! Sampling H
    hhsamp(itheta) = (gradpsi(psi, theta)**2 / profiles(3, irad)) / (bb0 * rr0)

    ! Sampling J
    jjsamp(itheta) = ( &
      gradpsi(psi, theta)**2 * profiles(3, irad) / (bb(psi, theta)**4) &
    ) * (bb0 / rr0**3)

    ! Sampling K
    kksamp(itheta) = ( &
      2d0 * g(psi) * dbbdtheta(psi, theta) / (bb(psi, theta)**3) &
    ) * (bb0 / rr0)

    ! Sampling L
    llsamp(itheta) = ( &
      (profiles(4, irad) + bb(psi, theta)**2) * profiles(3, irad) &
      / bb(psi, theta)**4 &
    ) * (bb0 / rr0)

    ! Sampling N
    nnsamp(itheta) = (profiles(4, irad) * kksamp(itheta)**2 / llsamp(itheta)) &
      / (bb0 * rr0)
  end do

  ! Perform FFT on the operator samples to get the matrix elements.
  ! The "/ ntheta" factors in the following lines are for normalization.
  ! There may be an additional factor needed to be divided by, such as 2 * pi.
  ! Check your fftr1d function, make sure after normalization, performing it
  ! on cos(theta) would give 1/2 on the first (right after the zeroth)
  ! Fourier coefficient.

  ! Matrix elements for H
  sampfft = fftr1d(hhsamp)
  fftcoefs(0:nfftcoef, 1, irad) = sampfft(0:nfftcoef) / ntheta

  ! Matrix elements for J
  sampfft = fftr1d(jjsamp)
  fftcoefs(0:nfftcoef, 2, irad) = sampfft(0:nfftcoef) / ntheta

  ! Matrix elements for K
  sampfft = fftr1d(kksamp)
  fftcoefs(0:nfftcoef, 3, irad) = sampfft(0:nfftcoef) / ntheta

  ! Matrix elements for L
  sampfft = fftr1d(llsamp)
  fftcoefs(0:nfftcoef, 4, irad) = sampfft(0:nfftcoef) / ntheta

  ! Matrix elements for N
  sampfft = fftr1d(nnsamp)
  fftcoefs(0:nfftcoef, 5, irad) = sampfft(0:nfftcoef) / ntheta
end do

! Note that if in a parallel code, the writing part should be executed by only
! one process, unless you parallelize it properly.
! Also note that here profiles and fftcoefs are written in column-major order
! as Fortran deals with multi-dimension arrays in this manner.  If you use
! another language to generate alcon.dat, make sure you still write the data
! in column-major order.
open (falcondat, file = 'alcon.dat', status = 'replace')
write (falcondat, *) nrad, nfftcoef, nprofile, nfft
write (falcondat, *) profiles
write (falcondat, *) fftcoefs
close (falcondat)

end subroutine acdgen
