'''This file deals with equilibrium data.'''

# Copyright 2014 Wenjun Deng <wdeng@wdeng.info>
#
# This file is part of ALCON.
#
# ALCON is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ALCON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ALCON.  If not, see <http://www.gnu.org/licenses/>.

import sys
write = sys.stdout.write
import re
import numpy as np
from scipy.interpolate import interp1d

class EqData:
    '''Equilibrium data class'''

    # Number of profiles and FFT data sets
    nprofile = 5
    nfft = 5

    def __init__(self, input_):
        self.input_ = input_

        # Index for nprofile:
        # 0: radial coordinate (usually rho); 1: q; 2: g q + I; 3: beta;
        # 4: rho_M
        self.profiles = np.zeros((self.nprofile, input_.nrad))

        # Generate uniform radial grids
        self.profiles[0, :] = np.linspace(
            input_.radrange[0], input_.radrange[1], input_.nrad)

        # Correspondence between index for nfft and matrix:
        # 0: H; 1: J; 2: K; 3: L; 4: N
        # For details about these matrixes, see Eqs. (A.23)--(A.26) and (A.42)
        # in [Nulcear Fusion 52, 043006 (2012)]
        self.fftcoefs = np.zeros(
            (input_.noffdiag + 1, self.nfft, input_.nrad), dtype = complex)

        if input_.eqtype == 'sa':
            self.load_sa()
        #elif input_.eqtype == 'sppro':
            #self.load_sppro()
        else:
            self.load_acd()
    # End of def __init__(self, input_):

    def load_sa(self):
        '''Load a simple analytic equilibrium model described in Appendix E in
        [Nuclear Fusion 52, 023005 (2012)] (extended to second order), load
        profiles from profile_sa.dat.
        This second order equilibrium is:
        B = 1 - epsilon * cos(theta0) + epsilon**2 / 2 * (1/q**2 + cos(2 theta0)) + O(epsilon**3)
          = 1 - epsilon * cos(theta) + epsilon**2 / 2 * (1/q**2 + 1) + O(epsilon**3)
        g = 1 - epsilon**2 / 2 + O(epsilon**3)
        I = epsilon**2 / q + epsilon**4 / (2 q) + O(epsilon**5)
        delta = epsilon * sin(theta0) - epsilon**2 / 2 * sin(2 theta0) + O(epsilon**3)
        theta = theta0 - epsilon * sin(theta0) + epsilon**2 / 4 * sin(2 theta0) + O(epsilon**3)
        theta0 = theta + epsilon * sin(theta) + epsilon**2 / 4 * sin(2 theta) + O(epsilon**3)
        zeta = zeta0 + O(epsilon**5)'''
        nprofile_sa = 4

        with open('profile_sa.dat', 'r') as fprofile_sa:
            aminor = np.fromfile(
                fprofile_sa, dtype = float, count = 1, sep = ' ')
            nrad_sa, nprofile1_sa = np.fromfile(
                fprofile_sa, dtype = int, count = 2, sep = ' ')
            if self.input_.v >= 2:
                write('[load_sa] Info: aminor = {:g}\n'.format(aminor))
                write('[load_sa] Info: nrad_sa = {:d}\n'.format(nrad_sa))
                write('[load_sa] Info: nprofile_sa = {:d}\n'.format(
                    nprofile1_sa))
            if nprofile1_sa != nprofile_sa:
                raise ValueError(('[load_sa] Error: profile_sa.dat is not '
                    'compatible; nprofile_sa = {:d} (should be {:d}).'
                    ).format(nprofile1_sa, nprofile_sa))
            profiles_sa = np.loadtxt(fprofile_sa)
            if profiles_sa.shape != (nrad_sa, nprofile_sa):
                raise ValueError(('[load_sa] Error: profile_sa.dat is not '
                    'self-consistent; specified data dimension: ({:d}, {:d}); '
                    'actual data dimension: ({:d}, {:d}).'
                    ).format(nrad_sa, nprofile_sa,
                    profiles_sa.shape[0], profiles_sa.shape[1]))
        # End of with open('profile_sa.dat', 'r') as fprofile_sa:

        profiles_sa = profiles_sa.T

        # Interpolate profiles to output grid
        for iprofile in range(1, nprofile_sa):
            intp_profile = interp1d(
                profiles_sa[0, :], profiles_sa[iprofile, :], kind = 'cubic')
            self.profiles[iprofile, :] = intp_profile(self.profiles[0, :])

        # For debugging, simply copy raw profiles
        #self.profiles[0:4, :] = profiles_sa

        # Shift profiles[2:4, :] to profiles[3:5, :] for
        # the need to insert (g q + I) back to profiles[2, :]
        self.profiles[3:5, :] = self.profiles[2:4, :]
        # Add (g q + I) to profiles[2, :]
        # g = 1 - epsilon**2/2 + O(epsilon**3)
        # epsilon = r/R0 = acdprofile(1, :) * aminor
        # I = (epsilon**2 + epsilon**4/2)/q + O(epsilon**5)
        epsilons = self.profiles[0, :] * aminor
        self.profiles[2, :] = ((1.0 - epsilons**2 / 2.0) * self.profiles[1, :]
            + (epsilons**2 + epsilons**4 / 2.0) / self.profiles[1, :])

        #write('check profile data:\n')
        #for iprofile in range(self.nprofile):
            #write('[{:d}, 0:4]: {!s}\n'.format(
                #iprofile, self.profiles[iprofile, 0:4]))

        # Calculate fftcoefs, i.e., matrix elements in Eqs. (A.23)--(A.26) and
        #   (A.43) in [Nuclear Fusion 52, 043006 (2012)]
        #
        # Elements for matrix H, Eq. (A.23)
        # H = (epsilon/q)**2 / (g q + I)
        ifft = 0
        self.fftcoefs[0, ifft, :] = (
            (epsilons / self.profiles[1, :])**2 / self.profiles[2, :])

        # Elements for matrix J, Eq. (A.24)
        # J = (epsilon/q)**2 * (g q + I) / B^4
        # 1/B^4 = 1 + epsilon 4 cos(theta) + epsilon**2 * 2 * (5 cos(theta)**2 - 1 - 1/q**2) + O(epsilon**3)
        #       = (1 + epsilon**2 * (3 - 2 / q**2))
        #       + 2 epsilon (exp(i theta) + exp(-i theta))
        #       + 5/2 epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
        ifft = 1
        # First construct the part of (epsilon/q)**2 * (g q + I)
        for ifftcoef in range(min(2, self.input_.noffdiag) + 1):
            self.fftcoefs[ifftcoef, ifft, :] = (
                (epsilons / self.profiles[1, :])**2 * self.profiles[2, :])
        # Then multiply by 1/B^4 to finish constructing elements for J
        self.fftcoefs[0, ifft, :] *= (
            1.0 + epsilons**2 * (3.0 - 2.0 / self.profiles[1, :]**2))
        if self.input_.noffdiag >= 1:
            self.fftcoefs[1, ifft, :] *= 2.0 * epsilons
        if self.input_.noffdiag >= 2:
            self.fftcoefs[2, ifft, :] *= 2.5 * epsilons**2

        # Elements for matrix K, Eq. (A.25)
        # K = 2 g / B^3 * (partial B/partial theta)
        #   = 2 epsilon * sin(theta) + 3 epsilon**2 * sin(2 theta) + O(epsilon**3)
        #   = -i epsilon (exp(i theta) - exp(-i theta)) - i 3/2 epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
        ifft = 2
        if self.input_.noffdiag >= 1:
            self.fftcoefs[1, ifft, :] = -1j * epsilons
        if self.input_.noffdiag >= 2:
            self.fftcoefs[2, ifft, :] = -1.5j * epsilons**2

        # Elements for matrix L, Eq. (A.26)
        # 1 / B^2 = 1 + epsilon 2 cos(theta) + epsilon**2 (3 cos(theta)**2 - 1 - 1/q**2) + O(epsilon**3)
        #         = 1 + epsilon**2 (1/2 - 1/q**2)
        #         + epsilon (exp(i theta) + exp(-i theta))
        #         + 3/4 epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
        # 1 / B^4 = (1 + epsilon**2 * (3 - 2 / q**2))
        #         + 2 epsilon (exp(i theta) + exp(-i theta))
        #         + 5/2 epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
        # L = (g q + I) (beta# / B^4 + 1 / B^2)
        #   = (g q + I) * (
        #       beta# (1 + epsilon**2 * (3 - 2 / q**2)) + 1 + epsilon**2 (1/2 - 1/q**2)
        #     + (2 beta# + 1) epsilon (exp(i theta) + exp(-i theta))
        #     + (5/2 beta# + 3/4) epsilon**2 (exp(i 2 theta) + exp(-i 2 theta)) + O(epsilon**3)
        #     )
        ifft = 3
        for ifftcoef in range(min(2, self.input_.noffdiag) + 1):
            self.fftcoefs[ifftcoef, ifft, :] = self.profiles[2, :]
        self.fftcoefs[0, ifft, :] *= (
            self.profiles[3, :] 
            * (1.0 + epsilons**2 * (3.0 - 2.0 / self.profiles[1, :]**2))
            + 1.0 + epsilons**2 * (0.5 - 1.0 / self.profiles[1, :]**2))
        if self.input_.noffdiag >= 1:
            self.fftcoefs[1, ifft, :] *= epsilons * (
                2.0 * self.profiles[3, :] + 1.0)
        if self.input_.noffdiag >= 2:
            self.fftcoefs[2, ifft, :] *= epsilons**2 * (
                2.5 * self.profiles[3, :] + 0.75)

        # Elements for matrix N, Eq. (A.43)
        # N = 4 beta# g**2 / (g q + I) * (partial B/partial theta)**2 / ((beta# + B^2) B^2)
        #   = 4 beta# g**2 / (g q + I) * epsilon**2 * (
        #       sin(theta)**2 / (beta# + 1)
        #     + epsilon 2 (beta# + 2) sin(theta)**2 * cos(theta) / (beta# + 1)^2 + O(epsilon**2)
        #     )
        #   = beta# g**2 / (g q + I) * epsilon**2 * (
        #       2 / (beta# + 1)
        #     + epsilon (beta# + 2) / (beta# + 1)**2 * (exp(i theta) + exp(-i theta))
        #     - 1 / (beta# + 1) * (exp(i 2 theta) + exp(-i 2 theta))
        #     - epsilon (beta# + 2) / (beta# + 1)**2 * (exp(i 3 theta) + exp(-i 3 theta))
        #     )
        ifft = 4
        for ifftcoef in range(min(3, self.input_.noffdiag) + 1):
            self.fftcoefs[ifftcoef, ifft, :] = (self.profiles[3, :]
                * (1.0 - epsilons**2) / self.profiles[2, :] * epsilons**2)
        self.fftcoefs[0, ifft, :] *= 2.0 / (self.profiles[3, :] + 1.0)
        if self.input_.noffdiag >= 1:
            self.fftcoefs[1, ifft, :] *= (
                epsilons * (self.profiles[3, :] + 2.0)
                / (self.profiles[3, :] + 1.0)**2)
        if self.input_.noffdiag >= 2:
            self.fftcoefs[2, ifft, :] /= -(self.profiles[3, :] + 1.0)
        if self.input_.noffdiag >= 3:
            self.fftcoefs[3, ifft, :] *= (
                -epsilons * (self.profiles[3, :] + 2.0)
                / (self.profiles[3, :] + 1.0)**2)
        #write('check fft data:\n')
        #for ifft in range(self.nfft):
            #write('fftcoefs[0:4, {:d}, 0:4]: {!s}\n'.format(
                #ifft, self.fftcoefs[0:4, ifft, 0:4]))
    # End of def load_sa(self):

    def load_acd(self):
        '''Load an equilibrium from alcon.dat
        alcon.dat can be generated from GTC electromagnetic run with fload=0'''
        if self.input_.v >= 1:
            write('[load_acd] Info: loading equilibrium from alcon.dat...\n')
        with open('alcon.dat', 'r') as facd:
            nrad_acd, nfftcoef_acd, nprofile_acd, nfft_acd = np.fromfile(
                facd, dtype = int, count = 4, sep = ' ')
            if self.input_.v >= 2:
                write('[load_acd] Info: nrad_acd = {:d}\n'.format(nrad_acd))
                write('[load_acd] Info: nfftcoef_acd = {:d}\n'.format(
                    nfftcoef_acd))
                write('[load_acd] Info: nprofile_acd = {:d}\n'.format(
                    nprofile_acd))
                write('[load_acd] Info: nfft_acd = {:d}\n'.format(nfft_acd))
            if self.input_.noffdiag > nfftcoef_acd:
                write('[load_acd] Warning: specified noffdiag is larger '
                    'than nfftcoef in alcon.dat; not enough data to '
                    'calculate specified order of coupling.\n')
                write('You may want to re-generate your alcon.dat to '
                    'include higher m-harmonic amplitude data (increase '
                    '# of poloidal grid points for FFT and keep more '
                    'coefficients of the transformed array).\n')
            if nprofile_acd != self.nprofile or nfft_acd != self.nfft:
                raise ValueError(('[load_acd] Error: alcon.dat is not '
                    'compatible; nprofile = {:d} (should be {:d}); '
                    'nfft = {:d} (should be {:d})').format(
                    nprofile_acd, self.nprofile, nfft_acd, self.nfft))
            profiles_acd_raw = np.fromfile(facd, dtype = float,
                count = nprofile_acd * nrad_acd, sep = ' ')
            #fftcoefs_acd_raw = np.fromfile(facd, dtype = complex,
                #count = (nfftcoef_acd + 1) * nfft_acd * nrad_acd, sep = ' ')
            re_complex = re.compile(r'\(([^,\)]+),([^,\)]+)\)')
            fftcoefs_acd_raw = np.zeros(
                (nfftcoef_acd + 1) * nfft_acd * nrad_acd, dtype = complex)
            ifftcoef = 0
            for line in facd:
                matches = re_complex.findall(line)
                for match in matches:
                    if ifftcoef < (nfftcoef_acd + 1) * nfft_acd * nrad_acd:
                        fftcoefs_acd_raw[ifftcoef] = complex(
                            float(match[0]), float(match[1]))
                        ifftcoef += 1
                    else:
                        break
                if ifftcoef >= (nfftcoef_acd + 1) * nfft_acd * nrad_acd:
                    break
            if ifftcoef < (nfftcoef_acd + 1) * nfft_acd * nrad_acd:
                raise ValueError(('[load_acd] Error: expected {:d} complex '
                    'numbers for fftcoefs, but could only read {:d} from '
                    'alcon.dat').format(
                    (nfftcoef_acd + 1) * nfft_acd * nrad_acd, ifftcoef))
        # End of with open('alcon.dat', 'r') as facd:

        profiles_acd = profiles_acd_raw.reshape(
            (nrad_acd, nprofile_acd)).T
        fftcoefs_acd = fftcoefs_acd_raw.reshape(
            (nrad_acd, nfft_acd, nfftcoef_acd + 1)).T

        # Interpolate profiles and FFT coeficients to output grid
        for iprofile in range(1, self.nprofile):
            intp_profile = interp1d(
                profiles_acd[0, :], profiles_acd[iprofile, :], kind = 'cubic')
            self.profiles[iprofile, :] = intp_profile(self.profiles[0, :])
        for ifftcoef in range(min(nfftcoef_acd, self.input_.noffdiag) + 1):
            for ifft in range(self.nfft):
                intp_fftcoef = interp1d(profiles_acd[0, :],
                    fftcoefs_acd[ifftcoef, ifft, :], kind = 'cubic')
                self.fftcoefs[ifftcoef, ifft, :] \
                    = intp_fftcoef(self.profiles[0, :])

        # For debugging, take alcon.dat data directly without interpolation
        #self.profiles = profiles_acd
        #self.fftcoefs = fftcoefs_acd
    # End of def load_acd(self):
# End of class EqData:
