'''This file implements the eigenvalue solver for Alfven continuum.'''

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
import numpy as np
import scipy.sparse as spsp
from scipy.sparse.linalg import eigs

class Solver:
    '''Solver for continuum'''

    def __init__(self, input_, eqdata):
        self.input_ = input_
        self.eqdata = eqdata
        self.nmpol = input_.nmpol
        if input_.finitebeta == 3:
            # For full continuum solution, the matrix dimension is doubled
            self.nmpol *= 2
        self.omegas = [[] for i in range(self.nmpol)] # List of solutions

    def _init_progress(self, nrad):
        '''Initializes variables for progress tracking'''
        self.nprogress = 10
        self.progress = (np.arange(1.0, self.nprogress + 1) * nrad
            / self.nprogress - 1.1)
        self.iprogress = 0

    def _prepare(self, irad):
        '''Prepares matrixes for radial grid irad'''
        matA = spsp.lil_matrix(
            (self.nmpol, self.nmpol), dtype = complex)
        matB = spsp.lil_matrix(
            (self.nmpol, self.nmpol), dtype = complex)

        for impol in range(self.input_.nmpol):
            jmpol1 = max(0, impol - self.input_.noffdiag)
            jmpol2 = min(self.input_.nmpol,
                impol + self.input_.noffdiag + 1)
            if self.input_.finitebeta == 2:
                kmpol = impol
                lmpol1 = jmpol1
                lmpol2 = jmpol2
            else:
                kmpol = impol + self.input_.nmpol
                lmpol1 = jmpol1 + self.input_.nmpol
                lmpol2 = jmpol2 + self.input_.nmpol

            if self.input_.finitebeta >= 2:
                # beta# * G^dagger * M * G part in matA
                matA[kmpol, kmpol] = (self.eqdata.profiles[3, irad]
                    * (self.input_.ntor * self.eqdata.profiles[1, irad]
                    - (self.input_.mpolrange[0] + impol))**2
                    / self.eqdata.profiles[2, irad])

                # L part in matB
                if impol > jmpol1:
                    matB[kmpol, lmpol1 : kmpol] = self.eqdata.fftcoefs[
                        impol - jmpol1 : 0 : -1, 3, irad].conjugate()
                matB[kmpol, kmpol : lmpol2] = self.eqdata.fftcoefs[
                    0 : jmpol2 - impol, 3, irad]

            if self.input_.finitebeta == 2:
                continue

            kmpol = impol + self.input_.nmpol
            lmpol1 = jmpol1 + self.input_.nmpol
            lmpol2 = jmpol2 + self.input_.nmpol

            # G^dagger * H * G part in matA
            if impol > jmpol1:
                matA[impol, jmpol1 : impol] = (self.eqdata.fftcoefs[
                    impol - jmpol1 : 0 : -1, 0, irad].conjugate()
                    * (self.input_.ntor * self.eqdata.profiles[1, irad]
                    - (self.input_.mpolrange[0] + impol))
                    * (self.input_.ntor * self.eqdata.profiles[1, irad]
                    - (self.input_.mpolrange[0] + np.arange(jmpol1, impol))))
            matA[impol, impol : jmpol2] = (self.eqdata.fftcoefs[
                0 : jmpol2 - impol, 0, irad]
                * (self.input_.ntor * self.eqdata.profiles[1, irad]
                - (self.input_.mpolrange[0] + impol))
                * (self.input_.ntor * self.eqdata.profiles[1, irad]
                - (self.input_.mpolrange[0] + np.arange(impol, jmpol2))))
            # Operator *= is not yet supported for slice of lil_matrix,
            # so the following snippet is commented out
            #matA[impol, jmpol1 : impol] = self.eqdata.fftcoefs[
                #jmpol1 - impol : 0, 0, irad].conjugate()
            #matA[impol, impol : jmpol2] = self.eqdata.fftcoefs[
                #0 : jmpol2 - impol, 0, irad]
            #matA[impol, jmpol1 : jmpol2] *= (
                #(self.input_.ntor * self.eqdata.profiles[1, irad]
                #- (self.input_.mpolrange[0] + impol))
                #* (self.input_.ntor * self.eqdata.profiles[1, irad]
                #- (self.input_.mpolrange[0] + np.arange(jmpol1, jmpol2))))

            # For slow sound approximation, A = G^dagger * H * G + N
            # Here adds the N part
            if self.input_.finitebeta == 1:
                if impol > jmpol1:
                    matA[impol, jmpol1 : impol] += self.eqdata.fftcoefs[
                        impol - jmpol1 : 0 : -1, 4, irad].conjugate()
                matA[impol, impol : jmpol2] += self.eqdata.fftcoefs[
                    0 : jmpol2 - impol, 4, irad]

            # J part in matB
            if impol > jmpol1:
                matB[impol, jmpol1 : impol] = self.eqdata.fftcoefs[
                    impol - jmpol1 : 0 : -1, 1, irad].conjugate()
            matB[impol, impol : jmpol2] = self.eqdata.fftcoefs[
                0 : jmpol2 - impol, 1, irad]

            if self.input_.finitebeta == 3:
                # -beta# * K part in matA
                if impol > jmpol1:
                    matA[impol, lmpol1 : kmpol] = (-self.eqdata.fftcoefs[
                        impol - jmpol1 : 0 : -1, 2, irad].conjugate()
                        * self.eqdata.profiles[3, irad])
                matA[impol, kmpol : lmpol2] = (-self.eqdata.fftcoefs[
                    0 : jmpol2 - impol, 2, irad]
                    * self.eqdata.profiles[3, irad])

                # K part in matB
                if impol > jmpol1:
                    matB[kmpol, jmpol1 : impol] = self.eqdata.fftcoefs[
                        impol - jmpol1 : 0 : -1, 2, irad].conjugate()
                matB[kmpol, impol : jmpol2] = self.eqdata.fftcoefs[
                    0 : jmpol2 - impol, 2, irad]
        # End of for impol in range(self.input_.nmpol):

        if self.input_.finitebeta == 3:
            # scipy.sparse.linalg.eigs() requires B to be Hermitian in the
            # eigen equation A * x = eigval * B * x, but our B is in general
            # not Hermitian for finitebeta == 3.  Noting that B^dagger * B is
            # always Hermitian, we multiply our eigen equation by B^dagger, and
            # define B_ = B^dagger * B and A_ = B^dagger * A.  Now our eigen
            # equation becomes: A_ * x = eigval * B_ * x, where B_ is
            # Hermitian.

            # B^dagger
            matBd = matB.conjugate().T

            # A_ = B^dagger * A
            # B_ = B^dagger * B, so that B_ is Hermitian
            matA_ = matBd * matA
            matB_ = matBd * matB

            # Store A_ and B_ as preparation for scipy.sparse.linalg.eigs()
            self._matA = matA_.tocsr()
            self._matB = matB_.tocsr()
            #write('{!s}\n'.format(self._matA))
            #write('{!s}\n'.format(self._matB))

        else: # if self.input_.finitebeta == 3:
            # In this case, B is Hermitian, so store A and B directly,
            # no need to use A_ and B_
            self._matA = matA.tocsr()
            self._matB = matB.tocsr()
        # End of else: if self.input_.finitebeta == 3:

    # End of def _prepare(self, irad):

    def _solve_radrange(self, iradrange, conn = None):
        '''Solves for given rage of radial grid irad'''
        for irad in range(*iradrange):
            # Print out progress on root process
            # (radial index range starts from 0 for root process)
            if iradrange[0] == 0:
                if self.input_.v == 1:
                    if irad == 0:
                        write('[solve] Info: solving continua...')
                        sys.stdout.flush()
                    while float(irad) >= self.progress[self.iprogress]:
                        self.iprogress += 1
                        write('{:.0%}...'.format(
                            float(self.iprogress) / self.nprogress))
                        sys.stdout.flush()
                        if self.iprogress == self.nprogress:
                            write('\n')
                            break
                elif self.input_.v >= 2:
                    write(('[_solve_radrange] Info: solving radial grid: '
                        '{:d} / {:d}\n').format(irad, self.input_.nrad))

            # Prepare matrixes A and B
            self._prepare(irad)

            # Adjust ncon if neccessary
            ncon = max(1, min(self.input_.ncon, self.nmpol - 1))

            # Solve eigenvalue problem: A * x = eigval * B * x
            eigvals, eigvecs = eigs(self._matA, k = ncon,
                M = self._matB, sigma = self.input_.sigma)
            if self.input_.v >= 2:
                write(('[_solve_radrange][{:d}] '
                    'Info: eigenvalues: {!s}\n').format(irad, eigvals))
            if self.input_.v >= 3:
                write(('[_solve_radrange][{:d}] '
                    'Info: eigenvectors: {!s}\n').format(irad, eigvecs))

            # Loop over eigen-solutions and store useful ones
            for ieig in range(len(eigvals)):
                omega = np.sqrt(eigvals[ieig] / self.eqdata.profiles[4, irad])

                # Skip solution whose imaginary part is too large
                if abs(omega.imag) >= omega.real * self.input_.imreratiocutoff:
                    continue

                # Scale solution
                omega_scaled = omega.real * self.input_.omegascale

                # Skip solution if scaled solution is too large
                if (self.input_.omegacutoff > 0
                    and self.input_.omegacutoff < omega_scaled):
                    continue

                ieigvecabsmax = np.argmax(np.abs(eigvecs[:, ieig]))
                self.omegas[ieigvecabsmax].append(
                    (self.eqdata.profiles[0, irad], omega_scaled))
            # End of for ieig in range(len(eigvals)):
            if self.input_.v >= 2:
                write('[_solve_radrange][{:d}] Info: omegas: {!s}\n'.format(
                    irad, self.omegas))
        # End of for irad in range(*iradrange):

        # Send results back to root process through given connection
        if conn is not None:
            conn.send(self.omegas)
    # End of def _solve_radrange(self, iradrange, conn = None):

    def solve(self):
        if self.input_.n > 1:
            # Use multiprocessing module for parallelization
            import multiprocessing as mp

            # Radial ranges for processes
            radranges = (self.input_.nrad * np.arange(self.input_.n + 1)
                // self.input_.n)

            # Process objects for subprocesses
            pes = [None] * self.input_.n

            # Connectors for Pipes between each subprocess and root process
            conns = [None] * self.input_.n

            # Loop for creating subprocesses
            for ipe in range(1, self.input_.n):
                # Create Pipe for transferring solutions
                conns[ipe], conn_sub = mp.Pipe()

                # Create subprocess
                pes[ipe] = mp.Process(target = self._solve_radrange,
                    args = ((radranges[ipe], radranges[ipe + 1]), conn_sub))
                pes[ipe].start()

            # Root process's solving job
            self._init_progress(radranges[1])
            self._solve_radrange((radranges[0], radranges[1]))

            # Loop for collecting solutions from subprocesses
            for ipe in range(1, self.input_.n):
                # Receive solutions from subprocess
                omegas = conns[ipe].recv()

                # Merge subprocess's solutions to root process's
                for impol in range(self.nmpol):
                    self.omegas[impol] += omegas[impol]

                pes[ipe].join() # Clean subprocess
        else: # Of if self.input_.n > 1:
            self._init_progress(self.input_.nrad)
            self._solve_radrange((0, self.input_.nrad))
    # End of def solve(self):
# End of class Solver:
