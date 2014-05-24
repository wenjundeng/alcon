'''This file deals with input parameters.'''

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
import argparse

def parse(argv, version):
    '''Parses command line arguments as input parameters'''
    parser = argparse.ArgumentParser(
        description = 'ALCON (Alfven continuum solver) version {}'.format(
            version))

    # Equilibrium type
    parser.add_argument('--eqtype', type = str,
        choices = ['alcon.dat', 'sa'], default = 'alcon.dat',
        help = ('equilibrium type: alcon.dat: load numerical equilibrium '
            'data from alcon.dat file; sa: use a simple analytic '
            'equilibrium model described in Appendix E in '
            '[Nuclear Fusion 52, 023005 (2012)](http://wdeng.info/?p=76) '
            '(extended to second order), and load profiles from '
            'profile_sa.dat; '
            #'sppro: load numerical equilibrium from spdata.dat and '
            #'profile.dat (not implemented yet); '
            'default: %(default)s'))

    # Solver parameters
    parser.add_argument('--finitebeta', type = int,
        metavar = '<finite beta type>', default = 3,
        help = ('0: pure Alfven continua, no acoustic coupling; 1: Alfven '
            'continua with acoustic coupling with slow sound approximation '
            '[Physics of Fluids B: Plasma Physics 4, 3713 (1992)]; '
            '2: pure sound continua; 3: accurate Alfven continua with '
            'acoustic coupling; default: %(default)s'))
    parser.add_argument('--radrange', type = float, nargs = 2,
        metavar = ('<rad1>', '<rad2>'), default = [0.01, 1.0],
        help = ('solving radial coordinate (usually rho) range from rad1 '
            'to rad2; rho is square root of normalized toroidal '
            'flux: rho = sqrt(psi_tor / psi_tor_wall); if eqtype == 1, '
            'the radial coordinate is determined by alcon.dat; if '
            'eqtype == 100, the radial coordinate is rho, i.e., r/a '
            'in this simple equilibrium; if eqtype == 2, the radial '
            'coordinate is determined by profile.dat; note that if '
            'rad1 == 0, on the first radial grid, all matrix elements '
            'are 0, so there would be no output for the first radial '
            'grid; default: 0.01 1.0'))
    parser.add_argument('--nrad', type = int,
        metavar = '<# of radial grid points>', default = 2000,
        help = 'default: %(default)s')
    parser.add_argument('--ntor', type = int,
        metavar = '<toroidal mode number n>', default = 4,
        help = 'default: %(default)s')
    parser.add_argument('--mpolrange', type = int, nargs = 2,
        metavar = ('<m1>', '<m2>'), default = [-20, 50],
        help = 'poloidal mode number m range from m1 to m2; default: -20 50')
    parser.add_argument('--noffdiag', type = int,
        metavar = '<# of off-diagonal stripes>', default = 20,
        help = ('representing cut-off of m-harmonic coupling; 0: no '
            'coupling, trivial solution: \omega = v_A k_\parallel; 1: '
            'm-harmonic only couples with (m - 1) and (m + 1), TAE gap '
            'shows; 2: m-harmonic only couples with (m - 1), (m + 1), '
            '(m - 2), and (m + 2); default: %(default)s'))
    parser.add_argument('--ncon', type = int,
        metavar = '<# of lowest Alfven continua to solve>', default = 50,
        help = 'default: %(default)s')
    parser.add_argument('--imreratiocutoff', type = float,
        metavar = '<imaginary to real ratio cutoff for a solution>',
        default = 0.1, help = ('drop solution if the absolute value of its '
            'imaginary part to its real part is larger than this cutoff '
            'value due to too large error; because the Alfven continuuma '
            'eigenvalue equation is positive semi-definite, the imaginary '
            'part of the eigenvalue should be 0; default: %(default)s'))
    parser.add_argument('--sigma', type = float,
        metavar = '<sigma for shift-invert mode>', default = -0.5,
        help = ('this is for using the shift-invert mode to assist the '
            'eigenvalue solver to find the lowest continua; the '
            'default value is recommended unless ALCON has difficulty '
            'in converging the solution to the accurate value, '
            'in which case you can try to adjust this sigma; '
            'default: %(default)s'))

    # Output parameters
    parser.add_argument('--omegascale', type = float,
        metavar = '<scale factor for omega before output>', default = 1.0,
        help = ('1.0 for no scaling, in which case omega is normalized by '
            'v_{Ap}/R_0, i.e., Eq. (A.36) in [Nuclear Fusion 52, 043006 '
            '(2012)](http://wdeng.info/?p=117); default: %(default)s'))
    parser.add_argument('--omegacutoff', type = float,
        metavar = '<cutoff value for scaled omega>', default = -1.0,
        help = ('drop solution if the scaled omega is larger than this '
            'cutoff value; set to negative to disable cutting-off; '
            'default: %(default)s'))
    parser.add_argument('--dirout', type = str,
        metavar = '<directory for output>', default = './output',
        help = ('error will occur if <directory for output> preexists, '
            'in which case remove it first or '
            'use --erasedirout; default: %(default)s'))
    parser.add_argument('--erasedirout', action = 'store_true',
        help = 'erase output directory if it preexists')

    # Other parameters
    parser.add_argument('-n', type = int,
        metavar = '<# of processes>', default = 1,
        help = ('parallelize the solver by the specified # of processes; '
            'if specified as 1, no parallelization; default: %(default)s'))
    parser.add_argument('-v', '--verbose',
        action = 'append_const', const = 1, default = [1],
        help = ('display more information; can be used multiple times to '
            'raise verbosity level, e.g., -vv'))
    parser.add_argument('-q', '--quiet',
        action = 'append_const', const = -1, dest = 'verbose',
        help ='display less information')
    parser.add_argument('--version', action = 'version',
        version = '%(prog)s {}'.format(version))

    args = parser.parse_args(argv)

    # Calculate verbosity
    args.v = 0 if args.verbose is None else sum(args.verbose)

    # Calculate # of poloidal m numbers
    args.nmpol = args.mpolrange[1] - args.mpolrange[0] + 1

    if args.v >= 1:
        write('ALCON version {}\n'.format(version))
        write('[input.parse] Info: input parameters: {!s}\n'.format(args))
    return args
# End of def parse(version):
