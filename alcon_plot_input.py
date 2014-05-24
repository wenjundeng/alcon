'''This file deals with input parameters for alcon_plot.'''

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
        description = 'ALCON plotting utility version {}'.format(version))

    parser.add_argument('--dirout', type = str,
        metavar = '<directory of ALCON output>', default = './output',
        help = 'default: %(default)s')

    # Operations
    parser.add_argument('--matplotlib', nargs = '?', type = str,
        choices = ['plot', 'png', 'pdf'], const = 'plot',
        help = ('use matplotlib to generate a figure; '
            'plot: plot on screen; png: generate a png figure file; pdf: '
            'generate a pdf figure file; default: %(const)s'))
    parser.add_argument('--matlab', nargs = '?', type = str,
        choices = ['m'], const = 'm',
        help = ('generate a MATLAB script for plotting; '
            'm: generate a .m script; default: %(const)s'))
    parser.add_argument('--idl', nargs = '?', type = str,
        choices = ['pro'], const = 'pro',
        help = ('generate an IDL script for plotting; '
            'pro: generate a .pro script; default: %(const)s'))
    parser.add_argument('--gnuplot', nargs = '?', type = str,
        choices = ['eps', 'epso', 'mp', 'mpo'], const = 'eps',
        help = ('generate a gnuplot script for plotting; '
            'eps: generate a .gp script and '
            'a Makefile for making a .eps figure; epso: same as the '
            'operation of eps, except that the Makefile will make use of '
            'Omega Tools <http://wdeng.info/codes/omegatools/>; '
            'mp: same as the operation of eps, except that the .gp script '
            'will use metapost termianl for output (a modern LaTeX '
            'distribution is required), then "make" generates a .mps '
            'figure, "make alcon.eps" generates a .eps figure (requires '
            'mps2eps <http://www.ctan.org/tex-archive/support/mps2eps> '
            'or Omega Tools), "make alcon.pdf" generates a .pdf figure, '
            '"make alcon.png" generates a .png figure (requires convert '
            'in ImageMagick <http://www.imagemagick.org>; mpo: same as '
            'the operation of mp, except that the Makefile will make use '
            'of Omega Tools and make all labels in LaTeX math mode; '
            'default: %(const)s'))

    # Output parameters
    colors_default = [
        '#BFBF00', '#CCCCCC', '#000000', '#FF0000', '#3333FF',
        '#FFAA00', '#007300', '#E600E6', '#00D900', '#00D9D9',
        '#990000', '#66B366', '#000080', '#990099', '#0080FF']
    parser.add_argument('--colors', nargs = '+', type = str,
        metavar = '<hex color>', default = colors_default,
        help = ('provide colors in hex RGB format, e.g., #FF0000 '
            'indicating red; the 1st and 2nd colors are for unlabeled '
            'Alfven and sound continua, respectively (if solving only '
            'the sound continua, the 1st color will be used for '
            'unlabeled sound continua); the 3rd and '
            'subsequent colors are for labeled continua; '
            'default: {}').format(' '.join(colors_default)))
    parser.add_argument('--mpol-labeled-start', type = int,
        metavar = '<poloidal m number from which labeled continua start>',
        default = 5, help = 'default: %(default)s')
    parser.add_argument('--radrange', type = float, nargs = 2,
        metavar = ('<rad1>', '<rad2>'),
        help = ('plotting radial coordinate (usually rho) range from rad1 '
            'to rad2; for more about radial coordinate, run '
            '"python alcon.py -h" '
            'and refer to the option --radrange; default: automatic'))
    parser.add_argument('--omegarange', type = float, nargs = 2,
        metavar = ('<omega1>', '<omega2>'),
        help = ('plotting omega range from omega1 '
            'to omega2; default: automatic'))
    parser.add_argument('--radlabel', type = str,
        metavar = '<label for radial coordinate, i.e., x-axis>')
    parser.add_argument('--omegalabel', type = str,
        metavar = '<label for omega, i.e., y-axis>',
        help = ('note that these labels may be interpreted differently for '
            'different plotting operations; for example, when using '
            '"--matplotlib", LaTeX math expressions can be used within '
            'a pair of "$" signs; when using "--gnuplot mpo", these '
            'labels are completely in LaTeX math mode, so use \\textrm{} '
            'to type normal (non-math) text'))

    # Other parameters
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

    # If only 1 color specified, duplicate it as at least 2 colors are
    # needed in the list.
    if len(args.colors) == 1:
        args.colors *= 2

    if args.v >= 1:
        write('ALCON plotting utility version {}\n'.format(version))
        write('[input.parse] Info: input parameters: {!s}\n'.format(args))
    return args
# End of def parse(version):
