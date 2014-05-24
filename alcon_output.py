'''This file deals with output.'''

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
import os
import shutil

def prepare_dir(input_):
    '''Prepares output directory'''
    if input_.erasedirout:
        # Erase output directory
        if input_.v >= 1:
            write('[prepare_dir] Info: erasing output directory...\n')
        shutil.rmtree(input_.dirout, True)

    # Create output directory
    if input_.v >= 1:
        write('[prepare_dir] Info: creating output directory...\n')
    os.makedirs(input_.dirout)

def write_omegas(solver):
    '''Writes solved omegas to files'''
    if solver.input_.v >= 1:
        write('[write_omegas] Info: writing solutions to output files...\n')

    for impol in range(solver.nmpol):
        if solver.omegas[impol]:
            if impol < solver.input_.nmpol:
                mpol = impol + solver.input_.mpolrange[0]
                if solver.input_.finitebeta == 3:
                    fn_prefix = 'a_'
                else:
                    fn_prefix = ''
            else:
                mpol = (
                    impol - solver.input_.nmpol + solver.input_.mpolrange[0])
                fn_prefix = 's_'
            if mpol < 0:
                fn_m = 'm_{:02d}'.format(abs(mpol))
            else:
                fn_m = 'm{:03d}'.format(mpol)
            fn_omega = '{}n{:03d}{}.dat'.format(
                fn_prefix, solver.input_.ntor, fn_m)
            with open(
                os.path.join(solver.input_.dirout, fn_omega), 'w') as fomega:
                for rad_omega in solver.omegas[impol]:
                    fomega.write('{0[0]:g}\t{0[1]:g}\n'.format(rad_omega))
        # End of if solver.omegas[impol]:
    # End of for impol in range(solver.nmpol):
# End of def output(solver):
