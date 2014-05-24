#!/usr/bin/env python
'''This file is the main program of the ALCON plotting utility.'''

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

__version__ = '2014-05-24 13:22:42-04:00'

import sys
write = sys.stdout.write
import alcon_plot_input
import alcon_plot_core

def main(argv = None):
    # Parse input parameters
    if argv is None:
        argv = sys.argv[1:]
    input_ = alcon_plot_input.parse(argv, __version__)

    # Initialize core
    core = alcon_plot_core.Core(input_)

    # Check if data files are not found
    if not core.fns:
        raise IOError(
            'Error: ALCON output data files are not found in {}.'.format(
                input_.dirout))

    # Specify default operation if no operation specified
    if (input_.matplotlib is None
        and input_.matlab is None
        and input_.idl is None
        and input_.gnuplot is None):
        write('Warning: no operation specified; '
            'default to "--matplotlib plot".\n')
        input_.matplotlib = 'plot'

    # Process matplotlib operation
    if input_.matplotlib is not None:
        core.matplotlib()

    # Process matlab operation
    if input_.matlab is not None:
        core.matlab()

    # Process idl operation
    if input_.idl is not None:
        core.idl()

    # Process gnuplot operation
    if input_.gnuplot is not None:
        core.gnuplot()
# End of def main():

if __name__ == '__main__':
    main()
