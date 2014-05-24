#!/usr/bin/env python
'''This file is the main program of ALCON.'''

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

import wctimer

import alcon_input
import alcon_eqdata
import alcon_solver
import alcon_output

def main(argv = None):
    # Create timer object
    wct = wctimer.WCTimer(5)

    # Start the timer for total time
    wct.start(0)
    wct.name(0, 'Total')

    # Parse input parameters
    wct.name(1, 'Parse input')
    wct.start(1)
    if argv is None:
        argv = sys.argv[1:]
    input_ = alcon_input.parse(argv, __version__)
    wct.pause(1)

    # Load equilibrium data
    wct.name(2, 'Load eq. data')
    wct.start(2)
    eqdata = alcon_eqdata.EqData(input_)
    wct.pause(2)

    # Prepare output directory
    wct.name(4, 'Output')
    wct.start(4)
    alcon_output.prepare_dir(input_)
    wct.pause(4)

    # Initialize solver and solve continua
    wct.name(3, 'Solver')
    wct.start(3)
    solver = alcon_solver.Solver(input_, eqdata)
    solver.solve()
    wct.pause(3)

    # Write output data to files
    wct.start(4)
    alcon_output.write_omegas(solver)
    wct.pause(4)

    # Stop total timer and print out timer summary
    wct.pause(0)
    if input_.v >= 1:
        wct.print_summary(0)
# End of def main():

if __name__ == '__main__':
    main()
