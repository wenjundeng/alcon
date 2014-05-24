'''This file includes a wall-clock timer class.'''

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
import time

class WCTimer:
    '''Wall-clock timer for timing each part of the code'''

    def __init__(self, ntimer):
        self.ntimer = ntimer
        self.times_total = [0.0] * ntimer
        self.times_start = [0.0] * ntimer
        self.times_current = [0.0] * ntimer
        self.brunning = [False] * ntimer
        self.names = [''] * ntimer

    def name(self, itimer, name):
        '''Names a timer'''
        if len(name) > 18:
            write(('[WCTimer] Warning: length of name "{}" is larger than 18; '
                'print_summary() may give awkward formatting.\n').format(name))
        self.names[itimer] = name

    def start(self, itimer):
        '''Starts a timer'''
        self.times_start[itimer] = time.time()
        self.brunning[itimer] = True

    def pause(self, itimer):
        '''Pauses a timer and add duration to total time'''
        self.times_total[itimer] += time.time() - self.times_start[itimer]
        self.brunning[itimer] = False

    def _update_current(self):
        '''Updates times_current'''
        t = time.time()
        for itimer in range(self.ntimer):
            if self.brunning[itimer]:
                self.times_current[itimer] = (self.times_total[itimer]
                    + t - self.times_start[itimer])
            else:
                self.times_current[itimer] = self.times_total[itimer]

    def print_summary(self, itimer_percent = None):
        '''Prints out a summary of all timers'''
        self._update_current()
        write('[WCTimer] Summary:\n')
        nline = (self.ntimer + 3) // 4
        for iline in range(nline):
            for jtimer in range(4):
                itimer = iline * 4 + jtimer
                if itimer >= self.ntimer:
                    break
                write('{:>18s} '.format(self.names[itimer]))
            write('\n')
            for jtimer in range(4):
                itimer = iline * 4 + jtimer
                if itimer >= self.ntimer:
                    break
                if itimer_percent is None:
                    write('{:>18s} '.format(
                        self.sec2text(self.times_current[itimer])))
                else:
                    write('{:>9s} ({:6.1%}) '.format(
                        self.sec2text(self.times_current[itimer]),
                        self.times_current[itimer]
                        / self.times_current[itimer_percent]))
            write('\n')
    # End of def print_summary(self, itimer_percent = None):

    def sec2text(self, sec):
        '''Converts # of seconds to human readable text
        e.g., 5.2min, 3.2hr, 4.1day'''
        suffix = 'sec'
        t = sec
        if t > 60.0:
            t /= 60.0
            suffix = 'min'
            if t > 60.0:
                t /= 60.0
                suffix = 'hr'
                if t > 24.0:
                    t /= 24.0
                    suffix = 'day'
        # End of if t > 60.0:
        return '{:.1f}{}'.format(t, suffix)

# End of class WCTimer:
