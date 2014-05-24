'''This file provides core funtionalities for alcon_plot.'''

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
import glob
import re
import numpy as np

class Core:
    '''Core functionalities for alcon_plot'''

    def __init__(self, input_):
        self.input_ = input_

        # Regular expression for analyzing file names
        self._re_fn = re.compile(r'(.*)n[0-9]{3}m([_0-9][0-9]{2}).dat')

        # Get ALCON output file names
        self.fns = glob.glob(os.path.join(input_.dirout,
            '*n[0-9][0-9][0-9]m[_0-9][0-9][0-9].dat'))

        # List of analysis of file names
        self._afns = [self._analyze_fn(fn) for fn in self.fns]

        # Indexes of file names sorted by self._analyze_fn()
        self._is_fns = sorted(range(len(self.fns)),
            key = lambda i: self._afns[i])

    def _analyze_fn(self, fn):
        '''Analyzes filename given by fn, return whether it should be plotted in
        thick line and its poloidal m number'''
        fn_base = os.path.basename(fn)
        match = self._re_fn.match(fn_base)
        bthick = match.group(1) != 's_'
        if match.group(2)[0] == '_':
            mpol = -int(match.group(2)[1:])
        else:
            mpol = int(match.group(2))
        return (bthick, mpol)

    def matplotlib(self):
        '''Processes matplotlib operation'''
        import matplotlib.pyplot as plt

        if self.input_.v >= 1:
            write('[matplotlib] processing plot...\n')

        # Plot
        for ifn in self._is_fns:
            fn = self.fns[ifn]
            bthick = self._afns[ifn][0]
            mpol = self._afns[ifn][1]
            icolor = mpol - self.input_.mpol_labeled_start + 2
            if icolor < 2 or icolor >= len(self.input_.colors):
                color = (self.input_.colors[0] if bthick
                    else self.input_.colors[1])
                label = None
            else:
                color = self.input_.colors[icolor]
                if bthick:
                    label = '$m = {:d}$'.format(mpol)
                else:
                    label = None
            size = (5.0 if bthick else 1.0)
            data = np.loadtxt(fn)
            if len(data.shape) == 1:
                data = data.reshape((1, 2))

            ax = plt.subplot(111)
            if label is None:
                ax.plot(data[:, 0], data[:, 1], 'o',
                    markeredgecolor = color, markerfacecolor = color,
                    markersize = size)
            else:
                ax.plot(data[:, 0], data[:, 1], 'o',
                    markeredgecolor = color, markerfacecolor = color,
                    markersize = size, label = label)
        # End of for ifn in self._is_fns:

        # Set axis ranges and labels
        if self.input_.radrange is not None:
            ax.set_xlim(self.input_.radrange)
        if self.input_.omegarange is not None:
            ax.set_ylim(self.input_.omegarange)
        if self.input_.radlabel is not None:
            ax.set_xlabel(self.input_.radlabel)
        if self.input_.omegalabel is not None:
            ax.set_ylabel(self.input_.omegalabel)

        # Set legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        legend = ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),
            numpoints = 20)

        # Set legend font color
        for i, text in enumerate(legend.get_texts()):
            text.set_color(self.input_.colors[i + 2])

        # Save figure to file or plot to screen
        if (self.input_.matplotlib == 'pdf'
            or self.input_.matplotlib == 'png'):
            plt.savefig(os.path.join(
                self.input_.dirout, 'alcon.{}'.format(self.input_.matplotlib)),
                bbox_inches = 'tight')
        elif self.input_.matplotlib == 'plot':
            plt.show()
    # End of def matplotlib(self):

    def matlab(self):
        '''Processes MATLAB operation'''

        if self.input_.v >= 1:
            write('[matlab] generating script...\n')

        # Define colors
        colors_m = [str([int(sc[1 : 3], 16) / 255.0,
            int(sc[3 : 5], 16) / 255.0,
            int(sc[5 : 7], 16) / 255.0])
            for sc in self.input_.colors]
        sm = 'colors = [{}];\n\n'.format(';\n'.join(colors_m))

        # Load data files
        for ifn in self._is_fns:
            fn_base = os.path.basename(self.fns[ifn])
            sm += "load '{}'\n".format(fn_base)

        # Plot
        sm += 'figure\nhold on\n\n'
        slegend = 'legend(hm, ...\n'
        for ifn in self._is_fns:
            fn_base2 = os.path.basename(self.fns[ifn])[: -4]
            bthick = self._afns[ifn][0]
            mpol = self._afns[ifn][1]
            icolor = mpol - self.input_.mpol_labeled_start + 3
            if icolor < 3 or icolor > len(self.input_.colors):
                if bthick:
                    icolor = 1
                else:
                    icolor = 2
                hm = ''
            else:
                if bthick:
                    hm = 'hm({:d}) = '.format(icolor - 2)
                    slegend += "'m = {:d}', ...\n".format(mpol)
                else:
                    hm = ''
            if bthick:
                style = 'filled'
                size = 20
            else:
                style = '.'
                size = 1
            sm += ("{0}scatter({1}(:, 1), {1}(:, 2), {2:d}, "
                "colors({3:d}, :), '{4}');\n").format(
                hm, fn_base2, size, icolor, style)
        # End of for ifn in self._is_fns:
        slegend += "'Location', 'NorthEastOutside');\n"
        sm += slegend

        # Set axis ranges and labels
        if self.input_.radrange is not None:
            sm += 'xlim([{0[0]:g} {0[1]:g}]);\n'.format(self.input_.radrange)
        if self.input_.omegarange is not None:
            sm += 'ylim([{0[0]:g} {0[1]:g}]);\n'.format(self.input_.omegarange)
        if self.input_.radlabel is not None:
            sm += "xlabel('{}');\n".format(self.input_.radlabel)
        if self.input_.omegalabel is not None:
            sm += "ylabel('{}');\n".format(self.input_.omegalabel)

        with open(os.path.join(self.input_.dirout, 'alcon.m'), 'w') as fm:
            fm.write(sm)
    # End of def matlab(self):

    def idl(self):
        '''Processes IDL operation'''

        if self.input_.v >= 1:
            write('[idl] generating script...\n')

        spro = 'pro alcon, eps = eps\n\n'

        # Define colors
        colors = ['#FFFFFF', '#000000'] + self.input_.colors
        ncolor = len(colors)
        colors_pro = [
            '{:d}l'.format(int(sc[5 : 7] + sc[3 : 5] + sc[1 : 3], 16))
            for sc in colors]
        spro += 'colors = [{}]\n\n'.format(', '.join(colors_pro))

        # Define eps colors
        spro += ('if keyword_set(eps) then begin\n'
            '\told_device = !d.name\n'
            "\tset_plot, 'PS'\n"
            "\tdevice, file = 'alcon.eps', /encapsulated, /color, bits=8\n")
        for ic, sc in enumerate(colors):
            spro += '\ttvlct, {:d}, {:d}, {:d}, {:d}\n'.format(
                int(sc[1 : 3], 16), int(sc[3 : 5], 16), int(sc[5 : 7], 16), ic)
        spro += '\tcolortable = lindgen({:d})\nendif\n\n'.format(ncolor)

        # Load data files
        spro += 'lun = 124\n'
        if self.input_.omegarange is None:
            spro += 'arromegamax = dblarr({:d})\n'.format(len(self.fns))
        for iifn, ifn in enumerate(self._is_fns):
            fn_base = os.path.basename(self.fns[ifn])
            mpol = self._afns[ifn][1]
            spro += ("nlines = file_lines('{0}', /noexpand_path)\n"
                'data{1:d} = dblarr(2, nlines)\n'
                "openr, lun, '{0}'\n"
                'readf, lun, data{1:d}\n'
                'close, lun\n').format(fn_base, iifn)
            if self.input_.omegarange is None:
                spro += 'arromegamax({0:d}) = max(data{0:d}(1, *))\n'.format(
                    iifn)

        # Plot
        spro += ('\nt = findgen(17) * (2.0 * !pi / 16.0)\n'
            'usersym, 0.3 * cos(t), 0.3 * sin(t), /fill\n\n'
            '!p.psym = 8\n'
            '!p.symsize = 1.0\n'
            '!p.background = colors(0)\n\n')
        if self.input_.radrange is None:
            spro += 'xrange = [0.0, 1.0]\n'
        else:
            spro += 'xrange = [{0[0]:g}, {0[1]:g}]\n'.format(
                self.input_.radrange)
        if self.input_.omegarange is None:
            spro += ('omegamax = max(arromegamax)\n'
                'yrange = [0.0, omegamax]\n')
        else:
            spro += 'yrange = [{0[0]:g}, {0[1]:g}]\n'.format(
                self.input_.omegarange)

        # Plot the axes
        spro += ('plot, [-100.0], [-100.0], '
            'color = colors(1), xrange = xrange, yrange = yrange')
        if self.input_.radlabel:
            spro += ", xtitle = '{}'".format(self.input_.radlabel)
        if self.input_.omegalabel:
            spro += ", ytitle = '{}'".format(self.input_.omegalabel)
        spro += '\n'
        for iifn, ifn in enumerate(self._is_fns):
            fn_base = os.path.basename(self.fns[ifn])
            bthick = self._afns[ifn][0]
            mpol = self._afns[ifn][1]
            icolor = mpol - self.input_.mpol_labeled_start + 3
            if icolor < 3 or icolor >= len(colors):
                icolor = (1 if bthick else 2)

            ssymsize = (', symsize = 3.0' if bthick else '')
            spro += ('oplot, data{0:d}(0, *), data{0:d}(1, *), '
                'color = colors({1:d}){2}\n').format(iifn, icolor, ssymsize)

        # Labels
        spro += '\n; put label at the data average position\n'
        for iifn, ifn in enumerate(self._is_fns):
            bthick = self._afns[ifn][0]
            mpol = self._afns[ifn][1]
            icolor = mpol - self.input_.mpol_labeled_start + 3
            if icolor >= 3 and icolor < len(colors) and bthick:
                spro += (
                    "xyouts, mean(data{0:d}(0, *)), mean(data{0:d}(1, *)), "
                    "'m={1:d}', color = colors({2:d})\n").format(
                    iifn, mpol, icolor)

        spro += ('\nif keyword_set(eps) then begin\n'
            '\tdevice, /close\n'
            '\tset_plot, old_device\n'
            'endif\n\nend\n')

        with open(os.path.join(self.input_.dirout, 'alcon.pro'), 'w') as fpro:
            fpro.write(spro)
    # End of def idl(self):

    def gnuplot(self):
        '''Processes gnuplot operation'''

        if self.input_.v >= 1:
            write('[gnuplot] generating script...\n')

        # Set terminal and line styles
        if self.input_.gnuplot[:3] == 'eps':
            sgp = 'set terminal postscript enhanced color eps\n'
            sgp += "set output 'alcon.eps'\n"
            ps1 = 0.6
            ps2 = 0.2
            pt = 7
            offset_xtics = 0.4
            offset_ytics = 0.5
            offset_xlabel = 1.0
            offset_ylabel = 2.0
        else: # Implying self.input_.gnuplot == 'mp' or 'mpo'
            sgp = 'set terminal mp color\n'
            sgp += "set output 'alcon.mp'\n"
            ps1 = 0.6
            ps2 = 0.001
            pt = 9
            offset_xtics = 0.5
            offset_ytics = 0.5
            offset_xlabel = 1.1
            offset_ylabel = 2.5
        for ils, color in enumerate(
            self.input_.colors[0:1] + self.input_.colors[2:]):
            sgp += ("set style line {:d} "
                "lc rgbcolor '{}' pt {:d} ps {:f}\n").format(
                100 + ils, color, pt, ps1)
        for ils, color in enumerate(self.input_.colors[1:]):
            sgp += ("set style line {:d} "
                "lc rgbcolor '{}' pt {:d} ps {:f}\n").format(
                200 + ils, color, pt, ps2)

        # Set tics, labels, ranges, key
        sgp += ('set xtics offset 0, {:f}\n'
            'set ytics offset {:f}, 0\n').format(
            offset_xtics, offset_ytics)
        if self.input_.radlabel:
            sgp += "set xlabel '{}' offset 0, {:f}\n".format(
                self.input_.radlabel, offset_xlabel)
        if self.input_.omegalabel:
            sgp += "set ylabel '{}' offset {:f}, 0\n".format(
                self.input_.omegalabel, offset_ylabel)
        if self.input_.radrange is not None:
            sgp += 'set xrange [{0[0]:g} : {0[1]:g}]\n'.format(
                self.input_.radrange)
        if self.input_.omegarange is not None:
            sgp += 'set yrange [{0[0]:g} : {0[1]:g}]\n'.format(
                self.input_.omegarange)
        sgp += ('set key outside\n'
            'set key textcolor variable\n'
            'plot \\\n')

        # Plot
        splots = []
        for ifn in self._is_fns:
            fn_base = os.path.basename(self.fns[ifn])
            bthick = self._afns[ifn][0]
            mpol = self._afns[ifn][1]
            icolor = mpol - self.input_.mpol_labeled_start + 1
            if icolor < 1 or icolor >= len(self.input_.colors) - 1:
                icolor = (100 if bthick else 200)
                title = 'notitle'
            else:
                if bthick:
                    icolor += 100
                    title = "title 'm = {:d}'".format(mpol)
                else:
                    icolor += 200
                    title = 'notitle'

            splots += ["'{}' {} with points ls {:d}".format(
                fn_base, title, icolor)]
        # End of for ifn in self._is_fns:
        sgp += ', \\\n'.join(splots) + '\n'

        with open(os.path.join(self.input_.dirout, 'alcon.gp'), 'w') as fgp:
            fgp.write(sgp)

        # Generate Makefile
        bomegatools = (self.input_.gnuplot[-1] == 'o')
        if self.input_.gnuplot[:3] == 'eps':
            smf = '.PHONY : clean\n\n'
        else: # Implying self.input_.gnuplot == 'mp' or 'mpo'
            if bomegatools:
                smps = '\tmpost_mps $<\n'
                spdf = '\tmptopdf_mps $<\n'
            else:
                smps = ('\tmpost $<\n'
                    '\tmv alcon.0 $@\n')
                spdf = ('\tmptopdf $<\n'
                    '\tmv alcon-mps.pdf $@\n')
            smf = ('.PHONY : clean cleanall\n\n'
                'alcon.mps : alcon.mp\n'
                '{}'
                '\nalcon.eps : alcon.mps\n'
                '\tmps2eps $<\n\n'
                'alcon.pdf : alcon.mps\n'
                '{}'
                '\nalcon.png : alcon.eps\n'
                '\tconvert -density 400 -flatten $< $@\n\n').format(
                    smps, spdf)
        if bomegatools:
            if self.input_.gnuplot[:3] == 'eps':
                smf += 'alcon.eps : alcon.gp alcon.d\n'
            else:
                smf += 'alcon.mp : alcon.gp alcon.d\n'
        else: # Implying not to use Omega Tools
            if self.input_.gnuplot[:3] == 'eps':
                smf += 'alcon.eps : alcon.gp \\\n'
            else:
                smf += 'alcon.mp : alcon.gp \\\n'
            sfns = [os.path.basename(self.fns[ifn]) for ifn in self._is_fns]
            smf += ' {}\n'.format(' \\\n '.join(sfns))

        if self.input_.gnuplot == 'mpo':
            smf += '\tgnuplot_lmp $<\n\n'
        else:
            smf += '\tgnuplot $<\n\n'
        if bomegatools:
            smf += ('alcon.d : alcon.gp\n'
                '\tmakedependgp $<\n\n'
                '-include alcon.d\n\n')
        smf += 'clean :\n'
        sd = (' alcon.d' if bomegatools else '')
        if self.input_.gnuplot[:3] == 'eps':
            smf += '\trm -f alcon.eps{}\n'.format(sd)
        else:
            smf += ('\trm -f mpxerr.tex mpxerr.log alcon.0 alcon.mp~ '
                'alcon.mpx alcon_mp.log alcon.log{}\n'
                '\trm -f alcon.mp\n\n'
                'cleanall : clean\n'
                '\trm -f alcon.mps alcon.eps alcon.pdf alcon.png\n').format(sd)

        with open(os.path.join(self.input_.dirout, 'Makefile'), 'w') as fmf:
            fmf.write(smf)
    # End of def gnuplot(self):
# End of class Core:
