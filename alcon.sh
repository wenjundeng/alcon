#!/bin/bash
# This script serves as an input file feeding alcon.py

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

# To get the full description of each input parameter, run: python alcon.py -h
cmd=(python alcon.py
	--eqtype alcon.dat  # equilibrium type, available choices: alcon.dat, sa
	--finitebeta 3  # finite beta type, choose from 0 to 3
	--radrange 0.01 1.0  # solving radial coordinate range
	--nrad 2000  # number of radial grid points
	--ntor 4  # toroidal mode number n
	--mpolrange -20 50  # poloidal mode number m range
	--noffdiag 20  # number of off-diagonal stripes
	--ncon 50  # number of lowest continua to solve
	--imreratiocutoff 0.1  # imaginary to real ratio cutoff value
	--sigma -0.5  # sigma for shift-invert mode in matrix eigen solver
	--omegascale 1.0  # scale factor for omega before output
	--omegacutoff 1.0  # cut off value for scaled omega
	--dirout ./output  # directory for output
	#--erasedirout  # erase directory for output if preexists
	-n 1  # number of processes for parallel solving
	)
"${cmd[@]}"
