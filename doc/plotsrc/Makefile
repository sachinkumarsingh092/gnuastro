# Prepare all the figures for the GNU Astronomy Utilities manual.
#
# Original author:
#     Mohammad Akhlaghi <mohammad@akhlaghi.org>
# Contributing author(s):
# Copyright (C) 2015-2020, Free Software Foundation, Inc.
#
# GNU Astronomy Utilities is free software: you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# GNU Astronomy Utilities is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Astronomy Utilities. If not, see
# <http://www.gnu.org/licenses/>.

# Definitions of which rules to apply:
all: all.pdf
.PHONY: clean





# Make all the plots and convert all the figures
all.pdf: all.tex ./tex/*.tex ./conversions.sh

#       Make the tikz directory if it doesn't exist:
	mkdir -p tikz

#       Run the LaTeX script:
	latex -shell-escape all.tex
	dvipdf all.dvi
	rm -f *.bcf *.xml *.aux *.auxlock *.log *.out *.dvi

#       Put all the EPS files in their appropriate place:
	cp tikz/all-figure0.eps ../gnuastro-figures/iandtime.eps
	cp tikz/all-figure1.eps ../gnuastro-figures/samplingfreq.eps
	cp tikz/all-figure2.eps ../gnuastro-figures/flatplane.eps
	cp tikz/all-figure3.eps ../gnuastro-figures/sphereandplane.eps

#       Do all the conversions:
	./conversions.sh ../gnuastro-figures/





# Clean all generated plots and conversions:
#   1. Copy all the actual PNG files into this directory.
#        NOTE I:  All non-plot figures should be saved in PNG.
#        NOTE II: All non-plot PNG figures should be listed here.
#   2. Remove everything in the gnuastro-figures directory
#   3. Put the actual images back in that directory.
clean:
	rm -rf all.pdf tikz
	cp ../gnuastro-figures/epicycles.png ./
	rm ../gnuastro-figures/*
	mv epicycles.png ../gnuastro-figures/
