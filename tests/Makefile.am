## Process this file with automake to produce Makefile.inx
##
## Original author:
##     Mohammad Akhlaghi <akhlaghi@gnu.org>
## Contributing author(s):
## Copyright (C) 2015, Free Software Foundation, Inc.
##
## Gnuastro is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Gnuastro is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with gnuastro. If not, see <http://www.gnu.org/licenses/>.

if COND_HASGHOSTSCRIPT
  MAYBE_HASGHOSTSCRIPT = "yes"
endif
if COND_HASLIBJPEG
  MAYBE_HASLIBJPEG = "yes"
endif

AM_TESTS_ENVIRONMENT =                              \
export topsrc=$(top_srcdir);                        \
export haslibjpeg=$(MAYBE_HASLIBJPEG);              \
export hasghostscript=$(MAYBE_HASGHOSTSCRIPT);

TESTS = prepconf.sh mkprof/mosaic1.sh mkprof/mosaic2.sh			\
mkprof/mosaic3.sh mkprof/mosaic4.sh imgcrop/imgcat.sh			\
imgcrop/wcscat.sh imgcrop/xcyc.sh imgcrop/xcycnoblank.sh		\
imgcrop/section.sh imgcrop/radec.sh convertt/fitstotxt.sh		\
convertt/fitstojpeg.sh convertt/blankch.sh convertt/jpegtotxt.sh	\
convertt/fitstojpegcmyk.sh convertt/jpegtofits.sh			\
convertt/fitstopdf.sh convolve/spatial.sh convolve/frequency.sh		\
imgwarp/imgwarp_scale.sh imgwarp/homographic.sh mknoise/addnoise.sh

EXTRA_DIST = $(TESTS) mkprof/mkprofcat1.txt mkprof/mkprofcat2.txt	\
mkprof/mkprofcat3.txt mkprof/mkprofcat4.txt imgcrop/cat.txt

CLEANFILES = *.log *.txt *.jpg *.fits *.pdf *.eps .gnuastro/*.conf	\
.gnuastro
