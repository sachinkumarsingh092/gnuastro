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
## along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.

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
imgcrop/section.sh imgcrop/radec.sh imgcrop/imgpolygon.sh		\
imgcrop/imgoutpolygon.sh imgcrop/wcspolygon.sh convertt/fitstotxt.sh	\
convertt/fitstojpeg.sh convertt/blankch.sh convertt/jpegtotxt.sh	\
convertt/fitstojpegcmyk.sh convertt/jpegtofits.sh			\
convertt/fitstopdf.sh convolve/spatial.sh convolve/frequency.sh		\
imgwarp/imgwarp_scale.sh imgwarp/homographic.sh mknoise/addnoise.sh	\
mkprof/ellipticalmasks.sh header/write.sh header/print.sh		\
header/update.sh header/delete.sh imgstat/basicstats.sh			\
subtractsky/subtractsky.sh noisechisel/noisechisel.sh			\
mkcatalog/simple.sh arithmetic/snimage.sh arithmetic/onlynumbers.sh	\
cosmiccal/simpletest.sh

EXTRA_DIST = $(TESTS) mkprof/mkprofcat1.txt			\
mkprof/mkprofcat1_mask.txt mkprof/mkprofcat2.txt		\
mkprof/mkprofcat3.txt mkprof/mkprofcat4.txt imgcrop/cat.txt

CLEANFILES = *.log *.txt *.jpg *.fits *.pdf *.eps .gnuastro/*.conf	\
.gnuastro





## Dependency structure
## ====================
##
## When run with `-jN' option, the tests can be run in parallel, so we need
## to define their dependencies on each other here. The output of all the
## tests that should be used (is recognized by Automake and will not mess
## the system) is the name of the test file appended by a `.log'. It
## contains the output of the program on standard output and error. THIS IS
## IMPORTANT: So even if the program fails, the .log file is created. The
## check if the input for a test exists or not should be checked in the
## test that depends on it, it can't be done here in the Makefile.

mkprof/mosaic1.sh: prepconf.sh.log
mkprof/mosaic2.sh: prepconf.sh.log
mkprof/mosaic3.sh: prepconf.sh.log
mkprof/mosaic4.sh: prepconf.sh.log
imgcrop/imgcat.sh: prepconf.sh.log mkprof/mosaic1.sh.log
imgcrop/wcscat.sh: prepconf.sh.log mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log \
                                   mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
imgcrop/xcyc.sh: prepconf.sh.log mkprof/mosaic1.sh.log
imgcrop/xcycnoblank.sh: prepconf.sh.log mkprof/mosaic1.sh.log
imgcrop/section.sh: prepconf.sh.log mkprof/mosaic1.sh.log
imgcrop/radec.sh: prepconf.sh.log mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log \
                                  mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
imgcrop/imgpolygon.sh: prepconf.sh.log mkprof/mosaic1.sh.log
imgcrop/imgoutpolygon.sh: prepconf.sh.log mkprof/mosaic1.sh.log
imgcrop/wcspolygon.sh: prepconf.sh.log mkprof/mosaic1.sh.log                  \
                       mkprof/mosaic2.sh.log mkprof/mosaic3.sh.log            \
                       mkprof/mosaic4.sh.log
convertt/fitstotxt.sh: prepconf.sh.log mkprof/mosaic1.sh.log
convertt/fitstojpeg.sh: prepconf.sh.log mkprof/mosaic1.sh.log
convertt/blankch.sh: prepconf.sh.log mkprof/mosaic1.sh.log
convertt/jpegtotxt.sh: prepconf.sh.log convertt/fitstojpeg.sh.log
convertt/fitstojpegcmyk.sh: prepconf.sh.log mkprof/mosaic1.sh.log
convertt/jpegtofits.sh: prepconf.sh.log convertt/blankch.sh.log
convertt/fitstopdf.sh: prepconf.sh.log imgcrop/section.sh.log
convolve/spatial.sh: prepconf.sh.log mkprof/mosaic1.sh.log
convolve/frequency.sh: prepconf.sh.log mkprof/mosaic1.sh.log
imgwarp/imgwarp_scale.sh: prepconf.sh.log convolve/spatial.sh.log
imgwarp/homographic.sh: prepconf.sh.log convolve/spatial.sh.log
mknoise/addnoise.sh: prepconf.sh.log imgwarp/imgwarp_scale.sh.log
mkprof/ellipticalmasks.sh: prepconf.sh.log mknoise/addnoise.sh.log
header/write.sh: prepconf.sh.log convolve/spatial.sh.log
header/print.sh: prepconf.sh.log header/write.sh.log
header/update.sh: prepconf.sh.log header/write.sh.log
header/delete.sh: prepconf.sh.log header/write.sh.log
imgstat/basicstats.sh: prepconf.sh.log mknoise/addnoise.sh.log
subtractsky/subtractsky.sh: prepconf.sh.log mknoise/addnoise.sh.log
noisechisel/noisechisel.sh: prepconf.sh.log mknoise/addnoise.sh.log
mkcatalog/simple.sh: prepconf.sh.log noisechisel/noisechisel.sh.log
arithmetic/snimage.sh: prepconf.sh.log noisechisel/noisechisel.sh.log
arithmetic/onlynumbers.sh: prepconf.sh.log
cosmiccal/simpletest.sh: prepconf.sh.log
