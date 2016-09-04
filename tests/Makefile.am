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





# Translate conditions that came from `./configure' into variables that can
# be used here.
if COND_HASGHOSTSCRIPT
  MAYBE_HASGHOSTSCRIPT = "yes"
endif
if COND_HASLIBJPEG
  MAYBE_HASLIBJPEG = "yes"
endif





# Environment variables for the test scripts.
AM_TESTS_ENVIRONMENT =                              \
export topsrc=$(top_srcdir);                        \
export haslibjpeg=$(MAYBE_HASLIBJPEG);              \
export hasghostscript=$(MAYBE_HASGHOSTSCRIPT);





# Compilations that are to be done with `make check'.
LDADD = -lgnuastro
check_PROGRAMS = versionc versioncpp
versionc_SOURCES = lib/versionc.c
versioncpp_SOURCES = lib/versioncpp.cpp




# The actual test scripts that are run:
TESTS = $(check_PROGRAMS) prepconf.sh mkprof/mosaic1.sh mkprof/mosaic2.sh     \
  mkprof/mosaic3.sh mkprof/mosaic4.sh mkprof/radeccat.sh imgcrop/imgcat.sh    \
  imgcrop/wcscat.sh imgcrop/xcyc.sh imgcrop/xcycnoblank.sh                    \
  imgcrop/section.sh imgcrop/radec.sh imgcrop/imgpolygon.sh                   \
  imgcrop/imgoutpolygon.sh imgcrop/wcspolygon.sh convertt/fitstotxt.sh        \
  convertt/fitstojpeg.sh convertt/blankch.sh convertt/jpegtotxt.sh            \
  convertt/fitstojpegcmyk.sh convertt/jpegtofits.sh convertt/fitstopdf.sh     \
  convolve/spatial.sh convolve/frequency.sh imgwarp/imgwarp_scale.sh          \
  imgwarp/homographic.sh mknoise/addnoise.sh mkprof/ellipticalmasks.sh        \
  mkprof/inputascanvas.sh header/write.sh header/print.sh header/update.sh    \
  header/delete.sh imgstat/basicstats.sh subtractsky/subtractsky.sh           \
  noisechisel/noisechisel.sh mkcatalog/simple.sh mkcatalog/aperturephot.sh    \
  arithmetic/snimage.sh arithmetic/onlynumbers.sh arithmetic/where.sh         \
  cosmiccal/simpletest.sh table/asciitobinary.sh table/binarytoascii.sh





# Files to distribute along with the tests.
EXTRA_DIST = $(TESTS) during-dev.sh mkprof/mkprofcat1.txt                   \
  mkprof/ellipticalmasks.txt mkprof/inputascanvas.txt                       \
  mkprof/mkprofcat2.txt mkprof/mkprofcat3.txt mkprof/mkprofcat4.txt         \
  mkprof/radeccat.txt imgcrop/cat.txt table/asciitobinary.txt





# Files that must be cleaned with `make clean'.
CLEANFILES = *.log *.txt *.jpg *.fits *.pdf *.eps





# CLEANFILES is only for files, not directories. Therefore we are using
# Automake's extending rules to clean the temporary `.gnuastro' directory
# that was built by the `prepconf.sh' scripot. See "Extending Automake
# rules", and the "What Gets Cleaned" sections of the Automake manual.
clean-local:; rm -rf .gnuastro





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
mkprof/radeccat.sh: prepconf.sh.log
table/asciitobinary.sh: prepconf.sh.log
table/binarytoascii.sh: table/asciitobinary.sh.log
imgcrop/imgcat.sh: mkprof/mosaic1.sh.log
imgcrop/wcscat.sh: mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log     \
                   mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
imgcrop/xcyc.sh: mkprof/mosaic1.sh.log
imgcrop/xcycnoblank.sh: mkprof/mosaic1.sh.log
imgcrop/section.sh: mkprof/mosaic1.sh.log
imgcrop/radec.sh: mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log      \
                  mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
imgcrop/imgpolygon.sh: mkprof/mosaic1.sh.log
imgcrop/imgoutpolygon.sh: mkprof/mosaic1.sh.log
imgcrop/wcspolygon.sh: mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log \
                       mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
convertt/fitstotxt.sh: mkprof/mosaic1.sh.log
convertt/fitstojpeg.sh: mkprof/mosaic1.sh.log
convertt/blankch.sh: mkprof/mosaic1.sh.log
convertt/jpegtotxt.sh: convertt/fitstojpeg.sh.log
convertt/fitstojpegcmyk.sh: mkprof/mosaic1.sh.log
convertt/jpegtofits.sh: convertt/blankch.sh.log
convertt/fitstopdf.sh: imgcrop/section.sh.log
convolve/spatial.sh: mkprof/mosaic1.sh.log
convolve/frequency.sh: mkprof/mosaic1.sh.log
imgwarp/imgwarp_scale.sh: convolve/spatial.sh.log
imgwarp/homographic.sh: convolve/spatial.sh.log
mknoise/addnoise.sh: imgwarp/imgwarp_scale.sh.log
mkprof/ellipticalmasks.sh: mknoise/addnoise.sh.log
mkprof/inputascanvas.sh: mknoise/addnoise.sh.log
header/write.sh: convolve/spatial.sh.log
header/print.sh: header/write.sh.log
header/update.sh: header/write.sh.log
header/delete.sh: header/write.sh.log
imgstat/basicstats.sh: mknoise/addnoise.sh.log
subtractsky/subtractsky.sh: mknoise/addnoise.sh.log
noisechisel/noisechisel.sh: mknoise/addnoise.sh.log
mkcatalog/simple.sh: noisechisel/noisechisel.sh.log
mkcatalog/aperturephot.sh: noisechisel/noisechisel.sh.log          \
                           mkprof/inputascanvas.sh
arithmetic/snimage.sh: noisechisel/noisechisel.sh.log
arithmetic/where.sh: noisechisel/noisechisel.sh.log
arithmetic/onlynumbers.sh: prepconf.sh.log
cosmiccal/simpletest.sh: prepconf.sh.log
