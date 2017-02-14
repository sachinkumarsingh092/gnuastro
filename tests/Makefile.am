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





# `config.h' is built by `../lib/Makefile.am' and is only meant for
# outside users (internally `config.h' has everything. Thus, we need to add
# it to the include directories of the test programs in this directory.
AM_CPPFLAGS = -I../lib





# Environment variables for the test scripts.
AM_TESTS_ENVIRONMENT =                              \
export topsrc=$(top_srcdir);                        \
export haslibjpeg=$(MAYBE_HASLIBJPEG);              \
export hasghostscript=$(MAYBE_HASGHOSTSCRIPT);





# Compilations that are to be done with `make check'.
LDADD = -lgnuastro
check_PROGRAMS = versioncpp arraymanip multithread
arraymanip_SOURCES = lib/arraymanip.c
versioncpp_SOURCES = lib/versioncpp.cpp
multithread_SOURCES = lib/multithread.c




# The actual test scripts that are run:
TESTS = $(check_PROGRAMS) prepconf.sh mkprof/mosaic1.sh mkprof/mosaic2.sh  \
  mkprof/mosaic3.sh mkprof/mosaic4.sh mkprof/radeccat.sh crop/imgcat.sh    \
  crop/wcscat.sh crop/xcyc.sh crop/xcycnoblank.sh crop/section.sh          \
  crop/radec.sh crop/imgpolygon.sh crop/imgoutpolygon.sh                   \
  crop/wcspolygon.sh convertt/fitstotxt.sh convertt/fitstojpeg.sh          \
  convertt/blankch.sh convertt/jpegtotxt.sh convertt/fitstojpegcmyk.sh     \
  convertt/jpegtofits.sh convertt/fitstopdf.sh convolve/spatial.sh         \
  convolve/frequency.sh warp/warp_scale.sh warp/homographic.sh             \
  mknoise/addnoise.sh mkprof/ellipticalmasks.sh mkprof/clearcanvas.sh      \
  header/write.sh header/print.sh header/update.sh header/delete.sh        \
  statistics/basicstats.sh subtractsky/subtractsky.sh                      \
  noisechisel/noisechisel.sh mkcatalog/simple.sh mkcatalog/aperturephot.sh \
  arithmetic/snimage.sh arithmetic/onlynumbers.sh arithmetic/where.sh      \
  arithmetic/or.sh cosmiccal/simpletest.sh table/txt-to-fits-binary.sh     \
  table/fits-binary-to-txt.sh table/txt-to-fits-ascii.sh                   \
  table/fits-ascii-to-txt.sh





# Files to distribute along with the tests.
EXTRA_DIST = $(TESTS) during-dev.sh mkprof/mkprofcat1.txt                 \
  mkprof/ellipticalmasks.txt mkprof/clearcanvas.txt mkprof/mkprofcat2.txt \
  mkprof/mkprofcat3.txt mkprof/mkprofcat4.txt mkprof/radeccat.txt         \
  crop/cat.txt table/table.txt





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
table/txt-to-fits-binary.sh: prepconf.sh.log
table/fits-binary-to-txt.sh: table/txt-to-fits-binary.sh.log
table/txt-to-fits-ascii.sh: prepconf.sh.log
table/fits-ascii-to-txt.sh: table/txt-to-fits-ascii.sh.log
crop/imgcat.sh: mkprof/mosaic1.sh.log
crop/wcscat.sh: mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log     \
                   mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
crop/xcyc.sh: mkprof/mosaic1.sh.log
crop/xcycnoblank.sh: mkprof/mosaic1.sh.log
crop/section.sh: mkprof/mosaic1.sh.log
crop/radec.sh: mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log      \
                  mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
crop/imgpolygon.sh: mkprof/mosaic1.sh.log
crop/imgoutpolygon.sh: mkprof/mosaic1.sh.log
crop/wcspolygon.sh: mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log \
                       mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
convertt/fitstotxt.sh: mkprof/mosaic1.sh.log
convertt/fitstojpeg.sh: mkprof/mosaic1.sh.log
convertt/blankch.sh: mkprof/mosaic1.sh.log
convertt/jpegtotxt.sh: convertt/fitstojpeg.sh.log
convertt/fitstojpegcmyk.sh: mkprof/mosaic1.sh.log
convertt/jpegtofits.sh: convertt/blankch.sh.log
convertt/fitstopdf.sh: crop/section.sh.log
convolve/spatial.sh: mkprof/mosaic1.sh.log
convolve/frequency.sh: mkprof/mosaic1.sh.log
warp/warp_scale.sh: convolve/spatial.sh.log
warp/homographic.sh: convolve/spatial.sh.log
mknoise/addnoise.sh: warp/warp_scale.sh.log
mkprof/ellipticalmasks.sh: mknoise/addnoise.sh.log
mkprof/clearcanvas.sh: mknoise/addnoise.sh.log
header/write.sh: mkprof/mosaic1.sh.log
header/print.sh: header/write.sh.log
header/update.sh: header/write.sh.log
header/delete.sh: header/write.sh.log
statistics/basicstats.sh: mknoise/addnoise.sh.log
subtractsky/subtractsky.sh: mknoise/addnoise.sh.log
noisechisel/noisechisel.sh: mknoise/addnoise.sh.log
mkcatalog/simple.sh: noisechisel/noisechisel.sh.log
mkcatalog/aperturephot.sh: noisechisel/noisechisel.sh.log          \
                           mkprof/clearcanvas.sh
arithmetic/snimage.sh: noisechisel/noisechisel.sh.log
arithmetic/where.sh: noisechisel/noisechisel.sh.log
arithmetic/or.sh: noisechisel/noisechisel.sh.log
arithmetic/onlynumbers.sh: prepconf.sh.log
cosmiccal/simpletest.sh: prepconf.sh.log
