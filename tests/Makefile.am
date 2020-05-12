## Process this file with automake to produce Makefile.inx
##
## Original author:
##     Mohammad Akhlaghi <mohammad@akhlaghi.org>
## Contributing author(s):
## Copyright (C) 2015-2020, Free Software Foundation, Inc.
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





## Translate conditions that came from './configure' into variables that can
## be used here. This includes the programs and teh dependencies of the
## checks.
##
## Dependency structure
## ====================
##
## The tests depend on each other. One way to fix this is to call them in
## the order that their dependence is satisfied. But that is prone to
## errors and a lot of problem. The user also not want to build all the
## programs. Also, Make allows us to work with multiple threads (with
## '-jN') and in that case, it is vital to define the dependencies.
##
## The output of all the tests that should be used (is recognized by
## Automake and will not mess the system) is the name of the test file
## appended by a '.log'. It contains the output of the program on standard
## output and error. THIS IS IMPORTANT: So even if the program fails, the
## .log file is created. The check if the input for a test exists or not
## should be checked in the test that depends on it, it can't be done here
## in the Makefile.
if COND_CHECK_WITH_VALGRIND
  MAYBE_CHECK_WITH_PROGRAM = "valgrind"
endif
if COND_HASGHOSTSCRIPT
  MAYBE_HASGHOSTSCRIPT = "yes"
endif
if COND_HASLIBJPEG
  MAYBE_HASLIBJPEG = "yes"
endif
if COND_HASCXX
  MAYBE_CXX_PROGS    = versioncxx
  MAYBE_CXX_TESTS    = lib/versioncxx.sh
  versioncxx_SOURCES = lib/versioncxx.cpp
endif
if COND_ARITHMETIC
  MAYBE_ARITHMETIC_TESTS = arithmetic/snimage.sh arithmetic/onlynumbers.sh \
  arithmetic/where.sh arithmetic/or.sh arithmetic/connected-components.sh

  arithmetic/onlynumbers.sh: prepconf.sh.log
  arithmetic/connected-components.sh: noisechisel/noisechisel.sh.log
  arithmetic/snimage.sh: noisechisel/noisechisel.sh.log
  arithmetic/where.sh: noisechisel/noisechisel.sh.log
  arithmetic/or.sh: segment/segment.sh.log
endif
if COND_BUILDPROG
  MAYBE_BUILDPROG_TESTS = buildprog/simpleio.sh

  buildprog/simpleio.sh: mkprof/mosaic1.sh.log
endif
if COND_CONVERTT
  MAYBE_CONVERTT_TESTS = convertt/fitstotxt.sh convertt/fitstojpeg.sh	\
  convertt/blankch.sh convertt/jpegtotxt.sh convertt/fitstojpegcmyk.sh	\
  convertt/jpegtofits.sh convertt/fitstopdf.sh

  convertt/fitstotxt.sh: mkprof/mosaic1.sh.log
  convertt/fitstojpeg.sh: mkprof/mosaic1.sh.log
  convertt/blankch.sh: mkprof/mosaic1.sh.log
  convertt/jpegtotxt.sh: convertt/fitstojpeg.sh.log
  convertt/fitstojpegcmyk.sh: mkprof/mosaic1.sh.log
  convertt/jpegtofits.sh: convertt/blankch.sh.log
  convertt/fitstopdf.sh: crop/section.sh.log
endif
if COND_CONVOLVE
  MAYBE_CONVOLVE_TESTS = convolve/spatial.sh convolve/frequency.sh

  convolve/spatial.sh: mkprof/mosaic1.sh.log
  convolve/frequency.sh: mkprof/mosaic1.sh.log
endif
if COND_COSMICCAL
  MAYBE_COSMICCAL_TESTS = cosmiccal/simpletest.sh

  cosmiccal/simpletest.sh: prepconf.sh.log
endif
if COND_CROP
  MAYBE_CROP_TESTS = crop/imgcat.sh crop/wcscat.sh crop/imgcenter.sh    \
  crop/imgcenternoblank.sh crop/section.sh crop/wcscenter.sh            \
  crop/imgpolygon.sh crop/imgpolygonout.sh crop/wcspolygon.sh

  crop/imgcat.sh: mkprof/mosaic1.sh.log
  crop/wcscat.sh: mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log     \
                  mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
  crop/imgcenter.sh: mkprof/mosaic1.sh.log
  crop/imgcenternoblank.sh: mkprof/mosaic1.sh.log
  crop/section.sh: mkprof/mosaic1.sh.log
  crop/wcscenter.sh: mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log      \
                     mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
  crop/imgpolygon.sh: mkprof/mosaic1.sh.log
  crop/imgpolygonout.sh: mkprof/mosaic1.sh.log
  crop/wcspolygon.sh: mkprof/mosaic1.sh.log mkprof/mosaic2.sh.log \
                      mkprof/mosaic3.sh.log mkprof/mosaic4.sh.log
endif
if COND_FITS
  MAYBE_FITS_TESTS = fits/write.sh fits/print.sh fits/update.sh	\
  fits/delete.sh fits/copyhdu.sh

  fits/write.sh: mkprof/mosaic1.sh.log
  fits/print.sh: fits/write.sh.log
  fits/update.sh: fits/write.sh.log
  fits/delete.sh: fits/write.sh.log
  fits/copyhdu.sh: fits/write.sh.log mkprof/mosaic2.sh.log
endif
if COND_MATCH
  MAYBE_MATCH_TESTS = match/positions.sh match/merged-cols.sh

  match/positions.sh: prepconf.sh.log
  match/merged-cols.sh: prepconf.sh.log
endif
if COND_MKCATALOG
  MAYBE_MKCATALOG_TESTS = mkcatalog/detections.sh   \
  mkcatalog/objects-clumps.sh mkcatalog/aperturephot.sh

  mkcatalog/objects-clumps.sh: segment/segment.sh.log
  mkcatalog/detections.sh: arithmetic/connected-components.sh.log
  mkcatalog/aperturephot.sh: noisechisel/noisechisel.sh.log          \
                             mkprof/clearcanvas.sh.log
endif
if COND_MKNOISE
  MAYBE_MKNOISE_TESTS = mknoise/addnoise.sh mknoise/addnoise-3d.sh

  mknoise/addnoise.sh: warp/warp_scale.sh.log
  mknoise/addnoise-3d.sh: mkprof/3d-cat.sh.log
endif
if COND_MKPROF
  MAYBE_MKPROF_TESTS = mkprof/mosaic1.sh mkprof/mosaic2.sh         \
  mkprof/mosaic3.sh mkprof/mosaic4.sh mkprof/radeccat.sh           \
  mkprof/ellipticalmasks.sh mkprof/clearcanvas.sh mkprof/3d-cat.sh \
  mkprof/3d-kernel.sh

  mkprof/3d-cat.sh: prepconf.sh.log
  mkprof/mosaic1.sh: prepconf.sh.log
  mkprof/mosaic2.sh: prepconf.sh.log
  mkprof/mosaic3.sh: prepconf.sh.log
  mkprof/mosaic4.sh: prepconf.sh.log
  mkprof/radeccat.sh: prepconf.sh.log
  mkprof/3d-kernel.sh: prepconf.sh.log
  mkprof/ellipticalmasks.sh: mknoise/addnoise.sh.log
  mkprof/clearcanvas.sh: mknoise/addnoise.sh.log
endif
if COND_NOISECHISEL
  MAYBE_NOISECHISEL_TESTS = noisechisel/noisechisel.sh

  noisechisel/noisechisel.sh: mknoise/addnoise.sh.log
endif
if COND_SEGMENT
  MAYBE_SEGMENT_TESTS = segment/segment.sh

  segment/segment.sh: noisechisel/noisechisel.sh.log
endif
if COND_STATISTICS
  MAYBE_STATISTICS_TESTS = statistics/basicstats.sh statistics/estimate_sky.sh \
                           statistics/from-stdin.sh

  statistics/basicstats.sh: mknoise/addnoise.sh.log
  statistics/estimate_sky.sh: mknoise/addnoise.sh.log
  statistics/from-stdin.sh: prepconf.sh.log
endif
if COND_TABLE
  MAYBE_TABLE_TESTS = table/txt-to-fits-binary.sh		\
  table/fits-binary-to-txt.sh table/txt-to-fits-ascii.sh	\
  table/fits-ascii-to-txt.sh

  table/txt-to-fits-binary.sh: prepconf.sh.log
  table/fits-binary-to-txt.sh: table/txt-to-fits-binary.sh.log
  table/txt-to-fits-ascii.sh: prepconf.sh.log
  table/fits-ascii-to-txt.sh: table/txt-to-fits-ascii.sh.log
endif
if COND_WARP
  MAYBE_WARP_TESTS = warp/warp_scale.sh warp/homographic.sh

  warp/warp_scale.sh: convolve/spatial.sh.log
  warp/homographic.sh: convolve/spatial.sh.log
endif

# Script tests.
SCRIPT_TESTS = script/list-by-night.sh

# We want to have several FITS files as input for this script.
script/list-by-night.sh: mkcatalog/aperturephot.sh.log




# Environment variables for the test scripts.
AM_TESTS_ENVIRONMENT =                                   \
export mkdir_p="$(MKDIR_P)";                             \
export topsrc=$(top_srcdir);                             \
export topbuild=$(top_builddir);                         \
export haslibjpeg=$(MAYBE_HASLIBJPEG);                   \
export hasghostscript=$(MAYBE_HASGHOSTSCRIPT);           \
export check_with_program=$(MAYBE_CHECK_WITH_PROGRAM);





# Library checks
# ==============
#
# The Gnuastro library is checked by compiling programs and linking them
# with the library. As described in the last paragraph of the "Scripts
# based test suites" section of the Automake manual, all targets specified
# by 'check_PROGRAMS' are compiled prior to actually running the targets of
# 'TESTS'. So they do not need to be specified as any dependency, they will
# be present when the '.sh' based tests are run.

# The 'gnuastro/config.h' (needed by Gnuastro's library) is built by
# '../lib/Makefile.am' and is only meant for outside users (to be tested
# here). Thus (unlike the programs, which use 'config.h') we need to add
# the top build directory to the include search directories ('-I').
LDADD = -lgnuastro $(CONFIG_LDADD)
AM_LDFLAGS = -L\$(top_builddir)/lib
AM_CPPFLAGS = -I\$(top_srcdir)/lib -I\$(top_builddir)/lib

# Rest of library check settings.
check_PROGRAMS = multithread $(MAYBE_CXX_PROGS)
multithread_SOURCES = lib/multithread.c
lib/multithread.sh: mkprof/mosaic1.sh.log





# Final Tests
# ===========
TESTS = prepconf.sh lib/multithread.sh $(MAYBE_CXX_TESTS)                  \
  $(MAYBE_ARITHMETIC_TESTS) $(MAYBE_BUILDPROG_TESTS)                       \
  $(MAYBE_CONVERTT_TESTS) $(MAYBE_CONVOLVE_TESTS) $(MAYBE_COSMICCAL_TESTS) \
  $(MAYBE_CROP_TESTS) $(MAYBE_FITS_TESTS) $(MAYBE_MATCH_TESTS)             \
  $(MAYBE_MKCATALOG_TESTS) $(MAYBE_MKNOISE_TESTS) $(MAYBE_MKPROF_TESTS)    \
  $(MAYBE_NOISECHISEL_TESTS) $(MAYBE_SEGMENT_TESTS)                        \
  $(MAYBE_STATISTICS_TESTS) $(MAYBE_SUBTRACTSKY_TESTS)                     \
  $(MAYBE_TABLE_TESTS) $(MAYBE_WARP_TESTS) $(SCRIPT_TESTS)





# Files to distribute along with the tests.
EXTRA_DIST = $(TESTS) during-dev.sh buildprog/simpleio.c crop/cat.txt     \
  mkprof/3d-cat.txt match/positions-1.txt match/positions-2.txt           \
  mkprof/mkprofcat1.txt mkprof/ellipticalmasks.txt mkprof/clearcanvas.txt \
  mkprof/mkprofcat2.txt mkprof/mkprofcat3.txt mkprof/mkprofcat4.txt       \
  mkprof/radeccat.txt statistics/stdin-input.txt table/table.txt





# Files that must be cleaned with 'make clean'.
CLEANFILES = *.log *.txt *.jpg *.fits *.pdf *.eps simpleio





# CLEANFILES is only for files, not directories. Therefore we are using
# Automake's extending rules to clean the temporary '.gnuastro' directory
# that was built by the 'prepconf.sh' scripot. See "Extending Automake
# rules", and the "What Gets Cleaned" sections of the Automake manual.
clean-local:; rm -rf .gnuastro
