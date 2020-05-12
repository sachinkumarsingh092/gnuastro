# Make labeled regions on an image with blank pixels.
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing gnuastro section).
#
# Original author:
#     Mohammad Akhlaghi <mohammad@akhlaghi.org>
# Contributing author(s):
# Copyright (C) 2015-2020, Free Software Foundation, Inc.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.





# Preliminaries
# =============
#
# Set the variables (The executable is in the build tree). Do the
# basic checks to see if the executable is made or if the defaults
# file exists (basicchecks.sh is in the source tree).
img=psf.fits
prog=buildprog
execname=../bin/$prog/ast$prog
source=$topsrc/tests/$prog/simpleio.c





# Skip?
# =====
#
# If the dependencies of the test don't exist, then skip it. There are two
# types of dependencies:
#
#   - The executable was not made (for example due to a configure option).
#
#   - The input data was not made (for example the test that created the
#     data file failed).
if [ ! -f $execname ]; then echo "$execname not created.";  exit 77; fi
if [ ! -f $img      ]; then echo "$img does not exist.";    exit 77; fi
if [ ! -f $source   ]; then echo "$source does not exist."; exit 77; fi





# Actual test script
# ==================
#
# We want to use the 'libgnuastro.la' corresponding to this install, not
# the one (that is possibly) installed (hence the use of '--la').
#
# Except for 'gnuastro/config.h', all headers are installed in
# '$topsrc/lib' and 'gnuastro/config.h' is in "../lib/"
#
# 'check_with_program' can be something like Valgrind or an empty
# string. Such programs will execute the command if present and help in
# debugging when the developer doesn't have access to the user's system.
echo "Test Environment"
echo "----------------"
echo "CPPFLAGS: $CPPFLAGS"
echo "LDFLAGS: $LDFLAGS"
echo "----------------"
$check_with_program $execname $source $img 1 --la=../lib/libgnuastro.la \
                              -I$topsrc/lib -I../lib/
