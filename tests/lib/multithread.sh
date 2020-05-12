# Run the program to test reading a FITS file to memory and writing it on
# the command-line in multi-threaded mode.
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
execname=./multithread





# SKIP or FAIL?
# =============
#
# If the actual executable wasn't built, then this is a hard error and must
# be FAIL. But if the input doesn't exist, its not this test's fault. So
# just SKIP this test.
if [ ! -f $execname ]; then
    echo "$execname library program not compiled.";
    exit 99;
fi;
if [ ! -f $img      ]; then echo "$img does not exist.";   exit 77; fi;





# Actual test script
# ==================
#
# 'check_with_program' can be something like Valgrind or an empty
# string. Such programs will execute the command if present and help in
# debugging when the developer doesn't have access to the user's system.
$check_with_program $execname
