# List all FITS files in build directory by date.
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
#
# We will be adding noise to two images: the warped (smaller) and unwarped
# (larger) mock images. The warped one will be used by programs that don't
# care about the size of the image, but the larger one will be used by
# those that do: for example SubtractSky and NoiseChisel will be better
# tested on a larger image.
prog=sort-by-night
dep1=fits
dep1name=../bin/$dep1/ast$dep1
execname=../bin/script/astscript-$prog





# Skip?
# =====
#
# If the dependencies of the test don't exist, then skip it. There are two
# types of dependencies:
#
#   - The executable script was not made.
#   - The programs it use weren't made.
if [ ! -f $execname ]; then echo "$execname doesn't exist."; exit 77; fi
if [ ! -f $dep1name ]; then echo "$dep1name doesn't exist."; exit 77; fi





# Put a link of Gnuastro program(s) used into current directory. Note that
# other script tests may have already brought it.
ln -sf $dep1name ast$dep1



# Actual test script
# ==================
#
# 'check_with_program' can be something like Valgrind or an empty
# string. Such programs will execute the command if present and help in
# debugging when the developer doesn't have access to the user's system.
#
# Since we want the script to recognize the programs that it will use from
# this same build of Gnuastro, we'll add the current directory to PATH.
export PATH="./:$PATH"
$check_with_program $execname *.fits
