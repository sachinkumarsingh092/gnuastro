# Create a 3D kernel with MakeProfiles.
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing gnuastro section).
#
# Original author:
#     Mohammad Akhlaghi <akhlaghi@gnu.org>
# Contributing author(s):
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
prog=mkprof
execname=../bin/$prog/ast$prog





# Skip?
# =====
#
# If the dependencies of the test don't exist, then skip it. There are two
# types of dependencies:
#
#   - The executable was not made (for example due to a configure option).
#
#   - Catalog doesn't exist (problem in tarball release).
if [ ! -f $execname ]; then echo "$execname not created."; exit 77; fi





# Actual test script
# ==================
$execname $cat --kernel=gaussian-3d,2,3,0.5 --oversample=1    \
          --output=3d-kernel.fits
