# Add noise to an input cube.
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
#
# We will be adding noise to two images: the warped (smaller) and unwarped
# (larger) mock images. The warped one will be used by programs that don't
# care about the size of the image, but the larger one will be used by
# those that do: for example SubtractSky and NoiseChisel will be better
# tested on a larger image.
prog=mknoise
img=3d-cat.fits
execname=../bin/$prog/ast$prog





# Skip?
# =====
#
# If the dependencies of the test don't exist, then skip it. There are two
# types of dependencies:
#
#   - The executable was not made (for example due to a configure option),
#
#   - The input data was not made (for example the test that created the
#     data file failed).
if [ ! -f $execname ]; then echo "$execname doesn't exist."; exit 77; fi
if [ ! -f $img      ]; then echo "$img does not exist.";     exit 77; fi





# Actual test script
# ==================
export GSL_RNG_SEED=1
export GSL_RNG_TYPE=ranlxs2
$execname --envseed $img
