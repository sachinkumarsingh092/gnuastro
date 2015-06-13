# Detect objects and clumps in an image using NoiseChisel.
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





# Preliminaries:
################
# Set the variabels (The executable is in the build tree). Do the
# basic checks to see if the executable is made or if the defaults
# file exists (basicchecks.sh is in the source tree).
prog=noisechisel
execname=../src/$prog/ast$prog





# If the executable was not made (the user chose to not install this
# package), skip this test:
if [ ! -f $execname ]; then
    exit 77
fi





# Actual test script:
#####################
img=convolve_spatial_noised.fits

# All of these conditions are put here because the image is relatively
# small, so we can't get good statistics. But that doesn't matter too
# much here, this is just a test to see if everything is working
# properly.
$execname $img --lmeshsize=100 --minbfrac=0.5 --detsnminarea=10   \
          --minnumfalse=50 --numnearest=5 --segsnminarea=15       \
          --checkdetection
