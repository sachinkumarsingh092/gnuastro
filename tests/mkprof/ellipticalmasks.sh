# Make two masks (one flat and one circumference) cat1.txt
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
prog=mkprof
execname=../src/$prog/ast$prog





# If the executable was not made (the user chose to not install this
# package), skip this test:
if [ ! -f $execname ]; then
    exit 77
fi



# Actual test script:
#####################
img=convolve_spatial_warped_noised.fits
cat=$topsrc/tests/$prog/mkprofcat1_mask.txt
$execname $cat $img --setconsttomin --replace --oversample=1 \
          --output="masks.fits"
