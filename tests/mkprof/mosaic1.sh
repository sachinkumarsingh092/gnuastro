# Create a mock image from cat1.txt:
#
#  - It includes one large and bright profile that is on the last
#    pixel of this image. The other tests will also build the same
#    profile with the absolute place fixed.
#
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing AstrUtils section).
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
execname=../src/$prog/astr$prog





# If the executable was not made (the user chose to not install this
# package), skip this test:
if [ ! -f $execname ]; then
    exit 77
fi





# Actual test script:
#####################
cat=$topsrc/tests/$prog/mkprofcat1.txt
$execname $cat --naxis1=100 --naxis2=100
mv 0.fits psf.fits
