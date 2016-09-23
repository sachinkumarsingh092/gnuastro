# Convert a FITS image into a JPEG file.
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
# Set the variabels (The executable is in the build tree). Do the
# basic checks to see if the executable is made or if the defaults
# file exists (basicchecks.sh is in the source tree).
img=psf.fits
prog=convertt
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
#
#   - libjpeg was not present on the system.
if [ ! -f $execname ] || [ ! -f $img ] || [ "x$haslibjpeg" != "xyes" ]; then
    exit 77;
fi






# Actual test script
# ==================
$execname blank blank blank $img --output=f2jcmyk.jpg
