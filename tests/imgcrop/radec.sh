# Crop a box based on the --xc and --yc options.
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing gnuastro section).
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
prog=imgcrop
execname=../src/$prog/ast$prog





# If the executable was not made (the user chose to not install this
# package), skip this test:
if [ ! -f $execname ]; then
    exit 77
fi





# Actual test script:
#####################
img=mkprofcat*.fits
$execname $img --ra=0.99917157 --dec=1.0008283 --wwidth=0.3 --output=imgcrop_radec.fits
