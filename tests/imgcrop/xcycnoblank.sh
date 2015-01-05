# Crop a box based on --xc and --yc, with no blank pixels.
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing AstrUtils section).





# Preliminaries:
################
# Set the variabels (The executable is in the build tree). Do the
# basic checks to see if the executable is made or if the defaults
# file exists (basicchecks.sh is in the source tree).
prog=imgcrop
execname=../src/$prog/astr$prog
source $topsrc/tests/basicchecks.sh





# Actual test script:
#####################
img=$HOME/Desktop/data/acs_I_030mas_042_sci.fits
$execname $img --xc=20 --yc=20 --noblank --output=imgcrop_xcycnb.fits
