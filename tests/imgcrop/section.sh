# Crop a box based on a specified section.
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
img=mkprofcat1.fits
$execname $img --section=-10:*+10,:250 --output=imgcrop_section.fits
