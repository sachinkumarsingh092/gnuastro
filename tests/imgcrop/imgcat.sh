# Crop from a catalog using x and y coordinates in Image mode.
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing AstrUtils section).





# Preliminaries:
################
# Set the variabels (The executable is in the build tree). Do the
# basic checks to see if the executable is made or if the defaults
# file exists (basicchecks.sh is in the source tree).
prog=astrimgcrop
execname=../src/imgcrop/$prog
source $topsrc/tests/basicchecks.sh





# Actual test script:
#####################
cat=$topsrc/tests/imgcrop/randcat.txt
img=$HOME/Desktop/data/acs_I_030mas_042_sci.fits
$execname $img $cat --imgmode --suffix=_imgcat.fits
