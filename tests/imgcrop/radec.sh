# Crop a box based on the --xc and --yc options.
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
img=$HOME/Desktop/data/*sci.fits
$execname $img --ra=150.06261 --dec=1.7911417 --output=imgcrop_radec.fits
