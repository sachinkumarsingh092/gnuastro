# Create a mock image from cat2.txt:
#
#  - It doesn't have any random profiles, only the large profile from
#    cat1.txt. The central position is set to be on the same real
#    place as in cat1.txt
#
#  - Also here, test the individual option.
#
#
# See the Tests subsection of the manual for a complete explanation
# (in the Installing AstrUtils section).


# Preliminaries:
################
# Set the variabels (The executable is in the build tree). Do the
# basic checks to see if the executable is made or if the defaults
# file exists (basicchecks.sh is in the source tree).
prog=mkprof
execname=../src/$prog/astr$prog
source $topsrc/tests/basicchecks.sh


# Actual test script:
#####################
cat=$topsrc/tests/$prog/mkprofcat2.txt
$execname $cat --naxis1=100 --naxis2=100 --crpix1=-99 --individual
