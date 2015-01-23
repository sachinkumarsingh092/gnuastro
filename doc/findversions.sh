#! /bin/bash
# Read the version numbers of all the sub packages.
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
#
# All the file names containing the version numbers are input as
# command line arguments. The directory part of all the names should
# be given ast the first argument:
#
# findversions.sh DIRname FILEname(s)...
#
# We want the final file ($output) to remain intact until the versions
# of all the subpackages have been found and written. So we will first
# create a blank temporary file to save the one-by-one outputs, then
# in the end (after it has successfully finished), copy that temporary
# file into the output file.
output="subpackageversions.texi"
tmp="tmp.texi"



# Remove the temporary file if it exists. It is very important that it
# be newly made in the next loop, because all the outputs will be
# appended to the file. We don't want anything remaining.
if [ -f $tmp ]; then
    rm $tmp
fi



echo
echo Current GNU Astronomy Utilities versions are:

#Read the subpackage version for each program:
for filename in $@
do
    # Save the output of grep (a line) as an array:
    a=($(grep "define SPACK_VERSION" $filename))
    b=${a[2]}   # Save the second token (word):
    c=${b#\"}   # Remove the initial "
    progversion=${c%\"}   # Remove the final "

    #Find the program name:
    a=($(grep "define SPACK " $filename))
    b=${a[2]}   # Save the second token (word):
    c=${b#\"}   # Remove the initial "
    d=${c%\"}   # Remove the final "
    progname=${d#astr} # Remove the initial "astr"

    echo astr$progname $progversion

    # Save the result to a file:
    echo @set ${progname^^}_VERSION $progversion >> $tmp
done

echo

#Replace the final file with the temporary file:
mv $tmp $output
