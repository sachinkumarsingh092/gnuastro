# This is not actually a test of any of the programs, it just brings
# in all the configuration files into a locally created .gnuastro
# directory for all the tests to use.
#
# Original author:
#     Mohammad Akhlaghi <akhlaghi@gnu.org>
# Contributing author(s):
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.


# Make the .gnuastro directory if not already created.
if [ ! -d .gnuastro ]; then
    mkdir .gnuastro
fi


#For each program bring in the configuration file:
for prog in imgcrop mkprof convertt convolve
do
    if [ ! -f .gnuastro/ast$prog.conf ]; then
	cp $topsrc/src/$prog/ast$prog.conf .gnuastro/
    fi
done
