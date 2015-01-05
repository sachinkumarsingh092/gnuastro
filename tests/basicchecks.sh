# This file is included in all the tests to check if the executable
# was actually made or it the defaults file with the default value for
# the program that is to be tested exists.
#
# This file is not run by its self.


echo `pwd`


# If the executable was not made (the user chose to not install this
# package), skip this test:
if [ ! -f $execname ]; then
    exit 77
fi





# Check if the defaults directory has been correctly made, if not,
# make it and put the appropriate default file in it.
if [ ! -d .astrutils ]; then
    mkdir .astrutils
fi
if [ ! -f .astrutils/astr$prog.conf ]; then
    cp $topsrc/src/$prog/astr$prog.conf .astrutils/
fi
