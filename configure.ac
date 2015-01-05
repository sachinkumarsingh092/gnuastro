# Process this file with autoconf to produce a configure script.
#
# Hand written file: used as input into GNU Autotools (autoconf).
#
# This is the autoconf instructions to build NoiseChisel. NoiseChisel
# is a program to detect and segment objects heavily immersed in
# noise. It is free software released under GNU GPLv3+.
#
# Copyright (C) 2014 Mohammad Akhlaghi
# Tohoku University Astronomical Institute, Sendai, Japan.
# http://astr.tohoku.ac.jp/~akhlaghi/

# Definitions:
AC_PREREQ([2.69])
AC_INIT([GNU AstrUtils], [0.1], [makhlaghi@gmail.com],
        [astrutils], [http://www.gnu.org/software/astrutils/])
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects gnu])
AC_CONFIG_SRCDIR([lib/fitsarrayvv.c])
AC_CONFIG_HEADERS([config.h])





# Checks for programs.
: ${CFLAGS=""}
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_AR





# Checks for libraries. The `absolutejunk=1` is put here so it doesn't
# automatically add the library to the LIBS variable that will be
# appended to all programs, some libraries will not be needed in some
# programs.
AC_CHECK_LIB([m], [sqrt], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([pthread], [pthread_attr_init], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([cfitsio], [ffopen], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([wcs], [wcspih], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([jpeg], [jpeg_stdio_dest], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([fftw3f], [fftwf_plan_dft_c2r_2d], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([gslcblas], [cblas_sdsdot], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lm])
AC_CHECK_LIB([gsl], [gsl_integration_qng], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lgslcblas -lm])





# Checks for header files.
AC_CHECK_HEADERS([argp.h error.h float.h stdlib.h stdint.h string.h sys/time.h unistd.h],
                 [], [AC_MSG_ERROR([Cannot continue.])])





# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_UINT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T





# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([gettimeofday memset sqrt strtol mkdir strchr])





# To add the help2man:
AC_PATH_PROG(HELP2MAN, help2man)





# Read the configuration date.
AC_MSG_CHECKING([for date of compilation])
AC_DEFINE_UNQUOTED([CONFIGDATE], ["`date +"%d %B %Y"`"],
                   [Date of compilation])
AC_DEFINE_UNQUOTED([CONFIGTIME], ["`date +"%H:%M"`"],
                   [Time of compilation])
AC_MSG_RESULT([done])






# If the Number of threads to use was specified as a configure command
# line option, use it. If not, check the number of threads on the
# system.
AC_ARG_WITH([numthreads],
            [AS_HELP_STRING([--with-numthreads=INT],
	                    [Specify the default number of CPU threads.])],
	    [], [with_numthreads=no])
AH_TEMPLATE([DP_NUMTHREADS], [Number of threads to use.])
AS_IF([test "x$with_numthreads" != xno],
      [AS_IF([ [[[ $with_numthreads =~ ^[1-9]+$ ]]] ],
             [AC_DEFINE_UNQUOTED([DP_NUMTHREADS], [$with_numthreads])],
             [AC_MSG_FAILURE([Value to --with-readline ($with_numthreads) is not a valid number])])],
      [AC_MSG_CHECKING([for number of available threads])
       AC_DEFINE_UNQUOTED([DP_NUMTHREADS], [`nproc`])
       AC_MSG_RESULT([$(nproc)])]
     )





# Set the one general parameters:
AC_DEFINE_UNQUOTED([CONF_POSTFIX], [".conf"], [Configuration file post fix.])
AC_DEFINE_UNQUOTED([CURDIRCONFIG_DIR], [".AC_PACKAGE_TARNAME/"],
                   [Local data dir.])
AC_DEFINE_UNQUOTED([USERCONFIG_DIR], ["/.local/etc/"],
                   [User data dir.])
AC_DEFINE_UNQUOTED([CONF_SHOWFMT], [" %-20s"],
                   [Configuration file name format.])






# Read arguments about which programs to install. After checking if
# the argument was actually called, remove any value the user might
# have given by setting them to "yes" if they are not "no". These
# options don't accept arguments.
ayes=false
AC_ARG_ENABLE([mockgals],
              [AS_HELP_STRING([--enable-mockgals],
  	            [Install MockGals and other enabled packages only.])],
	      [AS_IF([test "x$enable_mockgals" != xno],
                     [enable_mockgals=yes; ayes=true])],
              [enable_mockgals=notset])
AC_ARG_ENABLE([imgcrop],
              [AS_HELP_STRING([--enable-imagecrop],
  	            [Install ImageCrop and other enabled packages only.])],
	      [AS_IF([test "x$enable_imgcrop" != xno],
                     [enable_imgcrop=yes; ayes=true])],
              [enable_imgcrop=notset])




# If we had a "ano" variable to be "true" if there was a no, then we
# would get. So we see we have no need for such a variable.
#
# if [ ano == true ]; then
#    every "notset" becomes yes.        if [ ayes == true ]; then
# elif [ ayes == true ]; then                 every "notset" becomes "no"
#    every "notset" becomes no.   ==>   else
# else                                        every "notset" becomes "yes"
#    every "notset"  becomes yes.       fi
# fi
AS_IF([test $ayes = true ],
      [AS_IF([test $enable_mockgals = notset], [enable_mockgals=no])
       AS_IF([test $enable_imgcrop = notset], [enable_imgcrop=no])],

      [AS_IF([test $enable_mockgals = notset], [enable_mockgals=yes])
       AS_IF([test $enable_imgcrop = notset], [enable_imgcrop=yes])]
     )





# Make the enable_package values available for the Makefile:
AM_CONDITIONAL([COND_MOCKGALS], [test $enable_mockgals = yes])
AM_CONDITIONAL([COND_IMGCROP], [test $enable_imgcrop = yes])


# Tell autoconf what to work on:
AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 lib/Makefile
		 tests/Makefile
		 src/imgcrop/Makefile
		 src/mockgals/Makefile])





#Make the outputs:
AC_OUTPUT





# Print a bye-bye message.
AS_ECHO([])
AS_ECHO([=====================================])
AS_ECHO(["AC_PACKAGE_STRING is configured."])
AS_ECHO([])
AS_ECHO(["Ready to build, check and install:"])
AS_ECHO(["    make"])
AS_ECHO(["    make check"])
AS_ECHO(["    make install"])
AS_ECHO([=====================================])
AS_ECHO([])
