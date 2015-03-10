# Process this file with autoconf to produce a configure script.
#
# Hand written file: used as input into GNU Autotools (autoconf).
# This is part of GNU Astronomy Utilities (gnuastro) package.
#
# Copyright (C) 2013-2015 Mohammad Akhlaghi
# Tohoku University Astronomical Institute, Sendai, Japan.
# http://astr.tohoku.ac.jp/~akhlaghi/
#
# gnuastro is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# gnuastro is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with gnuastro. If not, see <http://www.gnu.org/licenses/>.

# Definitions:
AC_PREREQ([2.69])
AC_INIT([GNU Astronomy Utilities], [0.1], [bug-gnuastro@gnu.org],
        [gnuastro], [http://www.gnu.org/software/gnuastro/])
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects gnu])
AC_CONFIG_SRCDIR([lib/fitsarrayvv.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])




# Checks for programs.
: ${CFLAGS=""}
AC_PROG_CC
gl_EARLY
AM_PROG_AR
LT_INIT




# Checks for libraries. The `absolutejunk=1` is put here so it doesn't
# automatically add the library to the LIBS variable that will be
# appended to all programs, some libraries will not be needed in some
# programs.
AC_CHECK_LIB([m], [sqrt], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([pthread], [pthread_attr_init], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])])
AC_CHECK_LIB([cfitsio], [ffopen], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lpthread -lm])
AC_CHECK_LIB([wcs], [wcspih], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lcfitsio -lpthread -lm])
AC_CHECK_LIB([jpeg], [jpeg_stdio_dest], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lm])
AC_CHECK_LIB([gslcblas], [cblas_sdsdot], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lm])
AC_CHECK_LIB([gsl], [gsl_integration_qng], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lgslcblas -lm])

AC_CHECK_LIB([pthread], [pthread_barrier_wait],
             [AC_DEFINE([HAVE_PTHREAD_BARRIER], [1], [Has pthread_barrier])])




#Gnulib checks:
gl_INIT
# Checks for header files.
#AC_CHECK_HEADERS([argp.h error.h float.h stdlib.h stdint.h string.h sys/time.h unistd.h],
#                 [], [AC_MSG_ERROR([Cannot continue.])])






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
             [AC_MSG_FAILURE([Value to --with-numthreads ($with_numthreads) is not a valid number])])],
      [
       AC_PATH_PROG(has_nproc, nproc, [no])
       AS_IF([test "x$has_nproc" = xno],
             [AS_ECHO([])
              AS_ECHO([])
              AC_MSG_FAILURE([The executable `nproc' was not found on your system. It is part of GNU Coreutils and specifies the number of CPU threads available for programs on your system. You can either install GNU Coreutils or simply run configure again with the option: `--with-numthreads=NUMBER-OF-THREADS' (replace NUMBER-OF-THREADS with any positive integer you like, ideally, it is double the number of physical CPU cores on your system). All the programs in GNU Astronomy Utilities will then use NUMBER-OF-THREADS if they need to run on multiple threads. Please see the manual for more information.])
              ])
       AC_MSG_CHECKING([for number of available threads])
       AC_DEFINE_UNQUOTED([DP_NUMTHREADS], [`nproc`])
       AC_MSG_RESULT([$(nproc)])
      ]
     )





#Check necessary program(s):
anywarnings=no
AC_CHECK_PROG(has_ghostscript, gs, [yes], [no])
AS_IF([test "x$has_ghostscript" = "xno"], anywarnings=yes)
AM_CONDITIONAL([COND_HASGHOSTSCRIPT], [test "x$has_ghostscript" = "xyes"])




# Check if Gnulib tests should be done:
AC_ARG_ENABLE([gnulibcheck],
              [AS_HELP_STRING([--enable-gnulibcheck],
  	            [In `make check', also test GNU Gnulib.])],
	      [enable_gnulibcheck=yes], [enable_gnulibcheck=no])
AM_CONDITIONAL([COND_GNULIBCHECK], [test $enable_gnulibcheck = yes])





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
AC_ARG_ENABLE([convertt],
              [AS_HELP_STRING([--enable-convertt],
  	            [Install ConvertType and other enabled packages only.])],
	      [AS_IF([test "x$enable_convertt" != xno],
                     [enable_convertt=yes; ayes=true])],
              [enable_convertt=notset])
AC_ARG_ENABLE([convolve],
              [AS_HELP_STRING([--enable-convolve],
  	            [Install Convolve and other enabled packages only.])],
	      [AS_IF([test "x$enable_convolve" != xno],
                     [enable_convolve=yes; ayes=true])],
              [enable_convolve=notset])
AC_ARG_ENABLE([imgcrop],
              [AS_HELP_STRING([--enable-imgcrop],
  	            [Install ImageCrop and other enabled packages only.])],
	      [AS_IF([test "x$enable_imgcrop" != xno],
                     [enable_imgcrop=yes; ayes=true])],
              [enable_imgcrop=notset])
AC_ARG_ENABLE([mkprof],
              [AS_HELP_STRING([--enable-mkprof],
  	            [Install MakeMock and other enabled packages only.])],
	      [AS_IF([test "x$enable_mkprof" != xno],
                     [enable_mkprof=yes; ayes=true])],
              [enable_mkprof=notset])




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
      [AS_IF([test $enable_convertt = notset], [enable_convertt=no])
       AS_IF([test $enable_convolve = notset], [enable_convolve=no])
       AS_IF([test $enable_imgcrop = notset], [enable_imgcrop=no])
       AS_IF([test $enable_mkprof = notset], [enable_mkprof=no])],

      [AS_IF([test $enable_convertt = notset], [enable_convertt=yes])
       AS_IF([test $enable_convolve = notset], [enable_convolve=yes])
       AS_IF([test $enable_imgcrop = notset], [enable_imgcrop=yes])
       AS_IF([test $enable_mkprof = notset], [enable_mkprof=yes])]
     )





# Make the enable_package values available for the Makefile:
AM_CONDITIONAL([COND_CONVERTT], [test $enable_convertt = yes])
AM_CONDITIONAL([COND_CONVOLVE], [test $enable_convolve = yes])
AM_CONDITIONAL([COND_IMGCROP], [test $enable_imgcrop = yes])
AM_CONDITIONAL([COND_MKPROF], [test $enable_mkprof = yes])


# Tell autoconf what to work on:
AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 lib/Makefile
		 tests/Makefile
                 gnulib/lib/Makefile
                 src/mkprof/Makefile
		 src/imgcrop/Makefile
                 gnulib/tests/Makefile
                 src/convertt/Makefile
                 src/convolve/Makefile
                 ])





#Make the outputs:
AC_OUTPUT





# Print a bye-bye message.
PACKAGE_FULL_NAME="AC_PACKAGE_NAME (AC_PACKAGE_TARNAME) AC_PACKAGE_VERSION"
AS_ECHO([])
AS_ECHO([===================================================================])
AS_ECHO([===================================================================])
AS_ECHO(["$PACKAGE_FULL_NAME, is configured."])
AS_ECHO([])
AS_IF([test "x$anywarnings" = xyes],
      [
        AS_ECHO([Warnings:])
        AS_IF([test "x$has_ghostscript" = xno],
              [AS_ECHO(["  - GPL GhostScript, with the executable name \`gs' was not found"])
               AS_ECHO(["      in your PATH environment variable. If PDF outputs are desired"])
               AS_ECHO(["      from ConvertType, it will stop at the EPS output and cannot "])
               AS_ECHO(["      proceed to generate the PDF."])])
        AS_ECHO([])
        AS_ECHO([])
      ]
     )
AS_ECHO(["Please run the following commands to build, test"])
AS_ECHO(["and install AC_PACKAGE_TARNAME AC_PACKAGE_VERSION in your machine."])
AS_ECHO(["    make"])
AS_ECHO(["    make check"])
AS_ECHO(["    make install"])
AS_ECHO([===================================================================])
AS_ECHO([===================================================================])
AS_ECHO([])
