# Process this file with autoconf to produce a configure script.
#
# Hand written file: used as input into GNU Autotools (autoconf).
# This is part of GNU Astronomy Utilities (gnuastro) package.
#
# Original author:
#     Mohammad Akhlaghi <akhlaghi@gnu.org>
# Contributing author(s):
# Copyright (C) 2015, Free Software Foundation, Inc.
#
# Gnuastro is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Gnuastro is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.





# Definitions:
AC_PREREQ([2.69])
AC_INIT([GNU Astronomy Utilities],
         m4_esyscmd([bootstrapped/build-aux/git-version-gen    \
                                  .tarball-version --prefix ""]),
        [bug-gnuastro@gnu.org], [gnuastro],
        [http://www.gnu.org/software/gnuastro/])
AC_CONFIG_AUX_DIR([bootstrapped/build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects gnu])
AC_CONFIG_SRCDIR([lib/fitsarrayvv.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([bootstrapped/m4])





# Checks for programs.
: ${CFLAGS=""}
AC_PROG_CC
gl_EARLY
AM_PROG_AR
LT_INIT





# Include Autoconf archives' ax_pthread.m4 check
m4_include([bootstrapped/m4/ax_pthread.m4])
AX_PTHREAD([],[AC_MSG_ERROR([AC_PACKAGE_NAME Needs POSIX Threads (pthread)])])
CLIBS="$PTHREAD_LIBS $LIBS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
LDFLAGS="$LDFLAGS $PTHREAD_CFLAGS"
CC="$PTHREAD_CC"





# Checks for necessary libraries. The `absolutejunk=1` is put here so
# it doesn't automatically add the library to the LIBS variable that
# will be appended to all programs, some libraries will not be needed
# in some programs.
AC_CHECK_LIB([cfitsio], [ffopen], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lm])
AC_CHECK_LIB([wcs], [wcspih], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lcfitsio -lm])
AC_CHECK_LIB([gslcblas], [cblas_sdsdot], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lm])
AC_CHECK_LIB([gsl], [gsl_integration_qng], [absolutejunk=1],
             [AC_MSG_ERROR([Cannot continue.])], [-lgslcblas -lm])

AC_CHECK_LIB([wcs], [wcslib_version],
             [AC_DEFINE([HAVE_WCSLIBVERSION], [1], [Has wcslib_version])],
             [], [-lcfitsio -lm])
AC_CHECK_LIB([pthread], [pthread_barrier_wait],
             [AC_DEFINE([HAVE_PTHREAD_BARRIER], [1], [Has pthread_barrier])])






#Gnulib checks:
gl_INIT






# To add the help2man:
AC_PATH_PROG(HELP2MAN, help2man)





# Read the configuration date.
AC_MSG_CHECKING([for date of compilation])
AC_DEFINE_UNQUOTED([CONFIGDATE], ["`date +"%d %B %Y"`"],
                   [Date of compilation])
AC_DEFINE_UNQUOTED([CONFIGTIME], ["`date +"%H:%M"`"],
                   [Time of compilation])
AC_MSG_RESULT([done])







# Check libjpeg:
anywarnings=no
AC_CHECK_LIB([jpeg], [jpeg_stdio_dest],
             [has_libjpeg=yes], [has_libjpeg=no], [-lm])
AS_IF([test "x$has_libjpeg" = "xyes"],
      [AC_DEFINE([HAS_LIBJPEG], [], [Has libjpeg])],
      [anywarnings=yes])
AM_CONDITIONAL([COND_HASLIBJPEG], [test "x$has_libjpeg" = "xyes"])





# Check Ghostscript: "-dPDFFitPage" option to Ghostscript, used by
# ConvertType to convert from EPS to PDF, has been introduced in
# Ghostscript 9.10.  Make sure we have at least that version.
#
# Below, only when Ghostscript exists, we check its version and only if its
# version is larger than 9.10, does Gnuastro finally assume the existence
# of Ghostscript.
AC_CHECK_PROG(has_ghostscript, gs, [yes], [no])
AS_IF([test "x$has_ghostscript" = "xyes"],
      [AC_MSG_CHECKING(Ghostscript version)
       case `(gs --version) 2>/dev/null` in
         [1-8].* | 9.0[0-9]) has_ghostscript=no;;
       esac
       AC_MSG_RESULT( $(gs --version) )])
# Note: `has_ghostscript' can be set to `no' within the AS_IF above, so
# `anywarnings' cannot be an [RUN-IF-FALSE] argument to the AS_IF above.
AS_IF([test "x$has_ghostscript" = "xno"], [anywarnings=yes])
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
                     [enable_cognvolve=yes; ayes=true])],
              [enable_convolve=notset])
AC_ARG_ENABLE([cosmiccal],
              [AS_HELP_STRING([--enable-cosmiccal],
  	            [Install CosmicCalculator and other enabled packages only.])],
	      [AS_IF([test "x$enable_cosmiccal" != xno],
                     [enable_cosmiccal=yes; ayes=true])],
              [enable_cosmiccal=notset])
AC_ARG_ENABLE([header],
              [AS_HELP_STRING([--enable-header],
  	            [Install Header and other enabled packages only.])],
	      [AS_IF([test "x$enable_header" != xno],
                     [enable_header=yes; ayes=true])],
              [enable_header=notset])
AC_ARG_ENABLE([arithmetic],
              [AS_HELP_STRING([--enable-arithmetic],
  	            [Install ImageArithmetic and other enabled packages only.])],
	      [AS_IF([test "x$enable_arithmetic" != xno],
                     [enable_arithmetic=yes; ayes=true])],
              [enable_arithmetic=notset])
AC_ARG_ENABLE([imgcrop],
              [AS_HELP_STRING([--enable-imgcrop],
  	            [Install ImageCrop and other enabled packages only.])],
	      [AS_IF([test "x$enable_imgcrop" != xno],
                     [enable_imgcrop=yes; ayes=true])],
              [enable_imgcrop=notset])
AC_ARG_ENABLE([imgstat],
              [AS_HELP_STRING([--enable-imgstat],
  	            [Install ImageStatistics and other enabled packages only.])],
	      [AS_IF([test "x$enable_imgstat" != xno],
                     [enable_imgstat=yes; ayes=true])],
              [enable_imgstat=notset])
AC_ARG_ENABLE([imgwarp],
              [AS_HELP_STRING([--enable-imgwarp],
  	            [Install ImageWarp and other enabled packages only.])],
	      [AS_IF([test "x$enable_imgwarp" != xno],
                     [enable_imgwarp=yes; ayes=true])],
              [enable_imgwarp=notset])
AC_ARG_ENABLE([mkcatalog],
              [AS_HELP_STRING([--enable-mkcatalog],
  	            [Install MakeCatalog and other enabled packages only.])],
	      [AS_IF([test "x$enable_mkcatalog" != xno],
                     [enable_mkcatalog=yes; ayes=true])],
              [enable_mkcatalog=notset])
AC_ARG_ENABLE([mknoise],
              [AS_HELP_STRING([--enable-mknoise],
  	            [Install MakeNoise and other enabled packages only.])],
	      [AS_IF([test "x$enable_mknoise" != xno],
                     [enable_mknoise=yes; ayes=true])],
              [enable_mknoise=notset])
AC_ARG_ENABLE([mkprof],
              [AS_HELP_STRING([--enable-mkprof],
  	            [Install MakeProfile and other enabled packages only.])],
	      [AS_IF([test "x$enable_mkprof" != xno],
                     [enable_mkprof=yes; ayes=true])],
              [enable_mkprof=notset])
AC_ARG_ENABLE([noisechisel],
              [AS_HELP_STRING([--enable-noisechisel],
  	            [Install NoiseChisel and other enabled packages only.])],
	      [AS_IF([test "x$enable_noisechisel" != xno],
                     [enable_noisechisel=yes; ayes=true])],
              [enable_noisechisel=notset])
AC_ARG_ENABLE([subtractsky],
              [AS_HELP_STRING([--enable-subtractsky],
  	            [Install SubtractSky and other enabled packages only.])],
	      [AS_IF([test "x$enable_subtractsky" != xno],
                     [enable_subtractsky=yes; ayes=true])],
              [enable_subtractsky=notset])




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
       AS_IF([test $enable_cosmiccal = notset], [enable_cosmiccal=no])
       AS_IF([test $enable_header = notset], [enable_header=no])
       AS_IF([test $enable_arithmetic = notset], [enable_arithmetic=no])
       AS_IF([test $enable_imgcrop = notset], [enable_imgcrop=no])
       AS_IF([test $enable_imgstat = notset], [enable_imgstat=no])
       AS_IF([test $enable_imgwarp = notset], [enable_imgwarp=no])
       AS_IF([test $enable_mkcatalog = notset], [enable_mkcatalog=no])
       AS_IF([test $enable_mknoise = notset], [enable_mknoise=no])
       AS_IF([test $enable_mkprof = notset], [enable_mkprof=no])
       AS_IF([test $enable_noisechisel = notset], [enable_noisechisel=no])
       AS_IF([test $enable_subtractsky = notset], [enable_subtractsky=no])
       ],

      [AS_IF([test $enable_convertt = notset], [enable_convertt=yes])
       AS_IF([test $enable_convolve = notset], [enable_convolve=yes])
       AS_IF([test $enable_cosmiccal = notset], [enable_cosmiccal=yes])
       AS_IF([test $enable_header = notset], [enable_header=yes])
       AS_IF([test $enable_arithmetic = notset], [enable_arithmetic=yes])
       AS_IF([test $enable_imgcrop = notset], [enable_imgcrop=yes])
       AS_IF([test $enable_imgstat = notset], [enable_imgstat=yes])
       AS_IF([test $enable_imgwarp = notset], [enable_imgwarp=yes])
       AS_IF([test $enable_mkcatalog = notset], [enable_mkcatalog=yes])
       AS_IF([test $enable_mknoise = notset], [enable_mknoise=yes])
       AS_IF([test $enable_mkprof = notset], [enable_mkprof=yes])
       AS_IF([test $enable_noisechisel = notset], [enable_noisechisel=yes])
       AS_IF([test $enable_subtractsky = notset], [enable_subtractsky=yes])
       ]
     )





# Make the enable_package values available for the Makefile:
AM_CONDITIONAL([COND_CONVERTT], [test $enable_convertt = yes])
AM_CONDITIONAL([COND_CONVOLVE], [test $enable_convolve = yes])
AM_CONDITIONAL([COND_COSMICCAL], [test $enable_cosmiccal = yes])
AM_CONDITIONAL([COND_HEADER], [test $enable_header = yes])
AM_CONDITIONAL([COND_IMGCROP], [test $enable_imgcrop = yes])
AM_CONDITIONAL([COND_ARITHMETIC], [test $enable_arithmetic = yes])
AM_CONDITIONAL([COND_IMGSTAT], [test $enable_imgstat = yes])
AM_CONDITIONAL([COND_IMGWARP], [test $enable_imgwarp = yes])
AM_CONDITIONAL([COND_MKCATALOG], [test $enable_mkcatalog = yes])
AM_CONDITIONAL([COND_MKNOISE], [test $enable_mknoise = yes])
AM_CONDITIONAL([COND_MKPROF], [test $enable_mkprof = yes])
AM_CONDITIONAL([COND_NOISECHISEL], [test $enable_noisechisel = yes])
AM_CONDITIONAL([COND_SUBTRACTSKY], [test $enable_subtractsky = yes])





# Tell autoconf what to work on:
AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 lib/Makefile
		 tests/Makefile
                 src/mkprof/Makefile
                 src/header/Makefile
                 src/mknoise/Makefile
		 src/imgcrop/Makefile
           	 src/imgstat/Makefile
       		 src/imgwarp/Makefile
                 src/convertt/Makefile
                 src/convolve/Makefile
                 src/cosmiccal/Makefile
           	 src/mkcatalog/Makefile
                 src/arithmetic/Makefile
                 src/noisechisel/Makefile
                 src/subtractsky/Makefile
                 bootstrapped/lib/Makefile
                 bootstrapped/tests/Makefile
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
        AS_IF([test "x$has_libjpeg" = "xno"],
              [AS_ECHO(["  - libjpeg, could not be linked with in your library search path."])
               AS_ECHO(["      If JPEG outputs are desired from ConvertType, it will warn"])
               AS_ECHO(["      you and abort with an error."])])
        AS_IF([test "x$has_ghostscript" = "xno"],
              [AS_ECHO(["  - GPL GhostScript version 9.10 or later, with the executable"])
               AS_ECHO(["      name \`gs', was not found in your PATH environment variable."])
               AS_ECHO(["      If PDF outputs are desired from ConvertType, it will stop"])
               AS_ECHO(["      at the EPS output and cannot proceed to generate the PDF."])])
        AS_ECHO(["  All checks related to the above warnings will be skipped."])
        AS_ECHO([])
        AS_ECHO([])
      ]
     )
AS_ECHO(["Please run the following commands to build, test"])
AS_ECHO(["and install AC_PACKAGE_TARNAME AC_PACKAGE_VERSION on your machine."])
AS_ECHO(["    make"])
AS_ECHO(["    make check"])
AS_ECHO(["    make install"])
AS_ECHO([===================================================================])
AS_ECHO([===================================================================])
AS_ECHO([])
