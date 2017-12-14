# Process this file with autoconf to produce a configure script.
#
# Hand written file: used as input into GNU Autotools (autoconf).
# This is part of GNU Astronomy Utilities (gnuastro) package.
#
# We use some macros that are distributed as part of the GNU Autoconf
# Archive here. Those macros are all placed in the bootstrapped/m4/
# directory along with those that were imported from Gnulib. The
# ./bootstrap script includes that directory while its running Autoconf so
# we don't need to explicitly call that directory here.
#
# Original author:
#     Mohammad Akhlaghi <akhlaghi@gnu.org>
# Contributing author(s):
# Copyright (C) 2015 - 2016, Free Software Foundation, Inc.
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
                        .tarball-version --prefix "gnuastro_v"]),
        [bug-gnuastro@gnu.org], [gnuastro],
        [http://www.gnu.org/software/gnuastro/])
AC_CONFIG_AUX_DIR([bootstrapped/build-aux])
AM_INIT_AUTOMAKE([-Wall subdir-objects gnu])
AC_CONFIG_SRCDIR([lib/fits.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([bootstrapped/m4])





# Library version, see the GNU Libtool manual ("Library interface versions"
# section for the exact definition of each) for
GAL_CURRENT=3
GAL_REVISION=0
GAL_AGE=0
GAL_LT_VERSION="${GAL_CURRENT}:${GAL_REVISION}:${GAL_AGE}"
AC_SUBST(GAL_LT_VERSION)





# Checks for programs.
: ${CFLAGS=""}
AC_PROG_CC
AC_PROG_CXX
gl_EARLY
AM_PROG_AR
LT_INIT





# This macro will let the libraries know that we are now in the Gnuastro
# build system, not on the user's system. While we are building Gnuastro,
# we have the important installation information in `config.h'. But in the
# user's own programs, this information is defined in
# `gnuastro/config.h'. With this macro, the installed headers can decide
# if the latter should be included or not. Note that `gnuastro/config.h'
# is only built at installation time and doesn't exist when building
# Gnuastro. Therefore, this macro must not be defined in a user's program.
AC_DEFINE([IN_GNUASTRO_BUILD], [1], [In building, not usage])





# Generic compiler flags for all sub-directories.
CFLAGS="-Wall -O3 $CFLAGS"
LDFLAGS="-L\$(top_builddir)/lib $LDFLAGS"
CPPFLAGS="-I\$(top_srcdir)/lib $CPPFLAGS"





# Check for pthreads and add the appropriate compilation flags. AX_PTHREAD
# comes from the GNU Autoconf Archive's ax_pthread.m4.
AX_PTHREAD([],[AC_MSG_ERROR([AC_PACKAGE_NAME Needs POSIX Threads (pthread)])])
CLIBS="$PTHREAD_LIBS $LIBS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
LDFLAGS="$LDFLAGS $PTHREAD_CFLAGS"
CC="$PTHREAD_CC"





# See if the C++ compiler understands `-Qunused-arguments'. AX_PTHREAD adds
# this option to LDFLAGS when the C compiler knows this option. But as
# reported in bug #52490, it can happen that at the same time, the C++
# compiler doesn't recognize it.
cpp_Qunused_arguments_problem=no
for flg in $LDFLAGS; do
  AS_IF([test "$flg" = \-Qunused-arguments],
        [ AC_LANG(C++)
          AX_CHECK_COMPILE_FLAG([-Qunused-arguments],
                                [cpp_Qunused_arguments_problem=no],
                                [cpp_Qunused_arguments_problem=yes])
          AC_LANG(C) ])
done
AM_CONDITIONAL([COND_NO_CPP_QUNUSED_ARGUMENTS_PROBLEM],
               [test "x$cpp_Qunused_arguments_problem" = "xno"])





# Check if `malloc(0)' returns valid pointer
AC_FUNC_MALLOC





# Check the size of necessary system specific types.
AC_CHECK_SIZEOF([size_t])
AC_SUBST(SIZEOF_SIZE_T, [$ac_cv_sizeof_size_t])
AC_DEFINE_UNQUOTED([GAL_CONFIG_SIZEOF_SIZE_T], [$ac_cv_sizeof_size_t],
                   [On 32bit will be 4, on 64 bit, will be 8])
AC_CHECK_SIZEOF([long])
AC_SUBST(SIZEOF_LONG, [$ac_cv_sizeof_long])
AC_DEFINE_UNQUOTED([GAL_CONFIG_SIZEOF_LONG], [$ac_cv_sizeof_long],
                   [On 32bit will be 4, on 64 bit, will be 8])





# By default we assume no warnings
anywarnings=no





# Remove any occurance of the current directory `./', `.', or the full
# address of the current directory in PATH. The main problem is the
# `libtool' executable which Gnuastro builds internally in the top build
# directory. However, we also need to know if the system has libtool or
# not.
AC_MSG_CHECKING(if PATH contains current directory)
oldPATH=$PATH
currpwd=$(pwd)

# The first call to SED will remove any occurance of the current directory:
# `./', `.', or the full address.
#
#    NOTE 1: We cannot simply remove all `.'s, because hidden directories
#            (like the `~/.local' that is suggested for local
#            installations) will also be altered.
#
#    NOTE 2: An empty string in the list of strings (separated by `:')
#            means the current directory. This includes cases like: `::',
#            or a leading and trailing `:'. So after all the removals of
#            the current directory, we will remove all such cases.
#
#    NOTE 3: The SED separator can be any character immediately after `s',
#            it doesn't just have to be the commonly used `/'. Since `$pwd'
#            will possibly contain many `/'s, it is much more easier to use
#            a differen separator (`|' in this call to SED).
PATH=$(AS_ECHO([$PATH]) | $SED -e 's|'"$currpwd"'||g' \
                               -e 's|\.\.*//*||g'     \
                               -e 's|:\.\.*:|:|g'     \
                               -e 's|\.*$||'          \
                               -e 's|^\.*||'          \
                               -e 's|::*|:|g'         \
                               -e 's|^:||'            \
                               -e 's|:$||' )
AS_IF([test $oldPATH = $PATH],
      [ path_warning=no ],
      [ path_warning=yes; anywarnings=yes ])
AC_MSG_RESULT( $path_warning )





# Search for necessary libraries. After each library is found,
# AC_SEARCH_LIBS adds the -lLIBRARY flag to the LIBS variable which is then
# given to all the Makefiles. Each new flag is added to the left of the old
# one so order matters here. Note that the LIBS variable is also used in
# checking the next libraries, so the linking with their dependent
# libraries is done automatically with this order, and we don't have to
# explicitly set the dependency flags.
AC_SEARCH_LIBS(sqrt, m, [],
    [AC_MSG_ERROR([C math library not present, cannot continue.])])
AC_SEARCH_LIBS([cblas_sdsdot], [gslcblas], [],
    [AC_MSG_ERROR([GSL CBLAS not present, cannot continue.])])
AC_SEARCH_LIBS([gsl_integration_qng], [gsl], [],
    [AC_MSG_ERROR([GSL not found, cannot continue.])])

# Since version 0.42, if `libcurl' is installed, CFITSIO will link with it
# and thus it will be necessary to explicitly link with libcurl also. If it
# doesn't exist on the system, then CFITSIO won't link with it and there is
# no problem for Gnuastro either. So there is no need to stop the configure
# script if libcurl isn't found.
AC_SEARCH_LIBS([curl_global_init], [curl], [], [])
AC_SEARCH_LIBS([ffopen], [cfitsio], [],
    [AC_MSG_ERROR([CFITSIO not found, cannot continue.])])
AC_SEARCH_LIBS([wcspih], [wcs], [],
    [AC_MSG_ERROR([WCSLIB not found, cannot continue.])])

# These are secondary tests for more fine-grained control in libraries that
# have already been checked. We don't need to add them to the LIBS
# variable, so we are using AC_CHECK_LIB for these tests.
AC_CHECK_LIB([wcs], [wcslib_version], [has_wcslib_version=1],
             [has_wcslib_version=0], [-lcfitsio -lm])
AC_DEFINE_UNQUOTED([GAL_CONFIG_HAVE_WCSLIB_VERSION], [$has_wcslib_version],
                   [WCSLIB comes with wcslib_version])
AC_SUBST(HAVE_WCSLIB_VERSION, [$has_wcslib_version])

AC_CHECK_LIB([pthread], [pthread_barrier_wait], [has_pthread_barrier=1],
             [has_pthread_barrier=0])
AC_DEFINE_UNQUOTED([GAL_CONFIG_HAVE_PTHREAD_BARRIER], [$has_pthread_barrier],
                   [System has pthread_barrier])
AC_SUBST(HAVE_PTHREAD_BARRIER, [$has_pthread_barrier])






# Gnulib checks:
gl_INIT




# Gnulib checks for the proper name for the C99 equivalent `restrict'
# keyword and puts it in the `ac_cv_c_restrict' variable. If none exists,
# it will put a `no' inside of this variable. As described in the output
# `bootstrapped/m4/gnulib-common.m4', this is only necessary until Autoconf
# 2.70 is released. Afterwards, we can use AC_C_RESTRICT.
AS_IF([test "x$ac_cv_c_restrict" = "xno"],
      [gal_restrict_replace=], [gal_restrict_replace=$ac_cv_c_restrict])
AC_SUBST(RESTRICT_REPLACEMENT, [$gal_restrict_replace])





# To add the help2man:
AC_CHECK_PROG(has_help2man, help2man, [yes], [no])
AM_CONDITIONAL([COND_HASHELP2MAN], [test "x$has_help2man" = "xyes"])







# Check libjpeg:
AC_SEARCH_LIBS([jpeg_stdio_dest], [jpeg],
               [has_libjpeg=yes], [has_libjpeg=no])
AS_IF([test "x$has_libjpeg" = "xyes"],
      [AC_DEFINE([HAVE_LIBJPEG], [], [Has libjpeg])],
      [anywarnings=yes])
AM_CONDITIONAL([COND_HASLIBJPEG], [test "x$has_libjpeg" = "xyes"])





# Check libgit2:
AC_SEARCH_LIBS([git_libgit2_init], [git2],
               [has_libgit2=1], [has_libgit2=0])
AC_DEFINE_UNQUOTED([GAL_CONFIG_HAVE_LIBGIT2], [$has_libgit2],
                   [libgit2 is installed on the system])
AS_IF([test "x$has_libgit2" = "x1"], [], [anywarnings=yes])
AC_SUBST(HAVE_LIBGIT2, [$has_libgit2])





# Check the libtool executable on the system. Note that Gnuastro also ships
# with a version of Libtool. We don't want Gnuastro's Libtool, here we want
# to see if the system has libtool independent of Gnuastro so BuildProgram
# can use it later.
AC_CHECK_PROG(has_libtool, libtool, [yes], [no])

# If Libtool is present, make sure it is GNU Libtool
AS_IF([test "x$has_libtool" = "xyes"],
      [ AC_MSG_CHECKING(if libtool executable is GNU)
        AS_IF( libtool --version 2> /dev/null | grep GNU 2>&1 > /dev/null,
               [has_gnulibtool=yes], [has_gnulibtool=no])
        AC_MSG_RESULT( $has_gnulibtool )],
      [ has_gnulibtool=no ])

# When either the `libtool' executable isn't GNU or it doesn't exist, then
# look for `glibtool'.
AS_IF([test "x$has_gnulibtool" = "xyes"],
      [ gnulibtool_exec=libtool ],
      [ AC_CHECK_PROG(has_glibtool, glibtool, [yes], [no])
        AS_IF( [test "x$has_glibtool" = "xyes"],
               [has_gnulibtool=yes; gnulibtool_exec=glibtool],
               [has_gnulibtool=no; anywarnings=yes] ) ])

# Write the test results in the Makefiles and `config.h'.
AM_CONDITIONAL([COND_HASGNULIBTOOL], [test "x$has_gnulibtool" = "xyes"])
AC_DEFINE_UNQUOTED([GAL_CONFIG_GNULIBTOOL_EXEC], ["$gnulibtool_exec"],
                   [The executable to call GNU Libtool])




# Check Ghostscript: "-dPDFFitPage" option to Ghostscript, used by
# ConvertType to convert from EPS to PDF, has been introduced in
# Ghostscript 9.10.  Make sure we have at least that version.
#
# Below, only when Ghostscript exists, we check its version and only if its
# version is larger than 9.10, does Gnuastro finally assume the existence
# of Ghostscript. AX_COMPARE_VERSION comes from the GNU Autoconf Archive's
# ax_compare_version.m4.
AC_CHECK_PROG(has_ghostscript, gs, [yes], [no])
AS_IF([test "x$has_ghostscript" = "xyes"],
      [AC_MSG_CHECKING(Ghostscript version)
       gsversion=$(gs --version)
       AX_COMPARE_VERSION([9.10], [gt], [$gsversion], [has_ghostscript=no])
       AC_MSG_RESULT( $gsversion )])
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
AC_DEFINE_UNQUOTED([USERCONFIG_DIR], [".local/etc"], [User data dir.])
AC_DEFINE_UNQUOTED([CONF_SHOWFMT], [" %-20s"],
                   [Configuration file name format.])





# Read arguments about which programs to install. After checking if
# the argument was actually called, remove any value the user might
# have given by setting them to "yes" if they are not "no". These
# options don't accept arguments.
ayes=false
AC_ARG_ENABLE([arithmetic],
              [AS_HELP_STRING([--enable-arithmetic],
                    [Install Arithmetic and other enabled programs.])],
	      [AS_IF([test "x$enable_arithmetic" != xno],
                     [enable_arithmetic=yes; ayes=true])],
              [enable_arithmetic=notset])
AC_ARG_ENABLE([buildprog],
              [AS_HELP_STRING([--enable-buildprog],
                    [Install BuildProgram and other enabled programs.])],
	      [AS_IF([test "x$enable_buildprog" != xno],
                     [enable_buildprog=yes; ayes=true])],
              [enable_buildprog=notset])
AC_ARG_ENABLE([convertt],
              [AS_HELP_STRING([--enable-convertt],
                    [Install ConvertType and other enabled programs.])],
	      [AS_IF([test "x$enable_convertt" != xno],
                     [enable_convertt=yes; ayes=true])],
              [enable_convertt=notset])
AC_ARG_ENABLE([convolve],
              [AS_HELP_STRING([--enable-convolve],
                    [Install Convolve and other enabled programs.])],
	      [AS_IF([test "x$enable_convolve" != xno],
                     [enable_cognvolve=yes; ayes=true])],
              [enable_convolve=notset])
AC_ARG_ENABLE([cosmiccal],
              [AS_HELP_STRING([--enable-cosmiccal],
                    [Install CosmicCalculator and other enabled programs.])],
	      [AS_IF([test "x$enable_cosmiccal" != xno],
                     [enable_cosmiccal=yes; ayes=true])],
              [enable_cosmiccal=notset])
AC_ARG_ENABLE([crop],
              [AS_HELP_STRING([--enable-crop],
                    [Install Crop and other enabled programs.])],
	      [AS_IF([test "x$enable_crop" != xno],
                     [enable_crop=yes; ayes=true])],
              [enable_crop=notset])
AC_ARG_ENABLE([fits],
              [AS_HELP_STRING([--enable-fits],
                    [Install Fits and other enabled programs.])],
	      [AS_IF([test "x$enable_fits" != xno],
                     [enable_fits=yes; ayes=true])],
              [enable_fits=notset])
AC_ARG_ENABLE([match],
              [AS_HELP_STRING([--enable-match],
                    [Install Match and other enabled programs.])],
	      [AS_IF([test "x$enable_match" != xno],
                     [enable_match=yes; ayes=true])],
              [enable_match=notset])
AC_ARG_ENABLE([mkcatalog],
              [AS_HELP_STRING([--enable-mkcatalog],
                    [Install MakeCatalog and other enabled programs.])],
	      [AS_IF([test "x$enable_mkcatalog" != xno],
                     [enable_mkcatalog=yes; ayes=true])],
              [enable_mkcatalog=notset])
AC_ARG_ENABLE([mknoise],
              [AS_HELP_STRING([--enable-mknoise],
                    [Install MakeNoise and other enabled programs.])],
	      [AS_IF([test "x$enable_mknoise" != xno],
                     [enable_mknoise=yes; ayes=true])],
              [enable_mknoise=notset])
AC_ARG_ENABLE([mkprof],
              [AS_HELP_STRING([--enable-mkprof],
                    [Install MakeProfile and other enabled programs.])],
	      [AS_IF([test "x$enable_mkprof" != xno],
                     [enable_mkprof=yes; ayes=true])],
              [enable_mkprof=notset])
AC_ARG_ENABLE([noisechisel],
              [AS_HELP_STRING([--enable-noisechisel],
                    [Install NoiseChisel and other enabled programs.])],
	      [AS_IF([test "x$enable_noisechisel" != xno],
                     [enable_noisechisel=yes; ayes=true])],
              [enable_noisechisel=notset])
AC_ARG_ENABLE([statistics],
              [AS_HELP_STRING([--enable-statistics],
                    [Install Statistics and other enabled programs.])],
	      [AS_IF([test "x$enable_statistics" != xno],
                     [enable_statistics=yes; ayes=true])],
              [enable_statistics=notset])
AC_ARG_ENABLE([table],
              [AS_HELP_STRING([--enable-table],
                    [Install Table and other enabled programs.])],
	      [AS_IF([test "x$enable_table" != xno],
                     [enable_table=yes; ayes=true])],
              [enable_table=notset])
#AC_ARG_ENABLE([TEMPLATE],
#              [AS_HELP_STRING([--enable-TEMPLATE],
#                    [Install TEMPLATE and other enabled packages.])],
#	      [AS_IF([test "x$enable_TEMPLATE" != xno],
#                     [enable_TEMPLATE=yes; ayes=true])],
#              [enable_TEMPLATE=notset])
AC_ARG_ENABLE([warp],
              [AS_HELP_STRING([--enable-warp],
                    [Install Warp and other enabled programs.])],
	      [AS_IF([test "x$enable_warp" != xno],
                     [enable_warp=yes; ayes=true])],
              [enable_warp=notset])





# If we had a "ayes" variable to be "true" if there was a `yes'. So any
# program that is not explicitly requested must be ignored and vice versa
# (if no programs were explicitly requested, then enable all that weren't
# disabled).
AS_IF([test $ayes = true ],
      [
       AS_IF([test $enable_arithmetic = notset], [enable_arithmetic=no])
       AS_IF([test $enable_buildprog = notset], [enable_buildprog=no])
       AS_IF([test $enable_convertt = notset], [enable_convertt=no])
       AS_IF([test $enable_convolve = notset], [enable_convolve=no])
       AS_IF([test $enable_cosmiccal = notset], [enable_cosmiccal=no])
       AS_IF([test $enable_crop = notset], [enable_crop=no])
       AS_IF([test $enable_fits = notset], [enable_fits=no])
       AS_IF([test $enable_match = notset], [enable_match=no])
       AS_IF([test $enable_mkcatalog = notset], [enable_mkcatalog=no])
       AS_IF([test $enable_mknoise = notset], [enable_mknoise=no])
       AS_IF([test $enable_mkprof = notset], [enable_mkprof=no])
       AS_IF([test $enable_noisechisel = notset], [enable_noisechisel=no])
       AS_IF([test $enable_statistics = notset], [enable_statistics=no])
       AS_IF([test $enable_table = notset], [enable_table=no])
#      AS_IF([test $enable_TEMPLATE = notset], [enable_TEMPLATE=no])
       AS_IF([test $enable_warp = notset], [enable_warp=no])
       ],

      [
       AS_IF([test $enable_arithmetic = notset], [enable_arithmetic=yes])
       AS_IF([test $enable_buildprog = notset], [enable_buildprog=yes])
       AS_IF([test $enable_convertt = notset], [enable_convertt=yes])
       AS_IF([test $enable_convolve = notset], [enable_convolve=yes])
       AS_IF([test $enable_cosmiccal = notset], [enable_cosmiccal=yes])
       AS_IF([test $enable_crop = notset], [enable_crop=yes])
       AS_IF([test $enable_fits = notset], [enable_fits=yes])
       AS_IF([test $enable_match = notset], [enable_match=yes])
       AS_IF([test $enable_mkcatalog = notset], [enable_mkcatalog=yes])
       AS_IF([test $enable_mknoise = notset], [enable_mknoise=yes])
       AS_IF([test $enable_mkprof = notset], [enable_mkprof=yes])
       AS_IF([test $enable_noisechisel = notset], [enable_noisechisel=yes])
       AS_IF([test $enable_statistics = notset], [enable_statistics=yes])
       AS_IF([test $enable_table = notset], [enable_table=yes])
#      AS_IF([test $enable_TEMPLATE = notset], [enable_TEMPLATE=yes])
       AS_IF([test $enable_warp = notset], [enable_warp=yes])
       ]
     )





# BuildProgram depends on the presence of GNU Libtool, if it isn't present,
# then don't build it.
AS_IF([test "x$has_gnulibtool" = "xno"], [enable_buildprog=no])





# Make the enable_package values available for the Makefile:
AM_CONDITIONAL([COND_ARITHMETIC], [test $enable_arithmetic = yes])
AM_CONDITIONAL([COND_BUILDPROG], [test $enable_buildprog = yes])
AM_CONDITIONAL([COND_CONVERTT], [test $enable_convertt = yes])
AM_CONDITIONAL([COND_CONVOLVE], [test $enable_convolve = yes])
AM_CONDITIONAL([COND_COSMICCAL], [test $enable_cosmiccal = yes])
AM_CONDITIONAL([COND_CROP], [test $enable_crop = yes])
AM_CONDITIONAL([COND_FITS], [test $enable_fits = yes])
AM_CONDITIONAL([COND_MATCH], [test $enable_match = yes])
AM_CONDITIONAL([COND_MKCATALOG], [test $enable_mkcatalog = yes])
AM_CONDITIONAL([COND_MKNOISE], [test $enable_mknoise = yes])
AM_CONDITIONAL([COND_MKPROF], [test $enable_mkprof = yes])
AM_CONDITIONAL([COND_NOISECHISEL], [test $enable_noisechisel = yes])
AM_CONDITIONAL([COND_STATISTICS], [test $enable_statistics = yes])
AM_CONDITIONAL([COND_TABLE], [test $enable_table = yes])
#AM_CONDITIONAL([COND_TEMPLATE], [test $enable_TEMPLATE = yes])
AM_CONDITIONAL([COND_WARP], [test $enable_warp = yes])





# Tell autoconf what to work on: TEMPLATE cannot be put and then
# commented here like the cases above, so don't forget to add your new
# utility name here.
AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 lib/Makefile
		 tests/Makefile
                 bin/crop/Makefile
                 bin/fits/Makefile
                 bin/warp/Makefile
                 bin/table/Makefile
                 bin/match/Makefile
                 bin/mkprof/Makefile
                 bin/mknoise/Makefile
                 bin/convertt/Makefile
                 bin/convolve/Makefile
                 bin/buildprog/Makefile
                 bin/cosmiccal/Makefile
           	 bin/mkcatalog/Makefile
                 bin/arithmetic/Makefile
           	 bin/statistics/Makefile
                 bin/noisechisel/Makefile
                 bootstrapped/lib/Makefile
                 bootstrapped/tests/Makefile
                 ])





# Printing guiding messages. Autoconf will make the variable
# enable_guide_message from the first argument to AC_ARG_ENABLE. It will
# also give it a value. From the Autoconf manual, we see that
# `--disable-guide-message' is equivalent to a value of `no', while with no
# argument, the value will default to `yes'. In the last argument to
# AC_ARG_ENABLE, we also specify the default behavior (when it isn't given
# at all), here we want the default to be `yes'.
AC_ARG_ENABLE([guide-message],
              [AS_HELP_STRING([--disable-guide-message],
                   [No messages after each build step.])],
	      [], [enable_guide_message=yes])
AC_SUBST(GUIDEMESSAGE, [$enable_guide_message])




# Build the man page directory. Note that this is the cleanest and most
# portable way of building this directory. We don't want to do it in any of
# the Makefiles.
AC_CONFIG_COMMANDS([man page directory], [$MKDIR_P doc/man])





# Make the Makefiles:
AC_OUTPUT





# Print a message if necessary
AS_IF([test x$enable_guide_message = xyes],
[
  AS_ECHO([])
  AS_ECHO([=======================================================================])
  AS_ECHO([=======================================================================])
  AS_ECHO(["$PACKAGE_NAME (Gnuastro) $PACKAGE_VERSION is successfully"])
  AS_ECHO(["configured for this machine."])
  AS_ECHO([])
  AS_IF([test "x$anywarnings" = xyes],
      [
        AS_ECHO(["Configuration warning(s):"])
        AS_ECHO([])
        AS_IF([test "x$has_libjpeg" = "xno"],
              [AS_ECHO(["  - libjpeg, could not be linked with in your library search path."])
               AS_ECHO(["    If JPEG outputs are desired from ConvertType, it will warn"])
               AS_ECHO(["    you and abort with an error."])
               AS_ECHO([]) ])
        AS_IF([test "x$has_libgit2" = "x0"],
              [AS_ECHO(["  - libgit2, could not be linked with in your library search path."])
               AS_ECHO(["    When present, Git's describe output will be stored in the"])
               AS_ECHO(["    output files if Gnuastro's programs were called within a Git"])
               AS_ECHO(["    version controlled directory to help in reproducibility."])
               AS_ECHO([]) ])
        AS_IF([test "x$has_gnulibtool" = "xno"],
              [AS_ECHO(["  - GNU Libtool, could not be found in your search path."])
               AS_ECHO(["    Gnuastro's BuildProgram uses GNU libtool to link your source code"])
               AS_ECHO(["    with the various libraries (Gnuastro's dependencies). Therefore"])
               AS_ECHO(["    BuildProgram will not be built or installed. Please note that"])
               AS_ECHO(["    not having GNU libtool in your search path will not harm the rest"])
               AS_ECHO(["    of Gnuastro's building and installation. Gnuastro has its own"])
               AS_ECHO(["    internal implementation of GNU Libtool to build its self. This"])
               AS_ECHO(["    warning is only to let you know that BuildProgram will not be"])
               AS_ECHO(["    part of this build. The executable names searched were \`libtool'"])
               AS_ECHO(["    and \`glibtool'."])
               AS_ECHO([])
               AS_IF([test "x$has_libtool" = "xyes"],
                     [AS_ECHO(["    -- A libtool implementation was found, but it isn't GNU."])
                      AS_ECHO([])
                     ])
              ])
        AS_IF([test "x$has_ghostscript" = "xno"],
              [AS_ECHO(["  - GPL GhostScript version 9.10 or later, with the executable"])
               AS_ECHO(["    name \`gs', was not found in your PATH environment variable."])
               AS_ECHO(["    If PDF outputs are desired from ConvertType, it will abort"])
               AS_ECHO(["    with an EPS output which you can convert to PDF by other means."])
               AS_ECHO([]) ])
        # The last two scenarios described below are taken from
        # https://unix.stackexchange.com/questions/65700/is-it-safe-to-add-to-my-path-how-come
        AS_IF([test "x$path_warning" = "xyes"],
              [AS_ECHO(["  - Your PATH contains the current directory. This does not affect"])
               AS_ECHO(["    this build and installation of Gnuastro in any way, it is just to"])
               AS_ECHO(["    to remind you that this is a security risk."])
               AS_ECHO(["      It is a very serious security risk if it is closer to the start"])
               AS_ECHO(["    of your PATH: a malicious/wrong program might be run instead of"])
               AS_ECHO(["    a desired program, someone might find out you frequently mistype"])
               AS_ECHO(["    a command and install a matching one, someone might install a"])
               AS_ECHO(["    fake command with the name of one that is not installed. You can"])
               AS_ECHO(["    always run a program in the current directory by explicity adding"])
               AS_ECHO(["    a \`./' before it's name. Run the following command after"])
               AS_ECHO(["    installing Gnuastro to learn more about PATH:"])
               AS_ECHO(["        $ info gnuastro \"Installation directory\""])
               AS_ECHO([]) ])
        AS_ECHO(["  All checks related to the warning(s) above will be skipped."])
        AS_ECHO([])
      ]
     )
  AS_ECHO(["To build Gnuastro $PACKAGE_VERSION, please run:"])
  AS_ECHO([])
  AS_ECHO(["    make -j8"])
  AS_ECHO([])
  AS_ECHO(["(You can change the 8 to the number of CPU threads available.)"])
  AS_ECHO(["(Configure with '--disable-guide-message' for no messages.)"])
  AS_ECHO(["(Please be patient, some libraries can take a few minutes to compile.)"])
  AS_ECHO([=======================================================================])
  AS_ECHO([=======================================================================])
  AS_ECHO([])
])
