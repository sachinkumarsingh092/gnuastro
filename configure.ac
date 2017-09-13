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





# By default we assume no warnings
anywarnings=no





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





# Generic compiler flags for all subdirectories.
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
                   [Depends on bit length of system.])





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




















# The native types for binary arithmetic operations, see the manual for a
# detailed discussion. The initial list of which types to compile can be
# determined with the `--enable-bin-op-alltypes' option. If they have a
# value of 1, they will be compiled. It is possible to disable a previously
# compiled type with with `--disable-bin-op-TYPENAME', or
# `--enable-bin-op-TYPENAME=no'.
binop_alltypes=0
AC_ARG_ENABLE([bin-op-alltypes],
              [AS_HELP_STRING([--enable-bin-op-alltypes],
                    [Allow native binary operations for all types.])],
	      [AS_IF([test "x$enable_bin_op_alltypes" != xno],
                     [binop_alltypes=1], [binop_alltypes=0])], [])
AS_IF([test "x$binop_alltypes" != x0],
      [binop_uint8=1
       binop_int8=1
       binop_uint16=1
       binop_int16=1
       binop_uint32=1
       binop_int32=1
       binop_uint64=1
       binop_int64=1
       binop_float32=1
       binop_float64=1],
      [binop_uint8=1
       binop_int8=0
       binop_uint16=0
       binop_int16=0
       binop_uint32=0
       binop_int32=0
       binop_uint64=1
       binop_int64=1
       binop_float32=1
       binop_float64=1] )

AC_MSG_CHECKING(compilation of 8-bit unsigned int binary operators)
AC_ARG_ENABLE([bin-op-uint8],
              [AS_HELP_STRING([--enable-bin-op-uint8],
                    [Native binary operators on unsigned 8-bit int.])],
	      [AS_IF([test "x$enable_bin_op_uint8" != xno],
                     [binop_uint8=1], [binop_uint8=0])], [])
AS_IF([test "x$binop_uint8" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_UINT8], [$binop_uint8],
                   [Native binary operations on unsigned 8-bit int.])
AC_SUBST(HAVE_BIN_OP_UINT8, [$binop_uint8])
AC_MSG_RESULT($binoptprint)

AC_MSG_CHECKING(compilation of 8-bit signed int binary operators)
AC_ARG_ENABLE([bin-op-int8],
              [AS_HELP_STRING([--enable-bin-op-int8],
                    [Native binary operations on int8 data.])],
	      [AS_IF([test "x$enable_bin_op_int8" != xno],
                     [binop_int8=1], [binop_int8=0])], [])
AS_IF([test "x$binop_int8" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_INT8], [$binop_int8],
                   [Native binary operations on int8 data.])
AC_SUBST(HAVE_BIN_OP_INT8, [$binop_int8])
AC_MSG_RESULT($binoptprint)

AC_MSG_CHECKING(compilation of 16-bit unsigned int binary operators)
AC_ARG_ENABLE([bin-op-uint16],
              [AS_HELP_STRING([--enable-bin-op-uint16],
                    [Native binary operators on unsigned short data.])],
	      [AS_IF([test "x$enable_bin_op_uint16" != xno],
                     [binop_uint16=1], [binop_uint16=0])], [])
AS_IF([test "x$binop_uint16" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_UINT16], [$binop_uint16],
                   [Native binary operations on unsigned short data.])
AC_SUBST(HAVE_BIN_OP_UINT16, [$binop_uint16])
AC_MSG_RESULT($binoptprint)

AC_MSG_CHECKING(compilation of 16-bit signed int binary operators)
AC_ARG_ENABLE([bin-op-int16],
              [AS_HELP_STRING([--enable-bin-op-int16],
                    [Native binary operations on int16 data.])],
	      [AS_IF([test "x$enable_bin_op_int16" != xno],
                     [binop_int16=1], [binop_int16=0])], [])
AS_IF([test "x$binop_int16" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_INT16], [$binop_int16],
                   [Native binary operations on int16 data.])
AC_SUBST(HAVE_BIN_OP_INT16, [$binop_int16])
AC_MSG_RESULT($binoptprint)

AC_MSG_CHECKING(compilation of 32-bit unsigned int binary operators)
AC_ARG_ENABLE([bin-op-uint32],
              [AS_HELP_STRING([--enable-bin-op-uint32],
                    [Native binary operators on unsigned int32 data.])],
	      [AS_IF([test "x$enable_bin_op_uint32" != xno],
                     [binop_uint32=1], [binop_uint32=0])], [])
AS_IF([test "x$binop_uint32" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_UINT32], [$binop_uint32],
                   [Native binary operations on unsigned int32 data.])
AC_SUBST(HAVE_BIN_OP_UINT32, [$binop_uint32])
AC_MSG_RESULT($binoptprint)

AC_MSG_CHECKING(compilation of 32-bit signed int binary operators)
AC_ARG_ENABLE([bin-op-int32],
              [AS_HELP_STRING([--enable-bin-op-int32],
                    [Native binary operations on int32 data.])],
	      [AS_IF([test "x$enable_bin_op_int32" != xno],
                     [binop_int32=1], [binop_int32=0])], [])
AS_IF([test "x$binop_int32" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_INT32], [$binop_int32],
                   [Native binary operations on int32 data.])
AC_SUBST(HAVE_BIN_OP_INT32, [$binop_int32])
AC_MSG_RESULT($binoptprint)

AC_MSG_CHECKING(compilation of 64-bit unsigned int binary operators)
AC_ARG_ENABLE([bin-op-uint64],
              [AS_HELP_STRING([--enable-bin-op-uint64],
                    [Native binary operators on unsigned long data.])],
	      [AS_IF([test "x$enable_bin_op_uint64" != xno],
                     [binop_uint64=1], [binop_uint64=0])], [])
AS_IF([test "x$binop_uint64" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_UINT64], [$binop_uint64],
                   [Native binary operations on unsigned long data.])
AC_SUBST(HAVE_BIN_OP_UINT64, [$binop_uint64])
AC_MSG_RESULT($binoptprint)

AC_MSG_CHECKING(compilation of 64-bit signed int binary operators)
AC_ARG_ENABLE([bin-op-int64],
              [AS_HELP_STRING([--enable-bin-op-int64],
                    [Native binary operations on int64 data.])],
	      [AS_IF([test "x$enable_bin_op_int64" != xno],
                     [binop_int64=1], [binop_int64=0])], [])
AS_IF([test "x$binop_int64" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_INT64], [$binop_int64],
                   [Native binary operations on long64 data.])
AC_SUBST(HAVE_BIN_OP_INT64, [$binop_int64])
AC_MSG_RESULT($binoptprint)

AC_MSG_CHECKING(compilation of 32-bit floating point binary operators)
AC_ARG_ENABLE([bin-op-float32],
              [AS_HELP_STRING([--enable-bin-op-float32],
                    [Native binary operations on float32 data.])],
	      [AS_IF([test "x$enable_bin_op_float32" != xno],
                     [binop_float32=1], [binop_float32=0])], [])
AS_IF([test "x$binop_float32" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_FLOAT32], [$binop_float32],
                   [Native binary operations on float32 data.])
AC_SUBST(HAVE_BIN_OP_FLOAT32, [$binop_float32])
AC_MSG_RESULT($binoptprint)

AC_MSG_CHECKING(compilation of 64-bit floating point binary operators)
AC_ARG_ENABLE([bin-op-float64],
              [AS_HELP_STRING([--enable-bin-op-float64],
                    [Native binary operations on float64 data.])],
	      [AS_IF([test "x$enable_bin_op_float64" != xno],
                     [binop_float64=1], [binop_float64=0])], [])
AS_IF([test "x$binop_float64" != x0], [binoptprint=yes], [binoptprint=no])
AC_DEFINE_UNQUOTED([GAL_CONFIG_BIN_OP_FLOAT64], [$binop_float64],
                   [Native binary operations on float64 data.])
AC_SUBST(HAVE_BIN_OP_FLOAT64, [$binop_float64])
AC_MSG_RESULT($binoptprint)




















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
	      [AS_IF([test "x$enable_Fits" != xno],
                     [enable_fits=yes; ayes=true])],
              [enable_fits=notset])
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
      [
       AS_IF([test $enable_arithmetic = notset], [enable_arithmetic=no])
       AS_IF([test $enable_buildprog = notset], [enable_buildprog=no])
       AS_IF([test $enable_convertt = notset], [enable_convertt=no])
       AS_IF([test $enable_convolve = notset], [enable_convolve=no])
       AS_IF([test $enable_cosmiccal = notset], [enable_cosmiccal=no])
       AS_IF([test $enable_crop = notset], [enable_crop=no])
       AS_IF([test $enable_fits = notset], [enable_fits=no])
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





# Make the enable_package values available for the Makefile:
AM_CONDITIONAL([COND_ARITHMETIC], [test $enable_arithmetic = yes])
AM_CONDITIONAL([COND_BUILDPROG], [test $enable_buildprog = yes])
AM_CONDITIONAL([COND_CONVERTT], [test $enable_convertt = yes])
AM_CONDITIONAL([COND_CONVOLVE], [test $enable_convolve = yes])
AM_CONDITIONAL([COND_COSMICCAL], [test $enable_cosmiccal = yes])
AM_CONDITIONAL([COND_CROP], [test $enable_crop = yes])
AM_CONDITIONAL([COND_FITS], [test $enable_fits = yes])
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
        AS_IF([test "x$has_ghostscript" = "xno"],
              [AS_ECHO(["  - GPL GhostScript version 9.10 or later, with the executable"])
               AS_ECHO(["    name \`gs', was not found in your PATH environment variable."])
               AS_ECHO(["    If PDF outputs are desired from ConvertType, it will abort"])
               AS_ECHO(["    with an EPS output which you can convert to PDF by other means."])
               AS_ECHO([]) ])
        AS_ECHO(["  All checks related to the warning(s) will be skipped."])
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
