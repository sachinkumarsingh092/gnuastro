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
#     Mohammad Akhlaghi <mohammad@akhlaghi.org>
# Contributing author(s):
# Copyright (C) 2015-2020, Free Software Foundation, Inc.
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
GAL_CURRENT=10
GAL_REVISION=0
GAL_AGE=0
GAL_LT_VERSION="${GAL_CURRENT}:${GAL_REVISION}:${GAL_AGE}"
AC_SUBST(GAL_LT_VERSION)





# Checks for programs.
: ${CFLAGS=""}
: ${CXXFLAGS=""}
AC_PROG_CC
AC_PROG_CXX
gl_EARLY
AM_PROG_AR
LT_INIT





# This macro will let the libraries know that we are now in the Gnuastro
# build system, not on the user's system. While we are building Gnuastro,
# we have the important installation information in 'config.h'. But in the
# user's own programs, this information is defined in
# 'gnuastro/config.h'. With this macro, the installed headers can decide
# if the latter should be included or not. Note that 'gnuastro/config.h'
# is only built at installation time and doesn't exist when building
# Gnuastro. Therefore, this macro must not be defined in a user's program.
AC_DEFINE([IN_GNUASTRO_BUILD], [1], [In building, not usage])





# See if 'make check' should be made with Valgrind. This should be done
# before checking the compilation flags because it can re-write
# 'enable-debug').
AC_ARG_ENABLE(check-with-valgrind,
              [AS_HELP_STRING([--enable-check-with-valgrind],
                              [Run 'make check' programs within Valgrind])],
              [AS_IF([test "x$enable_check_with_valgrind" != xno],
                     [enable_check_with_valgrind=yes])],
              [enable_check_with_valgrind=no])
AS_IF([test "x$enable_check_with_valgrind" = "xyes"],
      [
        enable_debug=yes;
        AC_CHECK_PROG(has_valgrind, valgrind, [yes], [no])
        AS_IF([test "x$has_valgrind" = "xno"],
              [AC_MSG_ERROR([Valgrind not found. Please install it or don't use '--enable-check-with-valgrind'])])
      ])
AM_CONDITIONAL([COND_CHECK_WITH_VALGRIND], [test "x$enable_check_with_valgrind" = "xyes"])





# Set the compilation flags.
AC_ARG_ENABLE(debug,
              [AS_HELP_STRING([--enable-debug],
                              [No optimization, debug flags, no shared lib.])],
              [AS_IF([test "x$enable_debug" != xno], [enable_debug=yes])],
              [enable_debug=no])
AS_IF([test "x$enable_debug" = "xyes"],
      [cflags_add="-g -O0"; enable_shared=no],
      [cflags_add="-O3"])
CFLAGS="-Wall $cflags_add $CFLAGS"
CXXFLAGS="-Wall $cflags_add $CXXFLAGS"





# See if the C++ compiler was present. 'CXX' has already been set by
# 'AC_PROG_CXX' (above). According to the Autoconf manual: "if none of the
# [AC_PROG_CXX] checks succeed, then as a last resort [it will] set 'CXX'
# to 'g++'". Therefore, we can't rely on it to see if the compiler
# executable is actually usable. Therefore tests that rely on it will
# fail. Unfortunately some OSs (like Fedora), don't install 'g++' with
# 'gcc', so it may not always be present. To fix this, here we are just
# using 'AC_CHECK_PROG' to see if the 'CXX' executable is reachable in the
# search path and if it isn't, we'll disable the C++ check during 'make
# check' with an Automake conditional.
AC_CHECK_PROG(has_cxx, $CXX, "yes", "no")
AM_CONDITIONAL([COND_HASCXX], [test "x$has_cxx" = "xyes"])





# Check for pthreads and add the appropriate compilation flags. AX_PTHREAD
# comes from the GNU Autoconf Archive's ax_pthread.m4, see there for the
# documentation. Note that
AX_PTHREAD([],[AC_MSG_ERROR([AC_PACKAGE_NAME Needs POSIX Threads (pthread)])])
CXXFLAGS="$CXXFLAGS $PTHREAD_CFLAGS"
CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
LIBS="$PTHREAD_LIBS $LIBS"
CC="$PTHREAD_CC"





# See if the C++ compiler understands '-Qunused-arguments'. AX_PTHREAD adds
# puts this option in 'PTHREAD_CFLAGS' when the C compiler knows this
# option. We then pass it to CFLAGS and CXXFLAGS above. But as reported in
# bug #52490, it can happen that sometimes, the C++ compiler doesn't
# recognize it. So we need to do a separate check for C++.
cxxflags_tmp=
for flg in $CXXFLAGS; do
  AS_IF([test "$flg" = \-Qunused-arguments],
        [ AC_LANG(C++)
          AX_CHECK_COMPILE_FLAG([-Qunused-arguments],
                                [cxx_Qunused_arguments=yes;
                                 cxxflags_tmp="$cxxflags_tmp $flg"],
                                [cxx_Qunused_arguments=no])
          AC_LANG(C) ],
        [ cxxflags_tmp="$cxxflags_tmp $flg"])
done
CXXFLAGS="$cxxflags_tmp"





# Check if 'malloc(0)' returns valid pointer
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
AC_CHECK_SIZEOF([int])
AC_SUBST(SIZEOF_INT, [$ac_cv_sizeof_int])
AC_DEFINE_UNQUOTED([GAL_CONFIG_SIZEOF_INT], [$ac_cv_sizeof_int],
                   [It may be 16bit or 32bit])





# By default we assume no warnings
anywarnings=no





# Remove any occurance of the current directory './', '.', or the full
# address of the current directory in PATH. The main problem is the
# 'libtool' executable which Gnuastro builds internally in the top build
# directory. However, we also need to know if the system has libtool or
# not.
AC_MSG_CHECKING(if PATH contains current directory)
oldPATH=$PATH
currpwd=$(pwd)

# The first call to SED will remove any occurance of the current directory:
# './', '.', or the full address.
#
#    NOTE 1: We cannot simply remove all '.'s, because hidden directories
#            (like the '~/.local' that is suggested for local
#            installations) will also be altered.
#
#    NOTE 2: An empty string in the list of strings (separated by ':')
#            means the current directory. This includes cases like: '::',
#            or a leading and trailing ':'. So after all the removals of
#            the current directory, we will remove all such cases.
#
#    NOTE 3: The SED separator can be any character immediately after 's',
#            it doesn't just have to be the commonly used '/'. Since '$pwd'
#            will possibly contain many '/'s, it is much more easier to use
#            a differen separator ('|' in this call to SED).
PATH=$(AS_ECHO([$PATH]) | $SED -e 's|'"$currpwd"'||g' \
                               -e 's|\.\.*//*||g'     \
                               -e 's|:\.\.*:|:|g'     \
                               -e 's|\.*$||'          \
                               -e 's|^\.*||'          \
                               -e 's|::*|:|g'         \
                               -e 's|^:||'            \
                               -e 's|:$||' )
AS_IF([test $oldPATH = $PATH],
      [ path_warning=no                   ],
      [ path_warning=yes; anywarnings=yes ])
AC_MSG_RESULT( $path_warning )





# Libraries
# ---------
#
# After each library is found, AC_SEARCH_LIBS adds the -lLIBRARY flag to
# the LIBS variable which is then given to all the Makefiles. Each new flag
# is added to the left of the old one so order matters here. Note that the
# LIBS variable is also used in checking the next libraries, so the linking
# with their dependent libraries is done automatically with this order, and
# we don't have to explicitly set the dependency flags.
has_gsl=yes
has_libgit2=1
has_cmath=yes
has_wcslib=yes
has_cfitsio=yes
has_libtiff=yes
has_libjpeg=yes
has_gslcblas=yes
missing_mandatory=no
missing_optional_lib=no

# Keep the original LIBS to re-set in the end.
orig_LIBS="$LIBS"

# Order is important here.
AC_LIB_HAVE_LINKFLAGS([m], [], [#include <math.h>])
AS_IF([test "x$LIBM" = x],
      [missing_mandatory=yes; has_cmath=no],
      [LDADD="$LIBM $LDADD"; LIBS="$LIBM $LIBS"])

AC_LIB_HAVE_LINKFLAGS([gsl], [gslcblas], [
#include <gsl/gsl_rng.h>
void junk(void) { gsl_rng_env_setup(); } ])
AS_IF([test "x$LIBGSL" = x],
      [missing_mandatory=yes; has_gsl=no; has_gslcblas=no],
      [LDADD="$LTLIBGSL $LDADD"; LIBS="$LIBGSL $LIBS"])


# Since version 0.42, if 'libcurl' is installed, CFITSIO will link with it
# and thus it will be necessary to explicitly link with libcurl also. If it
# doesn't exist on the system, then CFITSIO won't link with it and there is
# no problem for Gnuastro either. So there is no need to stop the configure
# script if libcurl isn't found.
#
# For static builds, linking libcurl (also built in static-only mode) will
# depend on linking with zlib. So we'll need to also search for that
# also. But if libcurl isn't built in static-only mode, then it will have
# undefined symbols for many libraries (for example libpsl, libidn2,
# librtmp, libldap). So if you intend to make Gnuastro statically, then
# build Libcurl in static-only mode so you won't have to check for all
# these extra libraries here.
AC_LIB_HAVE_LINKFLAGS([z], [], [#include <zlib.h>])
AS_IF([test "x$LIBZ" = x], [],
      [LDADD="$LTLIBZ $LDADD"; LIBS="$LIBZ $LIBS"])


AC_LIB_HAVE_LINKFLAGS([curl], [], [#include <curl/curl.h>])
AS_IF([test "x$LIBCURL" = x], [],
      [LDADD="$LTLIBCURL $LDADD"; LIBS="$LIBCURL $LIBS"])


# Older versions of CFITSIO don't install a shared library, only a static
# one. Eventhough we add '-lm' to LDADD, on some systems, it complains
# about not finding basic math libraries. Therefore the configure script
# can't find CFITSIO, eventhough it exists. The solution is to manually add
# the math-library as a dependency of CFITSIO.
AC_LIB_HAVE_LINKFLAGS([cfitsio], [], [
#include <fitsio.h>
void junk(void) {
int status;
fitsfile *f;
ffopen(&f, "junk", READONLY, &status);} ])
AS_IF([test "x$LIBCFITSIO" = x],
      [missing_mandatory=yes; has_cfitsio=no],
      [LDADD="$LTLIBCFITSIO $LDADD"; LIBS="$LIBCFITSIO $LIBS"])


AC_LIB_HAVE_LINKFLAGS([wcs], [], [
#include <wcslib/wcshdr.h>
void junk(void) {
int nreject, nwcs;
struct wcsprm *wcs;
char *header="JUNK";
wcspih(header, 1, 0, 0, &nreject, &nwcs, &wcs);
} ])
AS_IF([test "x$LIBWCS" = x],
      [missing_mandatory=yes; has_wcslib=no],
      [LDADD="$LTLIBWCS $LDADD"; LIBS="$LIBWCS $LIBS"])


AC_ARG_WITH([libjpeg],
            [AS_HELP_STRING([--without-libjpeg],
                            [disable support for libjpeg])],
            [], [with_libjpeg=yes])
AS_IF([test "x$with_libjpeg" != xno],
      [ AC_LIB_HAVE_LINKFLAGS([jpeg], [], [
#include <stdio.h>
#include <stdlib.h>
#include <jpeglib.h>
void junk(void) {
  struct jpeg_decompress_struct cinfo;
  jpeg_create_decompress(&cinfo);
} ]) ])
AS_IF([test "x$LIBJPEG" = x],
      [missing_optional_lib=yes; has_libjpeg=no; anywarnings=yes],
      [LDADD="$LTLIBJPEG $LDADD"; LIBS="$LIBJPEG $LIBS"])
AM_CONDITIONAL([COND_HASLIBJPEG], [test "x$has_libjpeg" = "xyes"])


# Check libtiff: If an LZMA libray (part of the XZ Utils) is present,
# libtiff has probably been built with it. So we'll also need to link with
# the LZMA library. But if libtiff hasn't been linked with it and its
# present, there is no problem, the linker will just pass over it. So we
# don't need to stop the build if this fails.
AC_ARG_WITH([libtiff],
            [AS_HELP_STRING([--without-libtiff],
                            [disable support for libtiff])],
            [], [with_libtiff=yes])
AS_IF([test "x$with_libtiff" != xno],
      [ AC_LIB_HAVE_LINKFLAGS([lzma], [], [#include <lzma.h>])
        AC_LIB_HAVE_LINKFLAGS([tiff], [], [
#include <tiffio.h>
void junk(void) {TIFF *tif=TIFFOpen("junk", "r");} ])
      ])
AS_IF([test "x$LIBLZMA" = x], [],
      [LDADD="$LTLIBLZMA $LDADD"; LIBS="$LIBLZMA $LIBS"])
AS_IF([test "x$LIBTIFF" = x],
      [missing_optional_lib=yes; has_libtiff=no; anywarnings=yes],
      [LDADD="$LTLIBTIFF $LDADD"; LIBS="$LIBTIFF $LIBS"])
AM_CONDITIONAL([COND_HASLIBTIFF], [test "x$has_libtiff" = "xyes"])


# Check libgit2. Note that very old versions of libgit2 don't have the
# 'git_libgit2_init' function.
AC_ARG_WITH([libgit2],
            [AS_HELP_STRING([--without-libgit2],
                            [disable support for libgit2])],
            [], [with_libgit2=yes])
AS_IF([test "x$with_libgit2" != xno],
      [ AC_LIB_HAVE_LINKFLAGS([git2], [], [
#include <git2.h>
void junk(void) {git_libgit2_init();} ])
      ])
AS_IF([test "x$LIBGIT2" = x],
      [missing_optional_lib=yes; has_libgit2=0],
      [LDADD="$LTLIBGIT2 $LDADD"; LIBS="$LIBGIT2 $LIBS"])
AC_DEFINE_UNQUOTED([GAL_CONFIG_HAVE_LIBGIT2], [$has_libgit2],
                   [libgit2 is installed on the system])
AS_IF([test "x$has_libgit2" = "x1"], [], [anywarnings=yes])





# Higher-level library tests
# --------------------------
#
# Once we know that a library exsits, we need to check if it has some
# features or not. This must be done _after_ checking the existance of
# _all_ the libraries, because they may add elements to 'LIBS'/'LDADD' that
# causes possibly different versions of the libraries to be read.

# GSL's 'gsl_interp_steffen' isn't a function. So we'll need to use
# 'AC_LINK_IFELSE'. However, AC_LINK_IFELSE doesn't use 'LDADD', so we'll
# have to temporarily add 'LDADD' to LIBS, then set it back to the
# original.
AC_MSG_CHECKING(if GSL supports Steffen splines)
AC_LINK_IFELSE([AC_LANG_PROGRAM(
                   [[#include <gsl/gsl_interp.h>]],
                   [[const gsl_interp_type *itype=gsl_interp_steffen;]])],
               [AC_MSG_RESULT(yes)
                gsl_version_old=no],
               [AC_MSG_RESULT(no)
                gsl_version_old=yes; anywarnings=yes;])


# If the CFITSIO library has the 'fits_is_reentrant' function (it was added
# since version 3.30 of April 2012).
AC_CHECK_LIB([cfitsio], [fits_is_reentrant], [has_fits_is_reentrant=1],
             [has_fits_is_reentrant=0], [-lm])
AC_DEFINE_UNQUOTED([GAL_CONFIG_HAVE_FITS_IS_REENTRANT],
                   [$has_fits_is_reentrant],
                   [CFITSIO has the fits_is_reentrant function])
AC_SUBST(HAVE_FITS_IS_REENTRANT, [$has_fits_is_reentrant])


# If the WCS library has the 'wcslib_version' function.
AC_CHECK_LIB([wcs], [wcslib_version], [has_wcslib_version=1],
             [has_wcslib_version=0], [-lcfitsio -lm])
AC_DEFINE_UNQUOTED([GAL_CONFIG_HAVE_WCSLIB_VERSION], [$has_wcslib_version],
                   [WCSLIB comes with wcslib_version])
AC_SUBST(HAVE_WCSLIB_VERSION, [$has_wcslib_version])


# If the pthreads library has 'pthread_barrier_wait'.
AC_CHECK_LIB([pthread], [pthread_barrier_wait], [has_pthread_barrier=1],
             [has_pthread_barrier=0])
AC_DEFINE_UNQUOTED([GAL_CONFIG_HAVE_PTHREAD_BARRIER], [$has_pthread_barrier],
                   [System has pthread_barrier])
AC_SUBST(HAVE_PTHREAD_BARRIER, [$has_pthread_barrier])





# Programs
# --------

# Help2man:
AC_CHECK_PROG(has_help2man, help2man, [yes], [no])
AM_CONDITIONAL([COND_HASHELP2MAN], [test "x$has_help2man" = "xyes"])

# Check the libtool executable on the system. Note that Gnuastro also ships
# with a version of Libtool. We don't want Gnuastro's Libtool, here we want
# to see if the system has libtool independent of Gnuastro so BuildProgram
# can use it later. We also need to check some shells to run libtool within
# them.
AC_CHECK_PROG(has_libtool, libtool, [yes], [no])
AC_CHECK_PROG(has_bash,    bash,    [yes], [no])
AC_CHECK_PROG(has_zsh,     zsh,     [yes], [no])

# If Libtool is present, make sure it is GNU Libtool
AS_IF([test "x$has_libtool" = "xyes"],
      [ AC_MSG_CHECKING(if libtool executable is GNU)
        AS_IF( libtool --version 2> /dev/null | grep GNU 2>&1 > /dev/null,
               [has_gnulibtool=yes], [has_gnulibtool=no])
        AC_MSG_RESULT( $has_gnulibtool )],
      [ has_gnulibtool=no ])

# When either the 'libtool' executable isn't GNU or it doesn't exist, then
# look for 'glibtool'.
AS_IF([test "x$has_gnulibtool" = "xyes"],
      [ gnulibtool_exec=libtool ],
      [ AC_CHECK_PROG(has_glibtool, glibtool, [yes], [no])
        AS_IF( [test "x$has_glibtool" = "xyes"],
               [has_gnulibtool=yes; gnulibtool_exec=glibtool],
               [has_gnulibtool=no; anywarnings=yes] ) ])

# Older versions of GNU Libtool have problems with the 'dash' shell (a
# minimalist shell) and will crash (due to not having the '+=' operator),
# see bug #54430. So we need to check if the system's libtool and shell
# (called by 'sh' as in C's 'system' function within BuildProgram) can work
# with each other. If not, we need to search for Bash and tell BuildProgram
# to use Bash instead of the default. But older versions of Bash also don't
# support this operator, so we'll have to check with Bash is well.
libtool_shell="none"
AS_IF([test "x$has_gnulibtool" = "xyes"],
      [
        AC_MSG_CHECKING(for shell to use with libtool)

        # Make a C source file and run Libtool on it with the specified
        # shells.
        outname=libtool_shell_test
        cprog=libtool_shell_test.c
        AS_ECHO(["#include <stdio.h>"])                                 > $cprog
        AS_ECHO(["int main(void){printf(\"success\\n\"); return 0;}"]) >> $cprog
        ltargs="--quiet --tag=CC --mode=link $CC $cprog -O3 -o $outname"

        # Check the shells, starting with known shells and ultimately
        # trying with 'sh' (can be any shell).
        AS_IF([test "x$has_bash" = "xyes"],
              [AS_IF(bash -c "$gnulibtool_exec $ltargs" > /dev/null 2>&1,
                     [libtool_shell="bash"],
                     [bash_version=$BASH_VERSION]) ],
              [AS_IF([test "x$has_zsh" = "xyes"],
                     [AS_IF(zsh -c "$gnulibtool_exec $ltargs" > /dev/null 2>&1,
                            [libtool_shell="zsh"]) ],
                     [AS_IF(sh -c "$gnulibtool_exec $ltargs" > /dev/null 2>&1,
                            [libtool_shell="sh"]) ]) ])

        # Clean up: note that no output might have been generated (when no
        # proper shell was found). Therefore, for deleting the output file,
        # we'll call 'rm' with '-f' so it doesn't complain with an error in
        # such cases.
        rm $cprog
        rm -f $outname
        AC_MSG_RESULT($libtool_shell)
      ])

# If a good shell to call Libtool could be found, then Libtool is usable
# within a program (BuildProgram in this case).
AS_IF([test "x$libtool_shell" = "xnone"],
      [
        anywarnings=yes;
        usable_libtool=no;
        AS_IF([test "x$bash_version" = x], [junk=1],
              [AS_ECHO(["GNU Bash was found but not a suitable version: $bash_version"]) ])
      ],
      [
        usable_libtool=yes
        AC_DEFINE_UNQUOTED([GAL_CONFIG_GNULIBTOOL_SHELL], ["$libtool_shell"],
                           [Shell program to use with GNU Libtool])
        AC_DEFINE_UNQUOTED([GAL_CONFIG_GNULIBTOOL_EXEC], ["$gnulibtool_exec"],
                           [The executable to call GNU Libtool])
      ])





# Check Ghostscript: "-dPDFFitPage" option to Ghostscript, used by the
# library to convert from EPS to PDF, has been introduced in Ghostscript
# 9.10.  Make sure we have at least that version.
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
# Note: 'has_ghostscript' can be set to 'no' within the AS_IF above, so
# 'anywarnings' cannot be an [RUN-IF-FALSE] argument to the AS_IF above.
AS_IF([test "x$has_ghostscript" = "xno"], [anywarnings=yes])
AM_CONDITIONAL([COND_HASGHOSTSCRIPT], [test "x$has_ghostscript" = "xyes"])





# If any necessary dependency is missing inform the user and abort.
AS_IF([test "x$missing_mandatory" = "xyes"],
      [
        # List missing packages: print the GSL CBLAS message only if GSL is
        # present. Otherwise, it is just confusing for the users (CBLAS
        # will be installed with GSL). The CBLAS message is only
        # interesting if the GSL test has passed.
        AS_ECHO([""])
        AS_ECHO(["Missing MANDATORY dependencies (necessary to continue):"])
        AS_IF([test "x$has_cmath" = "xno"],
              [ AS_ECHO([" - C library (math): This may be the cause of all other failures."]) ])
        AS_IF([test "x$has_gsl" = "xno"],
              [ AS_ECHO([" - GNU Scientific Library (GSL): https://www.gnu.org/software/gsl"]) ],
              [ AS_IF([test "x$has_gslcblas" = "xno"],
                      [ AS_ECHO([" - The BLAS support of GNU Scientific Library (GSL). This should have"])
                        AS_ECHO(["   been installed along with GSL. Try re-installing GSL."]) ]) ])
        AS_IF([test "x$has_cfitsio" = "xno"],
              [ AS_ECHO([" - CFITSIO: https://heasarc.gsfc.nasa.gov/fitsio"]) ])
        AS_IF([test "x$has_wcslib" = "xno"],
              [ AS_ECHO([" - WCSLIB: http://www.atnf.csiro.au/people/mcalabre/WCS"]) ])

        # Optional dependencies:
        AS_IF([test "x$anywarnings" = xyes],
              [
                AS_ECHO([""])
                AS_ECHO(["OPTIONAL warnings/dependencies (for improved functionality):"])
                AS_IF([test "x$gsl_version_old" = "xyes"],
                      [ AS_ECHO([" - Old GSL version: https://www.gnu.org/s/gsl"])                ])
                AS_IF([test "x$has_libjpeg" = "xno"],
                      [ AS_ECHO([" - Missing Libjpeg (JPEG files): http://ijg.org"])              ])
                AS_IF([test "x$has_libtiff" = "xno"],
                      [ AS_ECHO([" - Missing Libtiff (TIFF files): http://libtiff.maptools.org"]) ])
                AS_IF([test "x$has_libgit2" = "x0"],
                      [ AS_ECHO([" - Missing Libgit2: https://libgit2.org"])                      ])
                AS_IF([test "x$usable_libtool" = "xno"],
                      [ AS_ECHO([" - Unusable GNU Libtool: https://www.gnu.org/s/libtool"])
                        AS_IF([test "x$has_gnulibtool" = "xyes"],
                              [ AS_ECHO(["  -- GNU Libtool is present, tested shells not supported."])  ])
                        AS_IF([test "x$has_libtool" = "xyes"],
                              [ AS_ECHO(["  -- A libtool implementation was found, but it isn't GNU."]) ])
                      ])
                AS_IF([test "x$has_ghostscript" = "xno"],
                      [ AS_ECHO([" - Missing GPL Ghostscript (v9.10 or later): https://www.ghostscript.com"])])
              ])

        # Suggestions on fixing the problem.
        AS_ECHO([""])
        AS_ECHO(["You can use your package manager for easy and fast installation of all"])
        AS_ECHO(["dependencies in one command. See the link below:"])
        AS_ECHO(["  https://www.gnu.org/s/gnuastro/manual/html_node/Dependencies-from-package-managers.html"])
        AS_ECHO([""])
        AS_ECHO(["If you have already installed a dependency (for example in '/install/path'),"])
        AS_ECHO(["but this script can't link with it, add the path to the LDFLAGS, CPPFLAGS and"])
        AS_ECHO(["LD_LIBRARY_PATH environment variables before running configure. For example"])
        AS_ECHO(["with the following commands (just correct the '/install/path' part)."])
        AS_ECHO(["  $ export LDFLAGS=\"\$LDFLAGS -L/install/path/lib\""])
        AS_ECHO(["  $ export CPPFLAGS=\"\$CPPFLAGS -I/install/path/include\""])
        AS_ECHO(["  $ export LD_LIBRARY_PATH=\"\$LD_LIBRARY_PATH:/install/path/lib\""])
        AS_ECHO([""])
        AS_ECHO(["[TIP:] Put these commands in your startup file (for example '~/.bashrc') to"])
        AS_ECHO(["avoid similar problems later. See the link below to learn more:"])
        AS_ECHO(["  https://www.gnu.org/s/gnuastro/manual/html_node/Installation-directory.html"])
        AS_ECHO([""])
        AS_ECHO([""])
        AC_MSG_ERROR([Mandatory dependency(s) missing, === SEE MESSAGE ABOVE ===.])
      ])





# Gnulib checks:
gl_INIT





# Check if Gnulib tests should be done:
AC_ARG_ENABLE([gnulibcheck],
              [AS_HELP_STRING([--enable-gnulibcheck],
  	            [In 'make check', also test GNU Gnulib.])],
	      [enable_gnulibcheck=yes], [enable_gnulibcheck=no])
AM_CONDITIONAL([COND_GNULIBCHECK], [test $enable_gnulibcheck = yes])





# Gnulib checks for the proper name for the C99 equivalent 'restrict'
# keyword and puts it in the 'ac_cv_c_restrict' variable. If none exists,
# it will put a 'no' inside of this variable. As described in the output
# 'bootstrapped/m4/gnulib-common.m4', this is only necessary until Autoconf
# 2.70 is released. Afterwards, we can use AC_C_RESTRICT.
AS_IF([test "x$ac_cv_c_restrict" = "xno"],
      [gal_restrict_replace=], [gal_restrict_replace=$ac_cv_c_restrict])
AC_SUBST(RESTRICT_REPLACEMENT, [$gal_restrict_replace])





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
AC_ARG_ENABLE([segment],
              [AS_HELP_STRING([--enable-segment],
                    [Install Segment and other enabled programs.])],
	      [AS_IF([test "x$enable_segment" != xno],
                     [enable_segment=yes; ayes=true])],
              [enable_segment=notset])
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





# If we had a "ayes" variable to be "true" if there was a 'yes'. So any
# program that is not explicitly requested must be ignored and vice versa
# (if no programs were explicitly requested, then enable all that weren't
# disabled).
AS_IF([test $ayes = true ],
      [
       AS_IF([test $enable_arithmetic = notset],  [enable_arithmetic=no])
       AS_IF([test $enable_buildprog = notset],   [enable_buildprog=no])
       AS_IF([test $enable_convertt = notset],    [enable_convertt=no])
       AS_IF([test $enable_convolve = notset],    [enable_convolve=no])
       AS_IF([test $enable_cosmiccal = notset],   [enable_cosmiccal=no])
       AS_IF([test $enable_crop = notset],        [enable_crop=no])
       AS_IF([test $enable_fits = notset],        [enable_fits=no])
       AS_IF([test $enable_match = notset],       [enable_match=no])
       AS_IF([test $enable_mkcatalog = notset],   [enable_mkcatalog=no])
       AS_IF([test $enable_mknoise = notset],     [enable_mknoise=no])
       AS_IF([test $enable_mkprof = notset],      [enable_mkprof=no])
       AS_IF([test $enable_noisechisel = notset], [enable_noisechisel=no])
       AS_IF([test $enable_segment = notset],     [enable_segment=no])
       AS_IF([test $enable_statistics = notset],  [enable_statistics=no])
       AS_IF([test $enable_table = notset],       [enable_table=no])
#      AS_IF([test $enable_TEMPLATE = notset],    [enable_TEMPLATE=no])
       AS_IF([test $enable_warp = notset],        [enable_warp=no])
       ],

      [
       AS_IF([test $enable_arithmetic = notset],  [enable_arithmetic=yes])
       AS_IF([test $enable_buildprog = notset],   [enable_buildprog=yes])
       AS_IF([test $enable_convertt = notset],    [enable_convertt=yes])
       AS_IF([test $enable_convolve = notset],    [enable_convolve=yes])
       AS_IF([test $enable_cosmiccal = notset],   [enable_cosmiccal=yes])
       AS_IF([test $enable_crop = notset],        [enable_crop=yes])
       AS_IF([test $enable_fits = notset],        [enable_fits=yes])
       AS_IF([test $enable_match = notset],       [enable_match=yes])
       AS_IF([test $enable_mkcatalog = notset],   [enable_mkcatalog=yes])
       AS_IF([test $enable_mknoise = notset],     [enable_mknoise=yes])
       AS_IF([test $enable_mkprof = notset],      [enable_mkprof=yes])
       AS_IF([test $enable_noisechisel = notset], [enable_noisechisel=yes])
       AS_IF([test $enable_segment = notset],     [enable_segment=yes])
       AS_IF([test $enable_statistics = notset],  [enable_statistics=yes])
       AS_IF([test $enable_table = notset],       [enable_table=yes])
#      AS_IF([test $enable_TEMPLATE = notset],    [enable_TEMPLATE=yes])
       AS_IF([test $enable_warp = notset],        [enable_warp=yes])
       ]
     )





# BuildProgram depends on the presence of GNU Libtool, if it isn't present,
# then don't build it.
AS_IF([test "x$usable_libtool" = "xno"], [enable_buildprog=no])





# Make the enable_package values available for the Makefile:
AM_CONDITIONAL([COND_ARITHMETIC],  [test $enable_arithmetic = yes])
AM_CONDITIONAL([COND_BUILDPROG],   [test $enable_buildprog = yes])
AM_CONDITIONAL([COND_CONVERTT],    [test $enable_convertt = yes])
AM_CONDITIONAL([COND_CONVOLVE],    [test $enable_convolve = yes])
AM_CONDITIONAL([COND_COSMICCAL],   [test $enable_cosmiccal = yes])
AM_CONDITIONAL([COND_CROP],        [test $enable_crop = yes])
AM_CONDITIONAL([COND_FITS],        [test $enable_fits = yes])
AM_CONDITIONAL([COND_MATCH],       [test $enable_match = yes])
AM_CONDITIONAL([COND_MKCATALOG],   [test $enable_mkcatalog = yes])
AM_CONDITIONAL([COND_MKNOISE],     [test $enable_mknoise = yes])
AM_CONDITIONAL([COND_MKPROF],      [test $enable_mkprof = yes])
AM_CONDITIONAL([COND_NOISECHISEL], [test $enable_noisechisel = yes])
AM_CONDITIONAL([COND_SEGMENT],     [test $enable_segment = yes])
AM_CONDITIONAL([COND_STATISTICS],  [test $enable_statistics = yes])
AM_CONDITIONAL([COND_TABLE],       [test $enable_table = yes])
#AM_CONDITIONAL([COND_TEMPLATE],   [test $enable_TEMPLATE = yes])
AM_CONDITIONAL([COND_WARP],        [test $enable_warp = yes])





# Reset LIBS to the initial value BEFORE generating the Makefiles (so the
# modified value doesn't get written into them). Then report the final
# linking flags and put them in the Makefiles.
LIBS="$orig_LIBS"
AC_SUBST(CONFIG_LDADD, [$LDADD])
AM_CONDITIONAL([COND_NORPATH], [test "x$enable_rpath" = "xno"])
AS_ECHO(["linking flags (LDADD) ... $LDADD"])





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
                 bin/script/Makefile
                 bin/mknoise/Makefile
                 bin/segment/Makefile
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
# '--disable-guide-message' is equivalent to a value of 'no', while with no
# argument, the value will default to 'yes'. In the last argument to
# AC_ARG_ENABLE, we also specify the default behavior (when it isn't given
# at all), here we want the default to be 'yes'.
AC_ARG_ENABLE([guide-message],
              [AS_HELP_STRING([--disable-guide-message],
                   [No messages after each build step.])],
	      [], [enable_guide_message=yes])
AC_SUBST(GUIDEMESSAGE, [$enable_guide_message])




# Build the man page directory. Note that this is the cleanest and most
# portable way of building this directory. We don't want to do it in any of
# the Makefiles.
AC_CONFIG_COMMANDS([man page directory], [$MKDIR_P doc/man])





# Prepare the Makefiles.
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

        AS_IF([test "x$gsl_version_old" = "xyes"],
              [dependency_notice=yes
               AS_ECHO(["  - GNU Scientific Library (GSL: https://www.gnu.org/s/gsl) version"])
               AS_ECHO(["    on this system doesn't have some features that can be useful in"])
               AS_ECHO(["    some parts of Gnuastro. This build won't crash, but Gnuastro will"])
               AS_ECHO(["    have less functionality afterwards. We thus recommend building"])
               AS_ECHO(["    and installing a more recent version of GSL (version >= 2.0,"])
               AS_ECHO(["    released in October 2015)."])
               AS_ECHO([]) ])

        AS_IF([test "x$has_libjpeg" = "xno"],
              [dependency_notice=yes
               AS_ECHO(["  - libjpeg (http://ijg.org), could not be linked with in your library"])
               AS_ECHO(["    search path, or is manually disabled. If JPEG inputs/outputs are"])
               AS_ECHO(["    requested, the respective tool will inform you and abort with an"])
               AS_ECHO(["    error."])
               AS_ECHO([]) ])

        AS_IF([test "x$has_libtiff" = "xno"],
              [dependency_notice=yes
               AS_ECHO(["  - libtiff (http://libtiff.maptools.org), could not be linked with in"])
               AS_ECHO(["    your library search path, or is manually disabled. If TIFF"])
               AS_ECHO(["    inputs/outputs are requested, the respective tool will inform"])
               AS_ECHO(["    you and abort with an error."])
               AS_ECHO([]) ])

        AS_IF([test "x$has_libgit2" = "x0"],
              [dependency_notice=yes
               AS_ECHO(["  - libgit2 (https://libgit2.org), could not be linked with in your"])
               AS_ECHO(["    library search path, or is manually disabled. When present, Git's"])
               AS_ECHO(["    describe output will be stored in the output files if Gnuastro's"])
               AS_ECHO(["    programs were called within a Gitversion controlled directory to"])
               AS_ECHO(["    help in reproducibility."])
               AS_ECHO([]) ])

        AS_IF([test "x$usable_libtool" = "xno"],
              [dependency_notice=yes
               AS_ECHO(["  - GNU Libtool (https://www.gnu.org/s/libtool) can't be used on this"])
               AS_ECHO(["    system (see below). Gnuastro's BuildProgram uses GNU libtool to"])
               AS_ECHO(["    link your source code with the various libraries (Gnuastro's"])
               AS_ECHO(["    dependencies). Therefore BuildProgram will not be built or installed."])
               AS_ECHO(["    Please note that not having GNU libtool in your search path will not"])
               AS_ECHO(["    harm the rest of Gnuastro's building and installation. Gnuastro has"])
               AS_ECHO(["    its own internal implementation of GNU Libtool to build its self. This"])
               AS_ECHO(["    warning is only to let you know that BuildProgram will not be"])
               AS_ECHO(["    part of this build. The executable names searched were 'libtool'"])
               AS_ECHO(["    and 'glibtool'. The shells searched were 'sh', 'bash' and 'zsh'."])
               AS_ECHO([])
               AS_IF([test "x$has_gnulibtool" = "xyes"],
                     [AS_ECHO(["    -- GNU Libtool is present, but couldn't be run in tested shells."])
                      AS_IF([test "x$bash_version" = x], [junk=1],
                            [AS_ECHO(["       (GNU Bash was found but not a suitable version: $bash_version)"]) ])
                      AS_ECHO([])],
                     [AS_IF([test "x$has_libtool" = "xyes"],
                            [AS_ECHO(["    -- A libtool implementation was found, but it isn't GNU."])
                             AS_ECHO([]) ]) ])
              ])

        AS_IF([test "x$has_ghostscript" = "xno"],
              [dependency_notice=yes
               AS_ECHO(["  - GPL GhostScript (https://www.ghostscript.com) version 9.10 or later,"])
               AS_ECHO(["    with the executable name 'gs', was not found in your PATH environment"])
               AS_ECHO(["    variable. If PDF outputs are desired, the respective tool it will abort"])
               AS_ECHO(["    with an EPS output which you can convert to PDF by other means."])
               AS_ECHO([]) ])

        # Notice for obtaining the optional dependencies using a package
        # manager.
        AS_IF([test "x$dependency_notice" = "xyes"],
              [AS_ECHO(["  You can use your package manager for easy and fast installation of all"])
               AS_ECHO(["  the mandatory and optional dependencies in one command. See the link"])
               AS_ECHO(["  below:"])
               AS_ECHO(["    https://www.gnu.org/s/gnuastro/manual/html_node/Dependencies-from-package-managers.html"])
               AS_ECHO([])
               AS_ECHO(["  All checks related to the warning(s) above will be skipped."])
               AS_ECHO([])
              ])

        # Notice about PATH: The last two scenarios described below are
        # taken from
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
               AS_ECHO(["    a './' before it's name. Run the following command after"])
               AS_ECHO(["    installing Gnuastro to learn more about PATH:"])
               AS_ECHO(["        $ info gnuastro \"Installation directory\""])
               AS_ECHO([]) ])
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
