/*********************************************************************
Functions to check and set command line argument values and files.
This is part of GNU Astronomy Utilities (gnuastro) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef CHECKSET_H
#define CHECKSET_H

#include <math.h>
#include <fitsio.h>


/**************************************************************/
/*********                 Macros                **************/
/**************************************************************/
#define CHECKCOLINCAT(INCOL,NAME) {					\
    size_t i;								\
    									\
    if( (INCOL) >= p->cs1 )						\
      error(EXIT_FAILURE, 0, "%s only has %lu columns while you "	\
	    "have requested column %lu (counting from zero) for "	\
	    "`--%s`.", p->up.catname, p->cs1, (INCOL), (NAME));		\
									\
    for(i=0;i<p->cs0;++i)						\
      if( isnan(p->cat[i*p->cs1+(INCOL)]) )				\
	error(EXIT_FAILURE, 0, "%s: Column %lu (--%s) in row %lu "	\
	      "could not be read as a number. See %s. Note that "	\
	      "counting starts from zero.",				\
	      p->up.catname, (INCOL), (NAME), i, TXTARRAYVVLOG);	\
  }




















/**************************************************************/
/********* Read arguments and check their values **************/
/**************************************************************/
void
intzeroorone(char *optarg, int *var, char *lo, char so, char* spack,
	     char *filename, size_t lineno);

void
intelzero(char *optarg, int *var, char *lo, char so, char *spack,
	  char *filename, size_t lineno);

void
intlzero(char *optarg, int *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno);

void
int4or8(char *optarg, int *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno);

void
intsmallerequalto(char *optarg, int *var, char *lo, char so, char *spack,
                  char *filename, size_t lineno, long maxvalue);

void
anylong(char *optarg, long *var, char *lo, char so, char *spack,
	char *filename, size_t lineno);

void
sizetelzero(char *optarg, size_t *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno);

void
sizetlzero(char *optarg, size_t *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno);

void
floatl0(char *optarg, float *var, char *lo, char so, char *spack,
	char *filename, size_t lineno);

void
floatl0s1(char *optarg, float *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno);

void
anyfloat(char *optarg, float *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno);

void
doublel0(char *optarg, double *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno);

void
anydouble(char *optarg, double *var, char *lo, char so, char *spack,
	  char *filename, size_t lineno);










/**************************************************************/
/**********          My String functions:          ************/
/**************************************************************/
int
stringhasspace(char *in);

char *
malloccat(char *inname, char *toappend);










/**************************************************************/
/********** Set file names and check if they exist ************/
/**************************************************************/
void
checkfile(char *filename);

void
checkremovefile(char *filename, int dontdelete);

int
dir0file1(char *name, int dontdelete);

void
automaticoutput(char *inname, char *suffix, int removedirinfo,
		int dontdelete, char **outname);

void
checkdirwriteaddslash(char **dirname);

#endif
