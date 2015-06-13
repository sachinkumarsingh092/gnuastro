/*********************************************************************
Functions to check and set command line argument values and files.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
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


#define CHECKMASKNAMEANDHDU(spack) {                                    \
    if(p->up.masknameset)                                               \
      {                                                                 \
        if(p->up.mhduset==0)                                            \
          error(EXIT_FAILURE, 0, "A mask image was specified (%s). "    \
                "However, no HDU is given for it. Please add a HDU "    \
                "for the mask with the `--mhdu' (`-H') option. If you " \
                "regularly use a mask, you may consider adding `mhdu' " \
                "to the %s configuration file. For more information, "  \
                "please see the `Configuration files' section of the "  \
                "%s manual by running ` info gnuastro ' on the "        \
                "command-line.",                                        \
                p->up.maskname, (spack), PACKAGE_NAME);                 \
        if(strcmp(p->up.inputname, p->up.maskname)==0)                  \
          {                                                             \
            if(strcmp(p->up.mhdu, p->cp.hdu)==0)                        \
              error(EXIT_FAILURE, 0, "The specified mask name and "     \
                    "input image name are the same while the input "    \
                    "image hdu name and mask hdu are also identical!"); \
          }                                                             \
      }                                                                 \
    else if(p->up.mhduset && strcmp(p->up.mhdu, p->cp.hdu))             \
      p->up.maskname=p->up.inputname;                                   \
    else                                                                \
      p->up.maskname=NULL;                                              \
  }


#define PRINTSTINGMAYBEWITHSPACE(name,string) {                         \
    if(stringhasspace(string))                                          \
      fprintf(fp, CONF_SHOWFMT"\"%s\"\n", name, string);                \
    else                                                                \
      fprintf(fp, CONF_SHOWFMT"%s\n", name, string);                    \
  }



















/**************************************************************/
/********* Read arguments and check their values **************/
/**************************************************************/
void
intzeroorone(char *optarg, int *var, char *lo, char so, char* spack,
	     char *filename, size_t lineno);

void
int4or8(char *optarg, int *var, char *lo, char so, char *spack,
        char *filename, size_t lineno);

void
intelzero(char *optarg, int *var, char *lo, char so, char *spack,
	  char *filename, size_t lineno);

void
intlzero(char *optarg, int *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno);

void
intsmallerequalto(char *optarg, int *var, char *lo, char so, char *spack,
                  char *filename, size_t lineno, long maxvalue);

void
longelzero(char *optarg, long *var, char *lo, char so, char *spack,
           char *filename, size_t lineno);

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
sizetpodd(char *optarg, size_t *var, char *lo, char so, char* spack,
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
doublele0(char *optarg, double *var, char *lo, char so, char* spack,
          char *filename, size_t lineno);

void
doublelvalue(char *optarg, double *var, char *lo, char so, char* spack,
             double value, char *filename, size_t lineno);

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

void
allocatecopyset(char *arg, char **copy, int *set);








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
