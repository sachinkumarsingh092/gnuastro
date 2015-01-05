/*********************************************************************
Functions to check and set command line argument values and files.
This is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "checkset.h"


#define FIXEDFORFILE        "The value of `%s` (given as %s)"
#define FIXEDFOROPTION      "The value of `--%s (-%c)` (given as %s)"
#define NOTEMSG_NOTNUMBER   "not readable as number."
#define NOTEMSG_EQZEROORONE "should be ==0 or ==1."
#define NOTEMSG_LEQZERO     "should be >=0."
#define NOTEMSG_LARGERZERO  "should be >0."
#define NOTEMSG_4OR8        "should be either 4 or 8."
#define NOTEMSG_SMALLERONE  "should be >=0 and <=1."


#define CHECKFULLNUMBER   if(*tailptr!='\0'){				\
    if(filename)							\
      error_at_line(EXIT_FAILURE, 0, filename, lineno,			\
		    FIXEDFORFILE" "NOTEMSG_NOTNUMBER, lo, optarg);	\
    else								\
      error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_NOTNUMBER,	\
	    lo, so, optarg);						\
  }





/****************************************************************
 *****************      Check the numbers    ********************
 ****************************************************************/
void
intzeroorone(char *optarg, int *var, char *lo, char so, char *spack,
	     char *filename, size_t lineno)
{
  long tmp;
  char *tailptr;

  *var=tmp=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
  if(tmp!=0 && tmp!=1)
    {
      if(filename)
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      FIXEDFORFILE" "NOTEMSG_EQZEROORONE, lo, optarg);
      else
	error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_EQZEROORONE,
	      lo, so, optarg);
    }
}





void
intelzero(char *optarg, int *var, char *lo, char so, char *spack,
	  char *filename, size_t lineno)
{
  long tmp;
  char *tailptr;
  *var=tmp=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
  if(tmp<0)
    {
      if(filename)
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      FIXEDFORFILE" "NOTEMSG_LEQZERO, lo, optarg);
      else
	error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_LEQZERO,
	      lo, so, optarg);
    }
}





void
intlzero(char *optarg, int *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno)
{
  long tmp;
  char *tailptr;
  *var=tmp=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
  if(tmp<=0)
    {
      if(filename)
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      FIXEDFORFILE" "NOTEMSG_LARGERZERO, lo, optarg);
      else
	error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_LARGERZERO,
	      lo, so, optarg);
    }
}





void
int4or8(char *optarg, int *var, char *lo, char so, char *spack,
	char *filename, size_t lineno)
{
  long tmp;
  char *tailptr;
  *var=tmp=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
  if(tmp!=4 && tmp!=8)
    {
      if(filename)
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      FIXEDFORFILE" "NOTEMSG_4OR8, lo, optarg);
      else
	error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_4OR8,
	      lo, so, optarg);
    }
}





void
anylong(char *optarg, long *var, char *lo, char so, char *spack,
	char *filename, size_t lineno)
{
  char *tailptr;
  *var=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
}





/* size_t is always >=0. But if the user gives a negative value, it
   will become an extremely large number (because size_t is
   unsigned). So we read the user's value in long, which accepts
   negative values, then we check that value before placing it in the
   size_t pointer.*/
void
sizetelzero(char *optarg, size_t *var, char *lo, char so, char* spack,
	    char *filename, size_t lineno)
{
  long tmp;
  char *tailptr;
  *var=tmp=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
  if(tmp<0)
    {
      if(filename)
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      FIXEDFORFILE" "NOTEMSG_LEQZERO, lo, optarg);
      else
	error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_LEQZERO,
	      lo, so, optarg);
    }
}




/* See explanation for sizetlzero above. */
void
sizetlzero(char *optarg, size_t *var, char *lo, char so, char* spack,
	   char *filename, size_t lineno)
{
  long tmp;
  char *tailptr;
  *var=tmp=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
  if(tmp<=0)
    {
      if(filename)
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      FIXEDFORFILE" "NOTEMSG_LARGERZERO, lo, optarg);
      else
	error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_LARGERZERO,
	      lo, so, optarg);
    }
}





void
floatl0(char *optarg, float *var, char *lo, char so, char* spack,
	char *filename, size_t lineno)
{
  float tmp;
  char *tailptr;
  *var=tmp=strtof(optarg, &tailptr);

  CHECKFULLNUMBER;
  if(tmp<=0)
    {
      if(filename)
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      FIXEDFORFILE" "NOTEMSG_LARGERZERO, lo, optarg);
      else
	error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_LARGERZERO,
	      lo, so, optarg);
    }
}





void
floatl0s1(char *optarg, float *var, char *lo, char so, char* spack,
	  char *filename, size_t lineno)
{
  float tmp;
  char *tailptr;
  *var=tmp=strtof(optarg, &tailptr);

  CHECKFULLNUMBER;
  if(tmp>1.0f || tmp<0)
    {
      if(filename)
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      FIXEDFORFILE" "NOTEMSG_SMALLERONE, lo, optarg);
      else
	error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_SMALLERONE,
	      lo, so, optarg);
    }
}




void
anyfloat(char *optarg, float *var, char *lo, char so, char *spack,
	 char *filename, size_t lineno)
{
  char *tailptr;
  *var=strtof(optarg, &tailptr);

  CHECKFULLNUMBER;
}




void
doublel0(char *optarg, double *var, char *lo, char so, char* spack,
	 char *filename, size_t lineno)
{
  double tmp;
  char *tailptr;
  *var=tmp=strtod(optarg, &tailptr);

  CHECKFULLNUMBER;
  if(tmp<=0)
    {
      if(filename)
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      FIXEDFORFILE" "NOTEMSG_LARGERZERO, lo, optarg);
      else
	error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_LARGERZERO,
	      lo, so, optarg);
    }
}





void
anydouble(char *optarg, double *var, char *lo, char so, char *spack,
	  char *filename, size_t lineno)
{
  char *tailptr;
  *var=strtod(optarg, &tailptr);

  CHECKFULLNUMBER;
}




















/**************************************************************/
/**********          My String functions:          ************/
/**************************************************************/
int
stringhasspace(char *in)
{
  do
    switch(*in)
      {
      case ' ': case '\t': case '\v':
	return 1;
      }
  while(*(++in)!='\0');
  return 0;
}





char *
malloccat(char *inname, char *toappend)
{
  char *out;
  size_t inl, apl;

  inl=strlen(inname);
  apl=strlen(toappend);

  errno=0;
  out=malloc(inl+apl+1);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "Allocating %lu bytes in malloccat",
	  inl+apl+1);

  strcpy(out, inname);
  strcat(out, toappend);
  return out;
}


















/**************************************************************/
/********** Set file names and check if they exist ************/
/**************************************************************/
/* Check if a file exists and report if it doesn't: */
void
checkfile(char *filename)
{
  FILE *tmpfile;
  errno=0;
  tmpfile = fopen(filename, "r");
  if(tmpfile)
    {
      errno=0;
      if(fclose(tmpfile)==EOF)
	error(EXIT_FAILURE, errno, filename);
    }
  else
    error(EXIT_FAILURE, errno, filename);
}





/* Check if a file exists. If so, remove it. */
void
checkremovefile(char *filename, int dontdelete)
{
  FILE *tmpfile;

  /* We want to make sure that we can open and write to this file. But
     the user might have asked to not delete the file, so the
     contents should not be changed. Therefore we have to open it with
     `r+`. */
  errno=0;
  tmpfile=fopen(filename, "r+");
  if (tmpfile)
    {
      /* Close the file and make sure that it should be deleted. */
      errno=0;
      if(fclose(tmpfile))
	error(EXIT_FAILURE, errno, filename);
      if(dontdelete)
	error(EXIT_FAILURE, 0, "%s already exists and you have "
	      "asked to not remove it with the `--dontdelete` "
	      "(`-D`) option.", filename);

      /* Delete the file: */
      errno=0;
      if(unlink(filename))
	error(EXIT_FAILURE, errno, filename);
    }
  /* If the file doesn't exist, there is no problem, we wanted to
     remove it any way! Any other kind of error should not be
     tolerated! */
  else if(errno!=ENOENT)
    error(EXIT_FAILURE, errno, filename);
}





/* Question that will be answered: is this name, the name of a
   potential (maybe non-existant) file or is it a directory? If it is
   neither, it will abort and tell the user. */
int
nameisafile(char *name, int dontdelete)
{
  struct stat nameinfo;

  errno=0;
  if(stat(name, &nameinfo)!=0)
    {
      if(errno==ENOENT)	/* ENOENT: No such file or directory. */
	{/* Make the file temporarily and see if everything is ok. */
	  checkremovefile(name, dontdelete);
	  return 1;		/* It will be a file, GOOD */
	}
      else			/* Some strange condition, ABORT */
	error(EXIT_FAILURE, errno, name);
    }

  if(S_ISDIR(nameinfo.st_mode))	/* It is a directory, BAD */
    return 0;
  else if (S_ISREG(nameinfo.st_mode)) /* It is a file, GOOD. */
    {
      /* Delete the file if it should be deleted. */
      if(dontdelete)
	error(EXIT_FAILURE, 0, "%s already exists and you have "
	      "asked to not remove it with the `--dontdelete` "
	      "(`-D`) option.", name);
      errno=0;
      if(unlink(name))
	error(EXIT_FAILURE, errno, "%s", name);
      return 1;
    }
  else 				/* Not a file or a dir, ABORT */
    error(EXIT_FAILURE, 0, "%s not a file or a directory.", name);

  return 0;
}





/* Change the ending of a name. If the name has a `.ext` postfix, then
   it will be changed to the string in `append`. If it doesn't,
   append will be added to the name. */
void
automaticoutput(char *inname, char *postfix, int removedirinfo,
		int dontdelete, char **outname)
{
  char *out;
  size_t i, l, offset=0;

  /* Check if the input file is actually a readable file or not! */
  checkfile(inname);

  /* Allocate space for the output name. It has been allocated before
     (it was either read from the command line or a default file), so
     free the space it already has.*/
  if(*outname) free(*outname);
  *outname=out=malloccat(inname, postfix);

  /* Put the input in the space and remove all elements including and
     after '.'. Note that if there is no '.' in the name, malloccat
     has already appended inname and postfix.*/
  l=strlen(inname);
  strcpy(out, inname);
  for(i=l;i!=0;--i)
    {
      /* We don't want to remove any '.' in a directory name so if a
	 '/' is encountered in our search from the end of the file
	 name, we won't continue. */
      if(out[i]=='/') break;
      else if(out[i]=='.')
	{
	  out[i]='\0';
	  strcat(out, postfix);
	  break;
	}
    }

  /* If it is desired to remove the directory information from the
     name, do it here. Some unused space will remain after removing
     the directory information, but that can be ignored, since it
     can't be too much. */
  if(removedirinfo)
    {
      l=strlen(out);
      for(i=l;i!=0;--i)	 	  /* Find the last forward slash.      */
	if(out[i]=='/')
	  {offset=i+1; break;}
      if(offset)
	for(i=offset;i<=l;++i)	/* <= because we want to shift the   */
	  out[i-offset]=out[i]; /* '\0' character in the string too. */
    }

  /* Remove the output file if it already exits. */
  checkremovefile(out, dontdelete);
}





/* Check if dirname is actually a real directory and that we can
   actually write inside of it. To insure all conditions an actual
   file will be made */
void
checkdirwriteaddslash(char **dirname)
{
  int file_d;
  char *tmpname, *indir=*dirname, buf[]="A test file.";

  /* Set the template for the temporary file: */
  if(indir[strlen(indir)-1]=='/')
    tmpname=malloccat(indir, "astrutilsXXXXXX");
  else
    tmpname=malloccat(indir, "/astrutilsXXXXXX");

  /* Make a temporary file name and try openning it. */
  errno=0;
  file_d=mkstemp(tmpname);
  if(file_d==-1)
    error(EXIT_FAILURE, errno, "Cannot write output in the directory %s",
	  indir);
  write(file_d, buf, strlen(buf));
  close(file_d);

  /* Delete the temporary file: */
  errno=0;
  if(unlink(tmpname)==-1)
    error(EXIT_FAILURE, errno, "Temporary file name made for testing "
	  "the %s directory (%s) could not be deleted", indir, tmpname);

  /* Remove the extra characters that were added for the random name. */
  tmpname[strlen(tmpname)-15]='\0';

  free(*dirname);
  *dirname=tmpname;
}
