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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include <fitsio.h>

#include "checkset.h"


#define FIXEDFORFILE        "The value of `%s` (given as %s)"
#define FIXEDFOROPTION      "The value of `--%s (-%c)` (given as %s)"
#define NOTEMSG_NOTNUMBER   "not readable as number"
#define NOTEMSG_EQZEROORONE "should be ==0 or ==1"
#define NOTEMSG_LEQZERO     "should be >=0"
#define NOTEMSG_LARGERZERO  "should be >0"
#define NOTEMSG_4OR8        "should be either 4 or 8"
#define NOTEMSG_SMALLERONE  "should be >=0 and <=1"


#define CHECKFULLNUMBER   if(*tailptr!='\0'){                           \
    if(filename)                                                        \
      error_at_line(EXIT_FAILURE, 0, filename, lineno,                  \
                    FIXEDFORFILE" "NOTEMSG_NOTNUMBER, lo, optarg);      \
    else                                                                \
      error(EXIT_FAILURE, 0, FIXEDFOROPTION" "NOTEMSG_NOTNUMBER,        \
            lo, so, optarg);                                            \
  }




/****************************************************************
 ************      Check and convert strings    *****************
 ****************************************************************/
/* See if a given string is a floating point number, if so, put it in
   the number.*/
int
gal_checkset_str_is_double(char *string, double *out)
{
  char *tailptr;
  double tmp=*out;
  *out=strtod(string, &tailptr);

  /* If the tail pointer (tailptr) is the string NULL character, then
     the string was a single double precision floating point
     number. However, if it is any other character, then put the
     initial value of out back inside of it and return 0. */
  if(*tailptr=='\0') return 1;
  *out=tmp;
  return 0;
}




















/****************************************************************
 *****************      Check the numbers    ********************
 ****************************************************************/
void
gal_checkset_int_zero_or_one(char *optarg, int *var, char *lo, char so,
                             char *spack, char *filename, size_t lineno)
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
gal_checkset_int_4_or_8(char *optarg, int *var, char *lo, char so,
                        char *spack, char *filename, size_t lineno)
{
  long tmp;
  char *tailptr;

  *var=tmp=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
  if(tmp!=4 && tmp!=8)
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
gal_checkset_int_el_zero(char *optarg, int *var, char *lo, char so,
                         char *spack, char *filename, size_t lineno)
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
gal_checkset_int_l_zero(char *optarg, int *var, char *lo, char so,
                        char *spack, char *filename, size_t lineno)
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
gal_checkset_int_smaller_equal_to(char *optarg, int *var, char *lo, char so,
                                  char *spack, char *filename, size_t lineno,
                                  long maxvalue)
{
  long tmp;
  char *tailptr;
  *var=tmp=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
  if(tmp>maxvalue)
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
gal_checkset_long_el_zero(char *optarg, long *var, char *lo, char so,
                          char *spack, char *filename, size_t lineno)
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
gal_checkset_any_long(char *optarg, long *var, char *lo, char so,
                      char *spack, char *filename, size_t lineno)
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
gal_checkset_sizet_el_zero(char *optarg, size_t *var, char *lo, char so,
                           char* spack, char *filename, size_t lineno)
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




/* See explanation for gal_checkset_sizet_l_zero above. */
void
gal_checkset_sizet_l_zero(char *optarg, size_t *var, char *lo, char so,
                          char* spack, char *filename, size_t lineno)
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





/* Positive and odd size_t. */
void
gal_checkset_sizet_p_odd(char *optarg, size_t *var, char *lo, char so,
                         char* spack, char *filename, size_t lineno)
{
  long tmp;
  char *tailptr;
  *var=tmp=strtol(optarg, &tailptr, 0);

  CHECKFULLNUMBER;
  if(tmp<0 || tmp%2==0)
    {
      if(filename)
        error_at_line(EXIT_FAILURE, 0, filename, lineno,
                      FIXEDFORFILE" should be >0 and odd.", lo, optarg);
      else
        error(EXIT_FAILURE, 0, FIXEDFOROPTION" should be >0 and odd",
              lo, so, optarg);
    }
}





void
gal_checkset_float_l_0(char *optarg, float *var, char *lo, char so,
                       char* spack, char *filename, size_t lineno)
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
gal_checkset_float_l_0_s_1(char *optarg, float *var, char *lo, char so,
                         char* spack, char *filename, size_t lineno)
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
gal_checkset_any_float(char *optarg, float *var, char *lo, char so,
                       char *spack, char *filename, size_t lineno)
{
  char *tailptr;
  *var=strtof(optarg, &tailptr);

  CHECKFULLNUMBER;
}




void
gal_checkset_double_l_0(char *optarg, double *var, char *lo, char so,
                        char* spack, char *filename, size_t lineno)
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
gal_checkset_double_el_0(char *optarg, double *var, char *lo, char so,
                         char* spack, char *filename, size_t lineno)
{
  double tmp;
  char *tailptr;
  *var=tmp=strtod(optarg, &tailptr);

  CHECKFULLNUMBER;
  if(tmp<0)
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
gal_checkset_double_l_value(char *optarg, double *var, char *lo, char so,
                            char* spack, double value, char *filename,
                            size_t lineno)
{
  float tmp;
  char *tailptr;
  *var=tmp=strtof(optarg, &tailptr);

  CHECKFULLNUMBER;
  if(tmp<=value)
    {
      if(filename)
        error_at_line(EXIT_FAILURE, 0, filename, lineno,
                      FIXEDFORFILE" should be > %.4f", lo, optarg, value);
      else
        error(EXIT_FAILURE, 0, FIXEDFOROPTION" should be > %.4f",
              lo, so, optarg, value);
    }
}





void
gal_checkset_any_double(char *optarg, double *var, char *lo, char so,
                        char *spack, char *filename, size_t lineno)
{
  char *tailptr;
  *var=strtod(optarg, &tailptr);

  CHECKFULLNUMBER;
}




















/**************************************************************/
/**********          Check fixed strings           ************/
/**************************************************************/
/* Check if the value to the `--type' option is recognized, if so set the
   integer value. */
void
gal_checkset_known_types(char *optarg, int *bitpix, char *filename,
                         size_t lineno)
{
  /* First check if the value is one of the accepted types. */
  if     (strcmp(optarg, "byte")==0)     *bitpix=BYTE_IMG;
  else if(strcmp(optarg, "short")==0)    *bitpix=SHORT_IMG;
  else if(strcmp(optarg, "long")==0)     *bitpix=LONG_IMG;
  else if(strcmp(optarg, "longlong")==0) *bitpix=LONGLONG_IMG;
  else if(strcmp(optarg, "float")==0)    *bitpix=FLOAT_IMG;
  else if(strcmp(optarg, "double")==0)   *bitpix=DOUBLE_IMG;
  else
    {
      if(filename)
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "given value of "
                      "the `type' option (`%s') is not recognized. It must "
                      "be `byte', `short', `long', `longlong', `float', or "
                      "`double'.", optarg);
      else
        error(EXIT_FAILURE, 0, "given value of the `--type' (`-T') option "
              "(`%s') is not recognized. It must be `byte', `short', `long' "
              "`longlong', `float', or `double'.", optarg);
    }
}



















/**************************************************************/
/**********          My String functions:          ************/
/**************************************************************/
int
gal_checkset_string_has_space(char *in)
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
gal_checkset_malloc_cat(char *inname, char *toappend)
{
  char *out;
  size_t inl, apl;

  inl=strlen(inname);
  apl=strlen(toappend);

  errno=0;
  out=malloc(inl+apl+1);
  if(out==NULL)
    error(EXIT_FAILURE, errno,
          "allocating %zu bytes in gal_checkset_malloc_cat", inl+apl+1);

  strcpy(out, inname);
  strcat(out, toappend);
  return out;
}




/* Copy the input string to the output (and also allocate the
   output. */
void
gal_checkset_allocate_copy(char *arg, char **copy)
{
  /* Allocate the necessary space: */
  errno=0;
  *copy=malloc(strlen(arg)+1);
  if(*copy==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes to copy %s",
          strlen(arg)+1, arg);
  strcpy(*copy, arg);
}




/* This function is mainly for reading in the arguments (from the
   command line or configuration files) that need to be copied. The
   set argument is for making sure that it has not already been set
   before, see the main.h files of any program. */
void
gal_checkset_allocate_copy_set(char *arg, char **copy, int *set)
{
  /* Incase *set==1, then you shouldn't do anything, just return. */
  if(*set) return;

  /* The variable was not copied, copy it: */
  errno=0;
  *copy=malloc(strlen(arg)+1);
  if(*copy==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes to copy %s",
          strlen(arg)+1, arg);
  strcpy(*copy, arg);
  *set=1;
}




















/**************************************************************/
/********** Set file names and check if they exist ************/
/**************************************************************/
/* Check if a file exists and report if it doesn't: */
void
gal_checkset_check_file(char *filename)
{
  FILE *tmpfile;
  errno=0;
  tmpfile = fopen(filename, "r");
  if(tmpfile)                        /* The file opened. */
    {
      if(fclose(tmpfile)==EOF)
        error(EXIT_FAILURE, errno, "%s", filename);
    }
  else
    error(EXIT_FAILURE, errno, "%s", filename);
}





/* Similar to `gal_checkset_check_file', but will report the result instead
   of doing it quietly. */
int
gal_checkset_check_file_report(char *filename)
{
  FILE *tmpfile;
  errno=0;
  tmpfile = fopen(filename, "r");
  if(tmpfile)                        /* The file opened. */
    {
      if(fclose(tmpfile)==EOF)
        error(EXIT_FAILURE, errno, "%s", filename);
      return 1;
    }
  else
    return 0;
}





/* Check if a file exists. If so, remove it. */
void
gal_checkset_check_remove_file(char *filename, int dontdelete)
{
  FILE *tmpfile;

  /* We want to make sure that we can open and write to this file. But
     the user might have asked to not delete the file, so the
     contents should not be changed. Therefore we have to open it with
     `r+`. */
  errno=0;
  tmpfile=fopen(filename, "r+");
  if (tmpfile)                        /* The file opened. */
    {
      /* Close the file and make sure that it should be deleted. */
      errno=0;
      if(fclose(tmpfile))
        error(EXIT_FAILURE, errno, "%s", filename);
      if(dontdelete)
        error(EXIT_FAILURE, 0, "%s already exists and you have "
              "asked to not remove it with the `--dontdelete` "
              "(`-D`) option", filename);

      /* Delete the file: */
      errno=0;
      if(unlink(filename))
        error(EXIT_FAILURE, errno, "%s", filename);
    }
  /* If the file doesn't exist, there is no problem, we wanted to
     remove it any way! Any other kind of error should not be
     tolerated! */
  else if(errno!=ENOENT)
    error(EXIT_FAILURE, errno, "%s", filename);
}





/* Check output file name: If a file exists or can exist and can be
   written to, this function will return 1. If not (for example it is
   a directory) it will return 0. Finally, if it exists but cannot be
   deleted, report an error and abort. */
int
gal_checkset_dir_0_file_1(char *name, int dontdelete)
{
  FILE *tmpfile;
  struct stat nameinfo;

  if(name==NULL)
    error(EXIT_FAILURE, 0, "a bug! The input to the "
          "gal_checkset_dir_0_file_1 function in checkset.c should not "
          "be NULL. Please contact us at "PACKAGE_BUGREPORT" so we can "
          "see what went wrong and fix it in future updates");

  errno=0;
  if(stat(name, &nameinfo)!=0)
    {
      if(errno==ENOENT)        /* ENOENT: No such file or directory. */
        {/* Make the file temporarily and see if everything is ok. */
          errno=0;
          tmpfile=fopen(name, "w");
          if (tmpfile)
            {
              fprintf(tmpfile, "Only to test write access.");
              errno=0;
              if(fclose(tmpfile))
                error(EXIT_FAILURE, errno, "%s", name);
              errno=0;
              if(unlink(name))
                error(EXIT_FAILURE, errno, "%s", name);
            }
          else
            error(EXIT_FAILURE, errno, "%s", name);
          return 1;                    /* It is a file name, GOOD */
        }
      else                             /* Some strange condition, ABORT */
        error(EXIT_FAILURE, errno, "%s", name);
    }

  if(S_ISDIR(nameinfo.st_mode))        /* It is a directory, BAD */
    return 0;
  else if (S_ISREG(nameinfo.st_mode))  /* It is a file, GOOD. */
    {
      gal_checkset_check_remove_file(name, dontdelete);
      return 1;
    }
  else                                 /* Not a file or a dir, ABORT */
    error(EXIT_FAILURE, 0, "%s not a file or a directory", name);

  error(EXIT_FAILURE, 0, "a bug! In gal_checkset_dir_0_file_1, (in "
        "checkset.c). The process should not reach the end of the "
        "function! Please contact us at "PACKAGE_BUGREPORT" so we can "
        "see what went wrong and fix it in future updates");
  return 0;
}





/* Allocate space and write the output name (outname) based on a given
   input name (inname). The suffix of the input name (if present) will
   be removed and the given suffix will be put in the end. */
void
gal_checkset_automatic_output(char *inname, char *suffix, int removedirinfo,
                              int dontdelete, char **outname)
{
  char *out;
  size_t i, l, offset=0;

  /* Merge the contents of the input name and suffix name (while also
     allocating the necessary space).*/
  out=gal_checkset_malloc_cat(inname, suffix);

  /* If there is actually a suffix, replace it with the (possibly) existing
     suffix. */
  if(suffix[0]!='\0')
    {
      /* Start from the end of the input array*/
      l=strlen(inname);
      for(i=l-1;i!=0;--i)
        {
          /* We don't want to touch anything before a `/' (directory
             names). We are only concerned with file names here. */
          if(out[i]=='/')
            {
              /* When `/' is the last input character, then the input is
                 clearly not a filename, but a directory name. In this
                 case, adding a suffix is meaningless (a suffix belongs to
                 a filename for Gnuastro's tools). So close the string
                 after the `/' and leave the loop. However, if the `/'
                 isn't the last input name charector, there is probably a
                 filename (without a "." suffix), so break from the
                 loop. No further action is required, since we initially
                 allocated the necessary space and concatenated the input
                 and suffix arrays. */
              if(i==l-1)
                out[i+1]='\0';
              break;
            }

          /* The input file names can be compressed names (for example
             `.fits.gz'). Currently the only compressed formats
             (decompressed within CFITSIO) are listed in
             `gal_fits_name_is_fits' and `gal_fits_suffix_is_fits'.*/
          else if(out[i]=='.' && !( ( out[i+1]=='g' && out[i+2]=='z' )
                                    || (out[i+1]=='f' && out[i+2]=='z' )
                                    || out[i+1]=='Z' ) )
            {
              out[i]='\0';
              strcat(out, suffix);
              break;
            }
        }
    }

  /* If it is desired to remove the directory information from the
     name, do it here. Some unused space will remain after removing
     the directory information, but that can be ignored, since it
     can't be too much. */
  if(removedirinfo)
    {
      l=strlen(out);
      for(i=l;i!=0;--i)         /* Find the last forward slash.      */
        if(out[i]=='/')
          {offset=i+1; break;}
      if(offset)
        for(i=offset;i<=l;++i)  /* <= because we want to shift the   */
          out[i-offset]=out[i]; /* '\0' character in the string too. */
    }

  /* Remove the created filename if it already exits. */
  gal_checkset_check_remove_file(out, dontdelete);

  /* Free the outname if it was already allocated before. */
  free(*outname);
  *outname=out;
}





/* Given a filename, this function will separate its directory name
   part. */
char *
gal_checkset_dir_part(char *input)
{
  char *out;
  size_t i, l;

  /* Find the first slash. */
  l=strlen(input);
  for(i=l;i!=0;--i)
    if(input[i]=='/')
      break;

  /* If there was no slash, then the current directory should be
     given: */
  if(i==0)
    {
      errno=0;
      out=malloc(3*sizeof *out);
      if(out==NULL)
        error(EXIT_FAILURE, errno, "%zu bytes for current directory "
              "gal_checkset_dir_part", 3*sizeof *out);
      strcpy(out, "./");
    }
  else
    {
      errno=0;
      out=malloc((l+1)*sizeof *out);
      if(out==NULL)
        error(EXIT_FAILURE, errno, "%zu bytes for gal_checkset_dir_part",
              (l+1)*sizeof *out);
      strcpy(out, input);
      out[i+1]='\0';
    }

  return out;
}





/* Given a file name, keep the non-directory part. Note that if there
   is no forward slash in the input name, the full input name is
   considered to be the notdir output.*/
char *
gal_checkset_not_dir_part(char *input)
{
  size_t i, l;
  char *out, *tmp=input;

  /* Find the first `/' to identify the directory */
  l=strlen(input);
  for(i=l;i!=0;--i)
    if(input[i]=='/')
      { tmp=&input[i+1]; break; }

  /* Get the length of the notdir name: */
  l=strlen(tmp);
  errno=0;
  out=malloc((l+1)*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for notdir", (l+1)*sizeof *out);

  strcpy(out, tmp);
  return out;
}





/* Check if dirname is actually a real directory and that we can
   actually write inside of it. To insure all conditions an actual
   file will be made */
void
gal_checkset_check_dir_write_add_slash(char **dirname)
{
  int file_d;
  char *tmpname, *indir=*dirname/*, buf[]="A test"*/;

  /* Set the template for the temporary file: */
  if(indir[strlen(indir)-1]=='/')
    tmpname=gal_checkset_malloc_cat(indir, "gnuastroXXXXXX");
  else
    tmpname=gal_checkset_malloc_cat(indir, "/gnuastroXXXXXX");

  /* Make a temporary file name and try openning it. */
  errno=0;
  file_d=mkstemp(tmpname);
  if(file_d==-1)
    error(EXIT_FAILURE, errno, "cannot write output in the directory %s",
          indir);
  /*
  errno=0;
  printf("\n\n%s\n\n", tmpname);
  if( write(file_d, buf, strlen(buf)) == -1 )
    error(EXIT_FAILURE, errno, "%s: writing to this temporary file to "
          "check the given `%s` directory", tmpname, indir);
  */
  errno=0;
  if( close(file_d) == -1 )
    error(EXIT_FAILURE, errno, "%s: Closing this temporary file to check "
          "the given `%s` directory", tmpname, indir);

  /* Delete the temporary file: */
  errno=0;
  if(unlink(tmpname)==-1)
    error(EXIT_FAILURE, errno, "%s: removing this temporary file made "
          "to check the given `%s directory`", tmpname, indir);

  /* Remove the extra characters that were added for the random name. */
  tmpname[strlen(tmpname)-14]='\0';

  free(*dirname);
  *dirname=tmpname;
}
