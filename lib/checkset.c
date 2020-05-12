/*********************************************************************
Functions to check and set command line argument values and files.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/fcntl.h>

#include <gnuastro/data.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>









/**************************************************************/
/**********               Environment              ************/
/**************************************************************/
/* The GSL random number generator (RNG) reads values from the
   environment. This function is designed to make the job easier for any
   program using GSL's RNG. If the user doesn't want to set the */
gsl_rng *
gal_checkset_gsl_rng(uint8_t envseed_bool, const char **name,
                     unsigned long int *seed)
{
  gsl_rng *rng;

  /* Let GSL read the environment and convert the type name (as string) to
     'gsl_rng_type'. After this function, 'gsl_rng_default' contains the
     generator's type and 'gsl_rng_default_seed' contains the (possibly)
     given seed.*/
  gsl_rng_env_setup();

  /* Allocate the random number generator based on the requested type and
     save its name. If no 'GSL_RNG_TYPE' is set, then use a fixed
     generator.*/
  rng=gsl_rng_alloc(secure_getenv("GSL_RNG_TYPE")
                    ? gsl_rng_default
                    : gsl_rng_ranlxs1);
  *name = gsl_rng_name(rng);

  /* Initialize the random number generator, depending on the
     'envseed_bool' argument. */
  *seed = ( envseed_bool
            ? gsl_rng_default_seed
            : gal_timing_time_based_rng_seed() );
  gsl_rng_set(rng, *seed);

  /* Return the GSL RNG structure. */
  return rng;
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
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes", __func__,
          inl+apl+1);

  strcpy(out, inname);
  strcat(out, toappend);
  return out;
}




/* Copy the input string to the output (and also allocate the
   output. */
void
gal_checkset_allocate_copy(const char *arg, char **copy)
{
  if(arg)
    {
      errno=0;
      *copy=malloc(strlen(arg)+1);
      if(*copy==NULL)
        error(EXIT_FAILURE, errno, "%s: %zu bytes to copy %s", __func__,
              strlen(arg)+1, arg);
      strcpy(*copy, arg);
    }
  else
    *copy=NULL;
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
    error(EXIT_FAILURE, errno, "%s: %zu bytes to copy %s", __func__,
          strlen(arg)+1, arg);
  strcpy(*copy, arg);
  *set=1;
}





/* The dataset may be alone in a file (for example a table in a text file)
   or it may an extension of a FITS file. In error messages in particular,
   we need to differentiate between the two. This function will check the
   filename and if it is FITS, it will return a string with the filename
   and HDU in parenthesis. If it isn't a FITS file, it will only return the
   filename. Note that the output needs to be freed, although when used in
   an error message, you can leave it to the system to free the
   space. There is no problem. */
char *
gal_checkset_dataset_name(char *filename, char *hdu)
{
  char *out;
  if( gal_fits_name_is_fits(filename) )
    {
      if( asprintf(&out, "%s (hdu %s)", filename, hdu)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
    }
  else
    gal_checkset_allocate_copy(filename, &out);
  return out;
}



















/**************************************************************/
/********** Set file names and check if they exist ************/
/**************************************************************/
/* Given a filename, this function will separate its directory name
   part. */
char *
gal_checkset_dir_part(char *filename)
{
  char *out;
  size_t i, l=strlen(filename);

  /* Find the first slash from the end. */
  for(i=l;i!=0;--i)
    if(filename[i]=='/')
      break;

  /* If there was no slash, then the current directory should be
     given: */
  if(i==0 && filename[0]!='/')
    gal_checkset_allocate_copy("./", &out);
  else
    {
      gal_checkset_allocate_copy(filename, &out);
      out[i+1]='\0';
    }

  return out;
}





/* Given a file name, keep the non-directory part. Note that if there
   is no forward slash in the input name, the full input name is
   considered to be the notdir output.*/
char *
gal_checkset_not_dir_part(char *filename)
{
  size_t i, l;
  char *out, *tmp=filename;

  /* Find the first '/' to identify the directory */
  l=strlen(filename);
  for(i=l;i!=0;--i)
    if(filename[i]=='/')
      { tmp=&filename[i+1]; break; }

  /* Get the length of the notdir name: */
  l=strlen(tmp);
  errno=0;
  out=malloc((l+1)*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%s: %zu bytes for notdir", __func__,
          (l+1)*sizeof *out);

  strcpy(out, tmp);
  return out;
}





/* Check if a file exists and report if it doesn't. */
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





/* Similar to 'gal_checkset_check_file', but will report the result instead
   of doing it quietly. */
int
gal_checkset_check_file_return(char *filename)
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





/* If a file doesn't exist and its directory is writable, return
   1. Otherwise, return 0. */
int
gal_checkset_writable_notexist(char *filename)
{
  int out=1;
  char *dir;
  FILE *tmpfile;

  /* If the filename is 'NULL' everything is ok (it doesn't exist)! In some
     cases, a NULL filename is interpretted to mean standard output. */
  if(filename==NULL) return 1;

  /* We want to make sure that we can open and write to this file. But
     the user might have asked to not delete the file, so the
     contents should not be changed. Therefore we have to open it with
     'r+'. */
  errno=0;
  tmpfile=fopen(filename, "r+");
  if (tmpfile)                        /* The file opened. */
    {
      /* Close the file. */
      errno=0;
      if(fclose(tmpfile))
        error(EXIT_FAILURE, errno, "%s", filename);

      /* The file exists, return 0. */
      out=0;
    }

  /* If the file doesn't exist, we just need to make sure if we have write
     permissions to its host directory. */
  else
    {
      /* Separate the directory part of the filename. */
      dir=gal_checkset_dir_part(filename);

      /* See if this directory is writable by this user. */
      errno=0;
      if( access(dir, W_OK) ) out=0;

      /* Clean up. */
      free(dir);
    }

  /* Return the final value. */
  return out;
}





/* Check if a file exists and can be opened. If the 'keep' value is
   non-zero, then the file will remain untouched, otherwise, it will be
   deleted (since most programs need to make a clean output). When the file
   is to be deleted and 'dontdelete' has a non-zero value, then the file
   won't be deleted, but the program will abort with an error, informing
   the user that the output can't be made. */
void
gal_checkset_writable_remove(char *filename, int keep, int dontdelete)
{
  char *dir;
  FILE *tmpfile;

  /* If the filename is 'NULL' everything is ok (it doesn't exist)! In some
     cases, a NULL filename is interpretted to mean standard output. */
  if(filename==NULL)
    return;

  /* We want to make sure that we can open and write to this file. But
     the user might have asked to not delete the file, so the
     contents should not be changed. Therefore we have to open it with
     'r+'. */
  errno=0;
  tmpfile=fopen(filename, "r+");
  if (tmpfile)                        /* The file opened. */
    {
      /* Close the file. */
      errno=0;
      if(fclose(tmpfile))
        error(EXIT_FAILURE, errno, "%s", filename);

      /* See if the file should be deleted. */
      if(keep==0)
        {
          /* Make sure it is ok to delete the file. */
          if(dontdelete)
            error(EXIT_FAILURE, 0, "%s already exists and you have "
                  "asked to not remove it with the '--dontdelete' "
                  "('-D') option", filename);

          /* Delete the file: */
          errno=0;
          if(unlink(filename))
            error(EXIT_FAILURE, errno, "%s", filename);
        }
    }

  /* If the file doesn't exist, we just need to make sure if we have write
     permissions to its host directory. */
  else
    {
      /* Separate the directory part of the filename. */
      dir=gal_checkset_dir_part(filename);

      /* Make sure this directory is writable by this user. */
      errno=0;
      if( access(dir, W_OK) )
        error(EXIT_FAILURE, errno, "cannot create any file(s) in the "
              "directory '%s'", dir);

      /* Clean up. */
      free(dir);
    }
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
    error(EXIT_FAILURE, 0, "%s: a bug! The input should not be NULL. "
          "Please contact us at %s so we can see what went wrong and "
          "fix it in future updates", __func__, PACKAGE_BUGREPORT);

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
      gal_checkset_writable_remove(name, 0, dontdelete);
      return 1;
    }
  else                                 /* Not a file or a dir, ABORT */
    error(EXIT_FAILURE, 0, "%s not a file or a directory", name);

  error(EXIT_FAILURE, 0, "%s: a bug! The process should not reach the end "
        "of the function! Please contact us at %s so we can see what went "
        "wrong and fix it in future updates", __func__, PACKAGE_BUGREPORT);
  return 0;
}





/* Allocate space and write the output name (outname) based on a given
   input name (inname). The suffix of the input name (if present) will
   be removed and the given suffix will be put in the end. */
char *
gal_checkset_automatic_output(struct gal_options_common_params *cp,
                              char *inname, char *suffix)
{
  char *out;
  size_t i, l, offset=0;

  /* Merge the contents of the input name and suffix name (while also
     allocating the necessary space).*/
  out=gal_checkset_malloc_cat(inname, suffix);

  /* If there is actually a suffix, replace it with the (possibly) existing
     suffix. */
  if(suffix)
    {
      /* Start from the end of the input array*/
      l=strlen(inname);
      for(i=l-1;i!=0;--i)
        {
          /* We don't want to touch anything before a '/' (directory
             names). We are only concerned with file names here. */
          if(out[i]=='/')
            {
              /* When '/' is the last input character, then the input is
                 clearly not a filename, but a directory name. In this
                 case, adding a suffix is meaningless (a suffix belongs to
                 a filename for Gnuastro's tools). So close the string
                 after the '/' and leave the loop. However, if the '/'
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
             '.fits.gz'). Currently the only compressed formats
             (decompressed within CFITSIO) are listed in
             'gal_fits_name_is_fits' and 'gal_fits_suffix_is_fits'.*/
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

  /* If we don't want the input directory information, remove them
     here. */
  if(!cp->keepinputdir)
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
  gal_checkset_writable_remove(out, cp->keep, cp->dontdelete);

  /* Return the resulting filename. */
  return out;
}





/* Check write-ability by trying to make a temporary file. Return 0 if the
   directory is writable, and 'errno' if it isn't. We won't be using
   'facccessat' because its not available on some systems (macOS 10.9 and
   earlier, see https://github.com/conda-forge/staged-recipes/pull/9723
   ). */
static int
checkset_directory_writable(char *dirname)
{
  int file_d;
  int errnum=0;
  char *tmpname;

  /* Set the template for the temporary file (accounting for the possible
     extra '/'). */
  if(dirname[strlen(dirname)-1]=='/')
    tmpname=gal_checkset_malloc_cat(dirname, "gnuastroXXXXXX");
  else
    tmpname=gal_checkset_malloc_cat(dirname, "/gnuastroXXXXXX");

  /* Make the temporary file (with the temporary name) and open it. */
  errno=0;
  file_d=mkstemp(tmpname);

  /* See if the file was opened (and thus created). */
  if(file_d==-1) errnum=errno;
  else
    { /* The file was created, close the file. */
      errno=0;
      if( close(file_d) == -1 ) errnum=errno;
      else
        {/* The file was closed, delete it. */
          errno=0;
          if(unlink(tmpname)==-1) errnum=errno;
        }
    }

  /* Return the final value. */
  return errnum;
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
    error(EXIT_FAILURE, errno, "cannot write output in the directory "
          "%s", indir);
  /*
  errno=0;
  printf("\n\n%s\n\n", tmpname);
  if( write(file_d, buf, strlen(buf)) == -1 )
    error(EXIT_FAILURE, errno, "%s: writing to this temporary file to "
          "check the given '%s' directory", tmpname, indir);
  */
  errno=0;
  if( close(file_d) == -1 )
    error(EXIT_FAILURE, errno, "%s: Closing this temporary file to check "
          "the given '%s' directory", tmpname, indir);

  /* Delete the temporary file: */
  errno=0;
  if(unlink(tmpname)==-1)
    error(EXIT_FAILURE, errno, "%s: removing this temporary file made "
          "to check the given '%s directory'", tmpname, indir);

  /* Remove the extra characters that were added for the random name. */
  tmpname[strlen(tmpname)-14]='\0';

  free(*dirname);
  *dirname=tmpname;
}





/* If the given directory exists and is writable, then nothing is done and
   this function returns 0. If it doesn't, it will be created. If it fails
   at creating the file, or the file isn't writable it returns a non-zero
   value: the errno, which can be directly used in 'error'. */
int
gal_checkset_mkdir(char *dirname)
{
  int errnum=0;
  struct stat st={0};

  /* See if the directory exists. */
  if( stat(dirname, &st) == -1 )
    { /* The directory doesn't exist, try making it. */
      errno=0;
      if( mkdir(dirname, 0700) == -1 )
        errnum=errno;
    }
  else
    /* The directory exists, see if its writable. */
    errnum=checkset_directory_writable(dirname);

  /* Return the final 'errno'. */
  return errnum;
}
