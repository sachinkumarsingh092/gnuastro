/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

#include <nproc.h>               /* From Gnulib.                   */

#include <gnuastro/fits.h>
#include <gnuastro/timing.h>     /* Includes time.h and sys/time.h */
#include <gnuastro/checkset.h>
#include <gnuastro/txtarray.h>
#include <gnuastro/commonargs.h>
#include <gnuastro/configfiles.h>

#include "main.h"

#include "ui.h"                  /* Needs main.h                   */
#include "args.h"                /* Needs main.h, includes argp.h. */


/* Set the file names of the places where the default parameters are
   put. */
#define CONFIG_FILE SPACK CONF_POSTFIX
#define SYSCONFIG_FILE SYSCONFIG_DIR "/" CONFIG_FILE
#define USERCONFIG_FILEEND USERCONFIG_DIR CONFIG_FILE
#define CURDIRCONFIG_FILE CURDIRCONFIG_DIR CONFIG_FILE










/**************************************************************/
/**************       Options and parameters    ***************/
/**************************************************************/
void
readconfig(char *filename, struct tableparams *p)
{
  FILE *fp;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  /*struct uiparams *up=&p->up;*/
  struct gal_commonparams *cp=&p->cp;
  char key='a';        /* Not used, just a place holder. */

  /* When the file doesn't exist or can't be opened, it is ignored. It
     might be intentional, so there is no error. If a parameter is
     missing, it will be reported after all defaults are read. */
  fp=fopen(filename, "r");
  if (fp==NULL) return;


  /* Allocate some space for `line` with `len` elements so it can
     easily be freed later on. The value of `len` is arbitarary at
     this point, during the run, getline will change it along with the
     pointer to line. */
  errno=0;
  line=malloc(len*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %lu bytes in readdefaults",
          len * sizeof *line);

  /* Read the tokens in the file:  */
  while(getline(&line, &len, fp) != -1)
    {
      /* Prepare the "name" and "value" strings, also set lineno. */
      GAL_CONFIGFILES_START_READING_LINE;




      /* Inputs: */
      if(strcmp(name, "hdu")==0)
        gal_checkset_allocate_copy_set(value, &cp->hdu, &cp->outputset);



      /* Outputs */
      else if(strcmp(name, "output")==0)
        gal_checkset_allocate_copy_set(value, &cp->output, &cp->outputset);



      /* Operating modes: */
      else if (strcmp(name, "information")==0)
        {
          if(p->up.informationset) continue;
          gal_checkset_int_zero_or_one(value, &p->up.information, name,
                                       key, SPACK, filename, lineno);
          p->up.informationset=1;
        }


      /* Read options common to all programs */
      GAL_CONFIGFILES_READ_COMMONOPTIONS_FROM_CONF


      else
        error_at_line(EXIT_FAILURE, 0, filename, lineno,
                      "`%s` not recognized.\n", name);
    }

  free(line);
  fclose(fp);
}





void
printvalues(FILE *fp, struct tableparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;


  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("hdu", cp->hdu);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  if(up->informationset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "information", up->information);

  GAL_CONFIGFILES_PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct tableparams *p)
{
  /*struct uiparams *up=&p->up;*/
  struct gal_commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("hdu");


  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}



















/**************************************************************/
/************     Read and Write column info    ***************/
/**************************************************************/
/* This function will read all the table information from a FITS table HDU
   and store them in arrays for use later. It is mainly good for getting
   all the information in a FITS table HDU. This function will only go
   through the header keywords once and does not depend on the ordering of
   the keywords, so it is much more efficient than having to ask for each
   column's information separately.*/
void
readallcolinfo(fitsfile *fitsptr, size_t ncols, int **otypecode,
               char ***otform, char ***ottype, char ***otunit)
{
  size_t index;
  char *tailptr, **c, **fc;
  int i, status=0, *typecode;
  char **tform, **ttype, **tunit;
  char keyname[FLEN_KEYWORD]="XXXXXXXXXXXXX", value[FLEN_VALUE];


  /* Allocate the arrays to keep the column information. Initialize the
     arrays with a NULL pointer to make sure that they are all found in the
     end if they are necessary.*/
  errno=0; typecode=*otypecode=malloc(ncols * sizeof *typecode);
  if(typecode==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for typecode",
          ncols * sizeof *typecode);
  errno=0; tform=*otform=malloc(ncols * sizeof *tform);
  if(tform==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for tform",
          ncols * sizeof *tform);
  fc=(c=tform)+ncols; do *c++=NULL; while(c<fc);
  errno=0; ttype=*ottype=malloc(ncols * sizeof *ttype);
  if(ttype==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for ttype",
          ncols * sizeof *ttype);
  fc=(c=ttype)+ncols; do *c++=NULL; while(c<fc);
  errno=0; tunit=*otunit=malloc(ncols * sizeof *tunit);
  if(tunit==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for tunit",
          ncols * sizeof *tunit);
  fc=(c=tunit)+ncols; do *c++=NULL; while(c<fc);


  /* Read all the keywords one by one and if they match, then put them in
     the correct value. Note that we are starting from keyword 9 because
     according to the FITS standard, the first 8 keys in a FITS table are
     reserved. */
  for(i=9; strcmp(keyname, "END"); ++i)
    {
      /* Read the next keyword. */
      fits_read_keyn(fitsptr, i, keyname, value, NULL, &status);

      /* Check the type of the keyword. */
      if(strncmp(keyname, "TFORM", 5)==0)
        {
          /* Currently we can only read repeat==1 cases. When no number
             exists before the defined capital letter, it defaults to 1,
             but if a number exists (for example `5D'), then the repeat is
             5 (there are actually five values here. Note that value[0] is
             a single quote.*/
          if(isdigit(value[1] && value[1]!='1'))
             error(EXIT_FAILURE, 0, "The repeat value of column %d is "
                   "%c, currently Table can only use columns with a repeat "
                   "of 1.", i+1, value[1]);


          /* The values to TFORM are only a single character, so start the
             pointer to copy at 1 and put the string terminator at 3. */
          value[2]='\0';
          index=strtoul(&keyname[5], &tailptr, 10)-1;
          if(index<ncols)
            {
              gal_checkset_allocate_copy(&value[1], &tform[index] );
              typecode[index]=gal_fits_tform_to_dtype(value[1]);
            }
        }
      else if(strncmp(keyname, "TTYPE", 5)==0)
        {
          /* All strings in CFITSIO start and finish with single quotation
             marks, CFITSIO puts them in itsself, so if we don't remove
             them here, we might have duplicates later, its easier to just
             remove them to have a simple string that might be used else
             where too (without the single quotes). */
          value[strlen(value)-1]='\0';
          index=strtoul(&keyname[5], &tailptr, 10)-1;
          if(index<ncols)
            gal_checkset_allocate_copy(&value[1], &ttype[index] );
        }
      else if(strncmp(keyname, "TUNIT", 5)==0)
        {
          /* similar to ttype, see above.*/
          value[strlen(value)-1]='\0';
          index=strtoul(&keyname[5], &tailptr, 10)-1;
          if(index<ncols)
            gal_checkset_allocate_copy(&value[1], &tunit[index] );
        }
    }


  /* Check if the mandatory TFORMn values are set: */
  for(i=0;i<ncols;++i)
    if(!tform[i])
      error(EXIT_FAILURE, 0, "TFORM%d could not be found in header", i+1);


  /* For a checkup:
  for(i=1;i<=*ncols;++i)
    printf("%d: %s, %s, %s\n", i, tform[i], ttype[i], tunit[i]);
  */
}




/* Print the column information. */
void
printinfo(struct tableparams *p)
{
  size_t i;
  char *typestring=NULL;

  printf("%s (hdu: %s)\n", p->up.fitsname, p->cp.hdu);
  printf("Column information\n");
  printf("---------------------------------------------------------\n");
  printf("%-5s%-25s%-20s%s\n", "No.", "Column name", "Data type",
         "Units");
  printf("---------------------------------------------------------\n");
  for(i=0;i<p->ncols;++i)
    {
      switch(p->typecode[i])
        {
        case TBIT:
          typestring="bit";
          break;
        case TBYTE:
          typestring="byte";
          break;
        case TLOGICAL:
          typestring="logicals";
          break;
        case TSTRING:
          typestring="string";
          break;
        case TSHORT:
          typestring="short";
          break;
        case TLONG:
          typestring="long";
          break;
        case TLONGLONG:
          typestring="longlong";
          break;
        case TFLOAT:
          typestring="float";
          break;
        case TDOUBLE:
          typestring="double";
          break;
        case TCOMPLEX:
          typestring="complex";
          break;
        case TDBLCOMPLEX:
          typestring="dblcomplex";
          break;
        case TSBYTE:
          typestring="signed byte";
          break;
        case TUINT:
          typestring="unsigned int";
          break;
        case TUSHORT:
          typestring="unsigned short";
          break;
        default:
          error(EXIT_FAILURE, 0, "%d (from TFORM%lu='%c') is not a "
                "recognized CFITSIO datatype.",
                p->typecode[i], i, p->tform[i][0]);
        }
      printf("%-5lu%-25s%-20s%s\n", i+1, p->ttype[i] ? p->ttype[i] : "---",
             typestring, p->tunit[i] ? p->tunit[i] : "---");
    }


  /* Print the number of rows: */
  printf("---------------------------------------------------------\n");
  printf("Number of rows: %lu\n", p->nrows);
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct tableparams *p)
{

  /* Set the FITS pointer and check the type of the fits file. */
  if(p->up.fitsname)
    {
      gal_fits_read_hdu(p->up.fitsname, p->cp.hdu, 1, &p->fitsptr);
      gal_fits_table_size(p->fitsptr, &p->nrows, &p->ncols);
    }
  else
    error(EXIT_FAILURE, 0, "Table is a new addition to Gnuastro and "
          "and under heavy development, it currently doesn't support "
          "anything other than a FITS binary table.");



  /* Print the column information and exit successfully if the
     `--information' option is given. */
  if(p->up.information)
    {
      if(p->up.fitsname)
        {
          readallcolinfo(p->fitsptr, p->ncols, &p->typecode,
                         &p->tform, &p->ttype, &p->tunit);
          printinfo(p);
          freeandreport(p);
          exit(EXIT_SUCCESS);
        }
      else
        error(EXIT_FAILURE, 0, "the `--information' (`-i') option is only "
              "defined for FITS tables");
    }
}



















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
void
preparearrays(struct tableparams *p)
{

}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct tableparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = num_processors(NPROC_CURRENT);
  cp->removedirinfo = 1;

  /* Initialize this utility's special variables. */
  up->txtname=up->fitsname=NULL;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");

  /* Add the user default values and save them if asked. */
  GAL_CONFIGFILES_CHECK_SET_CONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    GAL_CONFIGFILES_REPORT_PARAMETERS_SET;

  /* Do a sanity check. */
  sanitycheck(p);

  /* Make the array of input images. */
  preparearrays(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct tableparams *p)
{
  int status=0;

  /* Free the allocated arrays: */
  free(p->tform);
  free(p->ttype);
  free(p->tunit);
  free(p->cp.hdu);
  free(p->typecode);
  free(p->cp.output);

  /* Close the FITS file: */
  if(p->up.fitsname && fits_close_file(p->fitsptr, &status))
    gal_fits_io_error(status, NULL);
}
