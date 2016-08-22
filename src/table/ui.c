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
#include <regex.h>
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
  /*struct uiparams *up=&p->up;*/
  struct gal_commonparams *cp=&p->cp;
  char *line, *name, *value, *tstring;
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

      else if(strcmp(name, "column")==0)
        {
          gal_checkset_allocate_copy(value, &tstring);
          gal_linkedlist_add_to_stll(&p->up.columns, tstring);
        }

      else if(strcmp(name, "ignorecase")==0)
        {
          if(p->up.ignorecaseset) continue;
          gal_checkset_int_zero_or_one(value, &p->up.ignorecase, "ignorecase",
                                       key, SPACK, filename, lineno);
          p->up.ignorecaseset=1;
        }



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
  struct gal_linkedlist_stll *tmp;
  struct gal_commonparams *cp=&p->cp;


  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("hdu", cp->hdu);
  if(up->columns)
    for(tmp=up->columns;tmp!=NULL;tmp=tmp->next)
      GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("column", tmp->v);
  if(up->ignorecaseset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "ignorecase", up->ignorecase);


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

  /* The user doesn't just want to see the table information, they actually
     want to print something. So if no columns are specified, then print
     all columns. */
}



















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
/* FUnction to print regular expression error. This is taken from the GNU C
   library manual, with small modifications to fit out style, */
void
regexerrorexit(int errcode, regex_t *compiled, char *input)
{
  char *regexerrbuf;
  size_t length = regerror (errcode, compiled, NULL, 0);

  errno=0;
  regexerrbuf=malloc(length);
  if(regexerrbuf==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for regexerrbuf", length);
  (void) regerror(errcode, compiled, regexerrbuf, length);

  error(EXIT_FAILURE, 0, "Regular expression error: %s in value to "
        "`--column' (`-c'): `%s'", regexerrbuf, input);
}





/* If values were given to the columns option, use them to make a list of
   columns that must be output. Note that because regular expressions are
   also allowed as values to the column option, we have no idea how many
   columns must be printed at first, so we define a linked list to keep the
   column numbers for later.*/
void
outputcolumns(struct tableparams *p)
{
  size_t i;
  long tlong;
  regex_t *regex;
  int regreturn=0;
  char *tailptr, *colstring;
  struct gal_linkedlist_sll *colsll=NULL;

  /* Go through each given column string and take the appropriate step. */
  while(p->up.columns)
    {
      /* Pop out the top node in the string linked list. */
      gal_linkedlist_pop_from_stll(&p->up.columns, &colstring);


      /* First, see if this given column is an integer or a name/regex. If
         the string is an integer, then tailptr shoult point to the null
         character. If it points to anything else, it shows that we are not
         dealing with an integer (usable as a column number). So floating
         point values are also not acceptable. */
      tlong=strtol(colstring, &tailptr,0);
      if(*tailptr=='\0')
        {
          /* Make sure we are not dealing with a negative number! */
          if(tlong<0)
            error(EXIT_FAILURE, 0, "the column numbers given to the "
                  "`--column' (`-c') option must not be negative, you "
                  "have given a value of `%ld'", tlong);

          /* Check if the given value is not larger than the number of
             columns in the input catalog. */
          if(tlong>p->ncols)
            error(EXIT_FAILURE, 0, "%s (hdu: %s) has %lu columns, but "
                  "you have asked for column number %lu", p->up.fitsname,
                  p->cp.hdu, p->ncols, tlong);

          /* Everything seems to be fine, put this column number in the
             output column numbers linked list. Note that internally, the
             column numbers start from 0, not 1.*/
          gal_linkedlist_add_to_sll(&colsll, tlong-1);
        }
      else
        {
          /* First we need to make sure that the full column information is
             ready (so we can parse the values of the column names in
             p->ttype). Note that the parsing function to read all column
             information is not set by default. Note that this is only done
             once (for the first string value to the `--column' option).*/
          if(p->ttype==NULL)
            readallcolinfo(p->fitsptr, p->ncols, &p->typecode,
                           &p->tform, &p->ttype, &p->tunit);

          /* Allocate the regex_t structure: */
          errno=0; regex=malloc(sizeof *regex);
          if(regex==NULL)
            error(EXIT_FAILURE, errno, "%lu bytes for regex", sizeof *regex);

          /* Go through all the columns names and see if this matches
             them. But first we have to "compile" the string into the
             regular expression, see the "POSIX Regular Expression
             Compilation" section of the GNU C Library.

             About the case of the string: the FITS standard says: "It is
             _strongly recommended_ that every field of the table be
             assigned a unique, case insensitive name with this keyword..."
             So the column names can be case-sensitive.

             Here, we don't care about the details of a match, the only
             important thing is a match, so we are using the REG_NOSUB
             flag.*/
          regreturn=0;
          regreturn=regcomp(regex, colstring, ( p->up.ignorecase
                                                ? REG_NOSUB + REG_ICASE
                                                : REG_NOSUB ) );
          if(regreturn)
            regexerrorexit(regreturn, regex, colstring);


          /* With the regex structure "compile"d you can go through all the
             column names. Just note that column names are not mandatory in
             the FITS standard, so some (or all) columns might not have
             names, if so `p->ttype[i]' will be NULL. */
          for(i=0;i<p->ncols;++i)
            if(p->ttype[i] && regexec(regex, p->ttype[i], 0, 0, 0)==0)
                gal_linkedlist_add_to_sll(&colsll, i);

          /* Free the regex_t structure: */
          regfree(regex);
        }

      /* We don't need this user provided column string any more. */
      free(colstring);
    }

  /* Put the desired columns (in reverse order due to the nature of a
     linked list) into an array to read from later, then pop everything to
     the un-used `i' variable from the list (which will automatically free
     any allocate space). */
  gal_linkedlist_sll_to_array(colsll, &p->ocols, &p->nocols, 1);
  while(colsll) gal_linkedlist_pop_from_sll(&colsll, &i);
}





void
preparearrays(struct tableparams *p)
{
  size_t i;

  /* Set the columns that should be included in the output. If up->columns
     is set, then use it, otherwise, set all the columns for printing. */
  if(p->up.columns)
    outputcolumns(p);
  else
    {
      p->nocols=p->ncols;
      errno=0;
      p->ocols=malloc(p->nocols * sizeof *p->ocols);
      if(p->ocols==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for p->ocols",
              p->nocols * sizeof *p->ocols);
      for(i=0;i<p->nocols;++i) p->ocols[i]=i;
    }
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

  /* Initialize this utility's pointers to NULL. */
  up->columns=NULL;
  up->txtname=up->fitsname=NULL;
  p->tform=p->ttype=p->tunit=NULL;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");

  /* Add the user default values and save them if asked. */
  GAL_CONFIGFILES_CHECK_SET_CONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Reverse the columns linked list here (before possibly printing).*/
  gal_linkedlist_reverse_stll(&up->columns);

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
  size_t i;
  int status=0;

  /* Free the allocated arrays: */
  free(p->ocols);
  free(p->cp.hdu);
  free(p->typecode);
  free(p->cp.output);

  /* Free the internal pointers first, then the actual arrays: */
  for(i=0;i<p->ncols;++i)
    {
      if(p->tform) free(p->tform[i]);
      if(p->ttype) free(p->ttype[i]);
      if(p->tunit) free(p->tunit[i]);
    }
  free(p->tform);
  free(p->ttype);
  free(p->tunit);

  /* Close the FITS file: */
  if(p->up.fitsname && fits_close_file(p->fitsptr, &status))
    gal_fits_io_error(status, NULL);
}
