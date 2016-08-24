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
/* Check the value given for the fge option. */
void
checksetfge(char *optarg, int *fge, char *filename, size_t lineno)
{
  *fge=optarg[0];
  if( *fge!='f' && *fge!='g' && *fge!='g' )
    {
      if(filename)
        error_at_line(EXIT_FAILURE, 0, filename, lineno,
                      "the value of `fge' must only be one of the three "
                      "`f', `g', or `e' characters. You have given `%s'",
                      optarg);
      else
        error(EXIT_FAILURE, 0, "the value of `--fge' (`-f') must only be "
              "one of the three `f', `g', or `e' characters. You have "
              "given `%s'", optarg);
    }
}





void
checksetfitstabletype(char *optarg, int *fitstabletype, char *filename,
                   size_t lineno)
{
  if( !strcmp(optarg, "ascii") )
    *fitstabletype=ASCII_TBL;
  else if( !strcmp(optarg, "binary") )
    *fitstabletype=BINARY_TBL;
  else
    {
      if(filename)
        error_at_line(EXIT_FAILURE, 0, filename, lineno,
                      "The value to the `fitstabletype' must be either "
                      "`ascii', or `binary'. You have given `%s'", optarg);
      else
        error(EXIT_FAILURE, 0, "The value to the `--fitstabletype' (`-t') "
              "option must be either one of `ascii', or `binary'. You have, "
              "given `%s'", optarg);
    }
}





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
        gal_checkset_allocate_copy_set(value, &cp->hdu, &cp->hduset);

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

      else if (strcmp(name, "feg")==0)
        {
          if(p->up.fegset) continue;
          checksetfge(value, &p->up.feg, filename, lineno);
          p->up.fegset=1;
        }

      else if (strcmp(name, "sintwidth")==0)
        {
          if(p->up.sintwidthset) continue;
          gal_checkset_sizet_el_zero(value, &p->up.sintwidth, name,
                                       key, SPACK, filename, lineno);
          p->up.sintwidthset=1;
        }

      else if (strcmp(name, "lintwidth")==0)
        {
          if(p->up.lintwidthset) continue;
          gal_checkset_sizet_el_zero(value, &p->up.lintwidth, name,
                                       key, SPACK, filename, lineno);
          p->up.lintwidthset=1;
        }

      else if (strcmp(name, "floatwidth")==0)
        {
          if(p->up.floatwidthset) continue;
          gal_checkset_sizet_el_zero(value, &p->up.floatwidth, name,
                                       key, SPACK, filename, lineno);
          p->up.floatwidthset=1;
        }

      else if (strcmp(name, "doublewidth")==0)
        {
          if(p->up.doublewidthset) continue;
          gal_checkset_sizet_el_zero(value, &p->up.doublewidth, name,
                                       key, SPACK, filename, lineno);
          p->up.doublewidthset=1;
        }

      else if (strcmp(name, "strwidth")==0)
        {
          if(p->up.strwidthset) continue;
          gal_checkset_sizet_el_zero(value, &p->up.strwidth, name,
                                       key, SPACK, filename, lineno);
          p->up.strwidthset=1;
        }

      else if (strcmp(name, "floatprecision")==0)
        {
          if(p->up.floatprecisionset) continue;
          gal_checkset_sizet_el_zero(value, &p->up.floatprecision,
                                       name, key, SPACK, filename, lineno);
          p->up.floatprecisionset=1;
        }

      else if (strcmp(name, "doubleprecision")==0)
        {
          if(p->up.doubleprecisionset) continue;
          gal_checkset_sizet_el_zero(value, &p->up.doubleprecision,
                                       name, key, SPACK, filename, lineno);
          p->up.doubleprecisionset=1;
        }
      else if (strcmp(name, "fitstabletype")==0)
        {
          if(p->up.fitstabletypeset) continue;
          checksetfitstabletype(value, &p->fitstabletype, filename, lineno);
          p->up.fitstabletypeset=1;
        }



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
  fprintf(fp, "\n# Input:\n");
  if(cp->hduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("hdu", cp->hdu);
  if(up->columns)
    for(tmp=up->columns;tmp!=NULL;tmp=tmp->next)
      GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("column", tmp->v);
  if(up->ignorecaseset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "ignorecase", up->ignorecase);


  fprintf(fp, "\n# Output:\n");
  if(up->fegset)
    fprintf(fp, CONF_SHOWFMT"%c\n", "feg", up->feg);
  if(up->sintwidthset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "sintwidth", up->sintwidth);
  if(up->lintwidthset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "lintwidth", up->lintwidth);
  if(up->floatwidthset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "floatwidth", up->floatwidth);
  if(up->doublewidthset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "doublewidth", up->doublewidth);
  if(up->strwidthset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "strwidth", up->strwidth);
  if(up->floatprecisionset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "floatprecision", up->floatprecision);
  if(up->doubleprecisionset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "doubleprecision", up->doubleprecision);
  if(up->fitstabletypeset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "fitstabletype",
            p->fitstabletype==ASCII_TBL ? "ascii" : "binary");


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
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("hdu");

  if(up->fegset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("feg");
  if(up->sintwidthset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("sintwidth");
  if(up->lintwidthset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("lintwidth");
  if(up->floatwidthset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("floatwidth");
  if(up->doublewidthset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("doublewidth");
  if(up->strwidthset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("strwidth");
  if(up->floatprecisionset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("floatprecision");
  if(up->doubleprecisionset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("doubleprecision");
  if(up->fitstabletypeset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("fitstabletype");

  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/************     Read and Write column info    ***************/
/**************************************************************/

/*  */
  /* Allocate the arrays to keep the column information. Initialize the
     arrays with a NULL pointer to make sure that they are all found in the
     end if they are necessary. The values that are directly read from the
     input file are initialized to NULL, so they can be treated
     appropriately if they do not exist in the input (for example if they
     are mandatory, or if they need to be printed).*/
void
allocinputcolinfo(struct tableparams *p)
{
  char **c, **fc;
  size_t ncols=p->up.ncols;
  struct uiparams *up=&p->up;

  /* up->datatype keeps the type of data in each column as CFITSIO type
     macros. */
  errno=0;
  up->datatype=malloc(ncols * sizeof *up->datatype);
  if(up->datatype==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for up->datatype",
          ncols * sizeof *up->datatype);

  /* up->ttstr keeps the actual string used to specify the datatype. */
  errno=0;
  up->ttstr=malloc(ncols * sizeof *up->ttstr);
  if(up->ttstr==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for ttstr",
          ncols * sizeof *up->ttstr);
  fc=(c=up->ttstr)+ncols; do *c++=NULL; while(c<fc);

  /* up->tname keeps the name of the column. */
  errno=0;
  up->tname=malloc(ncols * sizeof *up->tname);
  if(up->tname==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for tname",
          ncols * sizeof *up->tname);
  fc=(c=up->tname)+ncols; do *c++=NULL; while(c<fc);

  /* up->tunit keeps the input units of the column. */
  errno=0;
  up->tunit=malloc(ncols * sizeof *up->tunit);
  if(up->tunit==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for tunit",
          ncols * sizeof *up->tunit);
  fc=(c=up->tunit)+ncols; do *c++=NULL; while(c<fc);
}





/* This function will read all the table information from a FITS table HDU
   and store them in arrays for use later. It is mainly good for getting
   all the information in a FITS table HDU. This function will only go
   through the header keywords once and does not depend on the ordering of
   the keywords, so it is much more efficient than having to ask for each
   column's information separately.*/
void
allfitscolinfo(struct tableparams *p)
{
  char *tailptr;
  int i, status=0;
  struct uiparams *up=&p->up;
  size_t index, ncols=p->up.ncols;
  char keyname[FLEN_KEYWORD]="XXXXXXXXXXXXX", value[FLEN_VALUE];


  /* First allocate the arrays to keep the input column information. */
  allocinputcolinfo(p);


  /* Read all the keywords one by one and if they match, then put them in
     the correct value. Note that we are starting from keyword 9 because
     according to the FITS standard, the first 8 keys in a FITS table are
     reserved. */
  for(i=9; strcmp(keyname, "END"); ++i)
    {
      /* Read the next keyword. */
      fits_read_keyn(p->fitsptr, i, keyname, value, NULL, &status);

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
              gal_checkset_allocate_copy(&value[1], &up->ttstr[index] );
              up->datatype[index]=gal_fits_tform_to_dtype(value[1]);
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
            gal_checkset_allocate_copy(&value[1], &up->tname[index] );
        }
      else if(strncmp(keyname, "TUNIT", 5)==0)
        {
          /* similar to tname, see above.*/
          value[strlen(value)-1]='\0';
          index=strtoul(&keyname[5], &tailptr, 10)-1;
          if(index<ncols)
            gal_checkset_allocate_copy(&value[1], &up->tunit[index] );
        }
    }


  /* Check if the mandatory TFORMn values are set: */
  for(i=0;i<ncols;++i)
    if(!up->ttstr[i])
      error(EXIT_FAILURE, 0, "TFORM%d could not be found in header", i+1);
}





/* Prepare column information from a text input file. Note that we are
   currently using Gnuastro's very simple txtarray library, which was only
   designed for a 2D array of floating point numbers. Later, we must update
   it to be more aware of the types of input columns and also accept
   non-number columns.*/
void
alltxtcolinfo(struct tableparams *p)
{
  size_t i;
  size_t ncols=p->up.ncols;

  /* Check if there were any strings in the array. If there were strings,
     then warn the user that we currently can't deal with them. */
  if(gal_checkset_check_file_report(GAL_TXTARRAY_LOG))
    error(EXIT_FAILURE, 0, "The input text file `%s' contained "
          "non-numerical values (probably strings). Please see `%s' for "
          "a listing of all such terms. Currently Table cannot operate on "
          "such files. We are working on correcting this issue.",
          p->up.txtname, GAL_TXTARRAY_LOG);

  /* Allocate the arrays to keep the input column information. */
  allocinputcolinfo(p);

  /* Set all the column types to double and leave the other fields
     blank.*/
  for(i=0;i<ncols;++i)
    {
      p->up.datatype[i]=TDOUBLE;
      gal_checkset_allocate_copy("D", &p->up.ttstr[i] );
      gal_checkset_allocate_copy("", &p->up.tname[i] );
      gal_checkset_allocate_copy("", &p->up.tunit[i] );
    }
}




/* Print the column information. */
void
printinfo(struct tableparams *p)
{
  size_t i;
  char *typestring=NULL;
  struct uiparams *up=&p->up;

  printf("---------------------------------------------------------\n");
  if(up->fitsname)
    printf("%s (hdu: %s)\n", p->up.fitsname, p->cp.hdu);
  else
    printf("%s\n", p->up.txtname);
  printf("%-5s%-25s%-15s%s\n", "No.", "Column name", "Units",
         "Data type");
  printf("---------------------------------------------------------\n");
  for(i=0;i<up->ncols;++i)
    {
      switch(up->datatype[i])
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
                "recognized CFITSIO datatype.", up->datatype[i],
                i, up->ttstr[i][0]);
        }
      printf("%-5lu%-25s%-15s%s\n", i+1, up->tname[i] ? up->tname[i] : "---",
             up->tunit[i] ? up->tunit[i] : "---", typestring);
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
  struct uiparams *up=&p->up;

  /* If the desired FITS output type is ASCII, then inform the user that
     this type is not yet supported. */
  if(up->fitsname && p->fitstabletype==ASCII_TBL)
    error(EXIT_FAILURE, 0, "output to ASCII type FITS tables is currently "
          "not supported (due to lack of need!). If you need this feature, "
          "please get in touch with us at %s so we can increase the "
          "priority of this feature.", PACKAGE_BUGREPORT);

  if(up->fitsname)
    {
      /* Set the FITS pointer and check the type of the fits file. */
      gal_fits_read_hdu(p->up.fitsname, p->cp.hdu, 1, &p->fitsptr);
      gal_fits_table_size(p->fitsptr, &p->nrows, &up->ncols);
      up->infitstabletype=gal_fits_table_type(p->fitsptr);
      allfitscolinfo(p);
    }
  else
    {
      /* Read the text file into an input array and make the basic column
         information. */
      gal_txtarray_txt_to_array(up->txtname, &up->txtarray, &p->nrows,
                                &up->ncols);
      alltxtcolinfo(p);
    }

  /* Print the column information and exit successfully if the
     `--information' option is given. */
  if(p->up.information)
    {
      printinfo(p);
      freeandreport(p);
      exit(EXIT_SUCCESS);
    }

  /* Check the status of the output file. If no output file is set, then
     the output will be printed to standard output on the terminal.*/
  p->outputtofits=p->outputtotxt=p->outputtostdout=0;
  if(p->cp.outputset)
    {
      /* First check if the file exists and remove it if it does. */
      gal_checkset_check_remove_file(p->cp.output, p->cp.dontdelete);

      /* Now set the type of output. */
      if(gal_fits_name_is_fits(p->cp.output))
        p->outputtofits=1;
      else
        p->outputtotxt=1;
    }
  else
    p->outputtostdout=1;
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
  long tlong;
  regex_t *regex;
  int regreturn=0;
  size_t i, inindex;
  char *tailptr, *colstring;
  struct uiparams *up=&p->up;
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
          if(tlong>up->ncols)
            error(EXIT_FAILURE, 0, "%s (hdu: %s) has %lu columns, but "
                  "you have asked for column number %lu", p->up.fitsname,
                  p->cp.hdu, up->ncols, tlong);

          /* Everything seems to be fine, put this column number in the
             output column numbers linked list. Note that internally, the
             column numbers start from 0, not 1.*/
          gal_linkedlist_add_to_sll(&colsll, tlong-1);
        }
      else
        {
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
             names, if so `p->tname[i]' will be NULL. */
          for(i=0;i<up->ncols;++i)
            if(up->tname[i] && regexec(regex, up->tname[i], 0, 0, 0)==0)
              gal_linkedlist_add_to_sll(&colsll, i);

          /* Free the regex_t structure: */
          regfree(regex);
        }

      /* We don't need this user provided column string any more. */
      free(colstring);
    }

  /* Based on the number of columns found above, allocate an array of
     `outcolumn' structures to keep the information for each column. */
  p->nocols=gal_linkedlist_num_in_sll(colsll);
  errno=0;
  p->ocols=malloc(p->nocols*sizeof *p->ocols);
  if(p->ocols==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->ocols",
          p->nocols*sizeof *p->ocols);

  /* Fill in the output column with the needed input table
     information. Note that a simple linked list is first-in-last-out, so
     we have to fill in the output columns in reverse order. Also, note
     that we are popping from the linked list keeping the indexs of the
     output columns and thus also freeing their allocated space. */
  i=p->nocols-1;
  while(colsll)
    {
      gal_linkedlist_pop_from_sll(&colsll, &inindex);
      p->ocols[i].datatype=up->datatype[inindex];
      p->ocols[i].inindex=inindex;
      --i;
    }
}





void
preparearrays(struct tableparams *p)
{
  size_t i;
  struct uiparams *up=&p->up;

  /* Reverse the columns linked list here (before possibly printing).*/
  gal_linkedlist_reverse_stll(&up->columns);

  /* Set the columns that should be included in the output. If up->columns
     is set, then use it, otherwise, set all the columns for printing. */
  if(p->up.columns)
    outputcolumns(p);
  else
    {
      p->nocols=up->ncols;
      errno=0;
      p->ocols=malloc(p->nocols * sizeof *p->ocols);
      if(p->ocols==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for p->ocols",
              p->nocols * sizeof *p->ocols);
      for(i=0;i<p->nocols;++i)
        {
          p->ocols[i].datatype=up->datatype[i];
          p->ocols[i].inindex=i;
        }
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
  p->ocols=NULL;
  up->columns=NULL;
  up->txtname=up->fitsname=NULL;
  up->ttstr=up->tname=up->tunit=NULL;

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
  size_t i, j;
  int status=0;
  char **rowofstrings;
  struct uiparams *up=&p->up;

  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(up->datatype);
  free(p->cp.output);

  /* Free the input column information: */
  for(i=0;i<up->ncols;++i)
    {
      if(up->ttstr) free(up->ttstr[i]);
      if(up->tname) free(up->tname[i]);
      if(up->tunit) free(up->tunit[i]);
    }
  free(up->ttstr);
  free(up->tname);
  free(up->tunit);

  /* Free the output column information: */
  for(i=0;i<p->nocols;++i)
    {
      free(p->ocols[i].nulval);
      if(p->ocols[i].datatype==TSTRING)
        {
          rowofstrings=(char **)(p->ocols[i].data);
          for(j=0;j<p->nrows;++j)
            free(rowofstrings[j]);
        }
      else
        free(p->ocols[i].data);
    }
  free(p->ocols);

  /* Close the FITS file: */
  if(p->up.fitsname && fits_close_file(p->fitsptr, &status))
    gal_fits_io_error(status, NULL);
}
