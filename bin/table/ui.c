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

#include <gnuastro/fits.h>
#include <gnuastro/table.h>
#include <gnuastro/txtarray.h>

#include <nproc.h>               /* From Gnulib.                   */
#include <timing.h>              /* Includes time.h and sys/time.h */
#include <checkset.h>
#include <commonargs.h>
#include <configfiles.h>

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
    error(EXIT_FAILURE, errno, "ui.c: %zu bytes in readdefaults",
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
          gal_linkedlist_add_to_stll(&p->columns, tstring);
        }

      else if(strcmp(name, "searchin")==0)
        {
          if(p->up.searchinset) continue;
          gal_checkset_allocate_copy(value, &p->up.searchin);
          p->up.searchinset=1;
        }

      else if(strcmp(name, "ignorecase")==0)
        {
          if(p->up.ignorecaseset) continue;
          gal_checkset_int_zero_or_one(value, &p->ignorecase, "ignorecase",
                                       key, SPACK, filename, lineno);
          p->up.ignorecaseset=1;
        }



      /* Outputs */
      else if(strcmp(name, "output")==0)
        gal_checkset_allocate_copy_set(value, &cp->output, &cp->outputset);

      else if (strcmp(name, "fitstabletype")==0)
        {
          if(p->up.fitstabletypeset) continue;
          gal_checkset_allocate_copy_set(value, &p->up.fitstabletype,
                                         &p->up.fitstabletypeset);
        }



      /* Operating modes: */
      else if (strcmp(name, "information")==0)
        {
          if(p->up.informationset) continue;
          gal_checkset_int_zero_or_one(value, &p->information, name,
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
  if(p->columns)
    for(tmp=p->columns;tmp!=NULL;tmp=tmp->next)
      GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("column", tmp->v);
  if(up->ignorecaseset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "ignorecase", p->ignorecase);
  if(up->searchinset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "searchin", up->searchin);


  fprintf(fp, "\n# Output:\n");
  if(up->fitstabletypeset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "fitstabletype", p->up.fitstabletype);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  if(up->informationset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "information", p->information);

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

  if(up->fitstabletypeset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("fitstabletype");
  if(up->searchinset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("searchin");

  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct tableparams *p)
{
  struct uiparams *up=&p->up;

  /* If we had a FITS filename, the output table type is automatically set
     with the `--fitstabletype' option. However, when the file is a text
     file, that option should be ignored and outtype should be set to
     the text file macro.*/
  if(gal_fits_name_is_fits(up->filename))
    {
      if( !strcmp(up->fitstabletype, "ascii") )
        p->outtype=GAL_TABLE_TYPE_AFITS;
      else if( !strcmp(up->fitstabletype, "binary") )
        p->outtype=GAL_TABLE_TYPE_BFITS;
      else
        error(EXIT_FAILURE, 0, "the value to the `--fitstabletype' "
              "option on the command line or `fitstabletype' variable in "
              "any of the configuration files must be either `ascii', or "
              "`binary'. You have given `%s'", up->fitstabletype);
    }
  else
    p->outtype=GAL_TABLE_TYPE_TXT;

  /* Set the searchin integer value. */
  p->searchin=gal_table_searchin_from_str(p->up.searchin);
}




















/**************************************************************/
/***************          Information         *****************/
/**************************************************************/
void
print_information_exit(struct tableparams *p)
{
  int tabletype;
  size_t i, numcols;
  gal_data_t *allcols;
  char *name, *unit, *comment;

  allcols=gal_table_info(p->up.filename, p->cp.hdu, &numcols,
                         &tabletype);

  /* Print the legend */
  printf("Column information for\n");
  printf("%s (hdu: %s)\n", p->up.filename, p->cp.hdu);
  printf("%-8s%-25s%-20s%-18s%s\n", "No.", "Name", "Units", "Type",
         "Comment");
  printf("%-8s%-25s%-20s%-18s%s\n", "---", "----", "-----", "----",
         "-------");

  /* For each column, print the information, then free them. */
  for(i=0;i<numcols;++i)
    {
      name    = allcols[i].name;       /* Just defined for easier     */
      unit    = allcols[i].unit;       /* readability. The compiiler  */
      comment = allcols[i].comment;    /* optimizer will remove them. */
      printf("%-8zu%-25s%-20s%-18s%s\n", i+1,
             name ? name : "N/A" ,
             unit ? unit : "N/A" ,
             gal_data_type_string(allcols[i].type),
             comment ? comment : "N/A");
      if(name)    free(name);
      if(unit)    free(unit);
      if(comment) free(comment);
    }

  /* Clean everything else up and return successfully. */
  free(allcols);
  freeandreport(p);
  exit(EXIT_SUCCESS);
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
void
preparearrays(struct tableparams *p)
{
  /* Pull out the columns that we are looking for. */
  gal_linkedlist_reverse_stll(&p->columns);
  gal_table_read_cols(p->up.filename, p->cp.hdu, p->columns,
                      p->searchin, p->ignorecase);


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
  p->columns=NULL;
  up->filename=NULL;

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

  /* If the user just wanted the information, just print them and exit. */
  if(p->information)
    print_information_exit(p);

  /* Make the array of input images. */
  preparearrays(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct tableparams *p)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
}
