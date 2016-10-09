/*********************************************************************
Header - View and manipulate a data file header
Header is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

#include <gnuastro/fits.h>
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
readconfig(char *filename, struct headerparams *p)
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



      /* Outputs */




      /* Operating modes: */
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
printvalues(FILE *fp, struct headerparams *p)
{
  /*struct uiparams *up=&p->up;*/
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
  GAL_CONFIGFILES_PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct headerparams *p)
{
  /*struct uiparams *up=&p->up;*/
  struct gal_commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("hdu");


  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct headerparams *p)
{
  if(p->delete || p->up.update || p->up.write || p->asis || p->comment
     || p->history || p->date || p->up.rename)
    p->onlyview=0;
  else
    p->onlyview=1;
}



















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
void
setuprename(struct headerparams *p)
{
  char *c;
  struct gal_linkedlist_stll *tmp;

  for(tmp=p->up.rename; tmp!=NULL; tmp=tmp->next)
    {
      /* `c' is created in case of an error, so the input value can be
         reported. */
      errno=0;
      c=malloc(strlen(tmp->v) + 1);
      if(c==NULL) error(EXIT_FAILURE, errno, "space for c in setuprename");
      strcpy(c, tmp->v);

      /* Tokenize the input. */
      gal_linkedlist_add_to_stll(&p->renamefrom, strtok(tmp->v, ", "));
      gal_linkedlist_add_to_stll(&p->renameto, strtok(NULL, ", "));
      if(p->renamefrom->v==NULL || p->renameto->v==NULL)
        error(EXIT_FAILURE, 0, "`%s' could not be tokenized in order to "
              "complete rename. There should be a space character "
              "or a comma (,) between the two keyword names. If you have "
              "used the space character, be sure to enclose the value to "
              "the `--rename' option in double quotation marks", c);
      free(c);
    }
  /*
  {
    struct gal_linkedlist_stll *tmp2=p->renameto;
    for(tmp=p->renamefrom; tmp!=NULL; tmp=tmp->next)
      {
        printf("%s to %s\n", tmp->v, tmp2->v);
        tmp2=tmp2->next;
      }
  }
  */
}





void
fillfitsheaderll(struct gal_linkedlist_stll *input,
                 struct gal_fits_key_ll **output)
{
  long l, *lp;
  void *fvalue;
  double d, *dp;
  struct gal_linkedlist_stll *tmp;
  int i=0, datatype, vfree;
  char *c, *cf, *start, *tailptr;
  char *original, *keyname, *value, *comment, *unit;

  for(tmp=input; tmp!=NULL; tmp=tmp->next)
    {
      i=0;
      tailptr=NULL;

      /* `c' is created in case of an error, so the input value can be
         reported. */
      errno=0;
      original=malloc(strlen(tmp->v)+1);
      if(original==NULL)
        error(EXIT_FAILURE, errno, "space for c in setuprename");
      strcpy(original, tmp->v);

      /* Tokenize the input. Note that strlen does not include the \0
         character. So we have added it with a 1. */
      cf=(c=start=tmp->v)+strlen(tmp->v)+1;
      keyname=value=comment=unit=NULL;
      do
        {
          switch(*c)
            {
            case ',': case '\0':
              *c='\0';
              if(start!=c)
                switch(i)
                  {
                  case 0:
                    keyname=start;
                    break;
                  case 1:
                    value=start;
                    break;
                  case 2:
                    comment=start;
                    break;
                  case 3:
                    unit=start;
                    break;
                  default:
                    error(EXIT_FAILURE, 0, "%s: only three commas should "
                          "be given in the write or update keyword "
                          "options. The general expected format is:\n"
                          "    KEYWORD,value,\"a comment string\",unit\n",
                          original);
                  }
              ++i;
              start=c+1;
              break;

            default:
              break;
            }
        }
      while(++c<cf);
      if(keyname==NULL)
        error(EXIT_FAILURE, 0, "the keyword in %s was not readable. "
              "The general expected format is:\n"
              "    KEYWORD,value,\"a comment string\",unit\n"
              "Any space characters around the the comma (,) characters "
              "will be seen as part of the respective token", original);
      /*
      printf("\n\n-%s-\n-%s-\n-%s-\n-%s-\n", keyname, value, comment, unit);
      */

      /* Find the datatype of the value: */
      errno=0;
      l=strtol(value, &tailptr, 10);
      if(*tailptr=='\0' && errno==0)
        {
          vfree=1;
          datatype=TLONG;
          errno=0;
          fvalue=lp=malloc(sizeof *lp);
          if(lp==NULL)
            error(EXIT_FAILURE, errno, "%zu bytes for long integer",
                  sizeof *lp);
          *lp=l;
        }
      else
        {
          errno=0;
          d=strtod(value, &tailptr);
          if(*tailptr=='\0' && errno==0)
            {
              vfree=1;
              datatype=TDOUBLE;
              errno=0;
              fvalue=dp=malloc(sizeof *dp);
              if(dp==NULL)
                error(EXIT_FAILURE, errno, "%zu bytes for double",
                      sizeof *dp);
              *dp=d;
            }
          else
            { fvalue=value; datatype=TSTRING; vfree=0; }
        }


      gal_fits_add_to_key_ll(output, datatype, keyname, 0,
                             fvalue, vfree, comment, 0, unit);
      free(original);
    }
}





void
preparearrays(struct headerparams *p)
{
  size_t len;
  char *ffname;
  int status=0, iomode;

  /* Add hdu to filename: */
  errno=0;
  len=strlen(p->up.inputname)+strlen(p->cp.hdu)+4;
  ffname=malloc(len*sizeof *ffname);
  if(ffname==NULL)
    error(EXIT_FAILURE, errno, "%zu characters", len);
  sprintf(ffname, "%s[%s#]", p->up.inputname, p->cp.hdu);

  /* Open the FITS file: */
  if(p->onlyview)
    iomode=READONLY;
  else
    iomode=READWRITE;
  if( fits_open_file(&p->fptr, ffname, iomode, &status) )
    gal_fits_io_error(status, "reading file");
  free(ffname);

  /* Separate the comma-separated values:  */
  if(p->up.rename)
    setuprename(p);
  if(p->up.update)
    fillfitsheaderll(p->up.update, &p->update);
  if(p->up.write)
    fillfitsheaderll(p->up.write, &p->write);

}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct headerparams *p)
{
  struct gal_commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = num_processors(NPROC_CURRENT);
  cp->removedirinfo = 1;

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
freeandreport(struct headerparams *p)
{
  int status=0;

  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);

  /* Close the FITS file: */
  if(fits_close_file(p->fptr, &status))
    gal_fits_io_error(status, NULL);

  if(p->wcs)
    wcsvfree(&p->nwcs, &p->wcs);
}
