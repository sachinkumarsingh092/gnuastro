/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <nproc.h>               /* From Gnulib.                   */
#include <timing.h>              /* Includes time.h and sys/time.h */
#include <checkset.h>
#include <commonargs.h>
#include <configfiles.h>
#include <fixedstringmacros.h>

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
readconfig(char *filename, struct imgarithparams *p)
{
  FILE *fp;
  char *tokeephdu;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
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
        {
          gal_checkset_allocate_copy(value, &tokeephdu);
          gal_linkedlist_add_to_stll(&p->hdus, tokeephdu);
        }

      else if (strcmp(name, "mask")==0)
        gal_checkset_allocate_copy_set(value, &up->maskname,
                                       &up->masknameset);

      else if (strcmp(name, "mhdu")==0)
        gal_checkset_allocate_copy_set(value, &up->mhdu, &up->mhduset);





      /* Outputs */
      else if(strcmp(name, "output")==0)
        gal_checkset_allocate_copy_set(value, &cp->output,
                                       &cp->outputset);
      else if(strcmp(name, "type")==0)
        {
          if(p->up.typeset) continue;
          gal_checkset_known_types(value, &p->outtype, filename, lineno);
          p->up.typeset=1;
        }


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
printvalues(FILE *fp, struct imgarithparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_linkedlist_stll *hdu;
  struct gal_commonparams *cp=&p->cp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image(s):\n");

  /* The order of the HDU linked list has already been corrected, so
     just print them as they were read in. */
  for(hdu=p->hdus; hdu!=NULL; hdu=hdu->next)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("hdu", hdu->v);

  if(up->masknameset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("mask", up->maskname);
  if(up->mhdu)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("mhdu", up->mhdu);

  /* Output: */
  fprintf(fp, "\n# Output:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);
  if(up->typeset)
    gal_configfiles_print_type(fp, p->outtype);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  GAL_CONFIGFILES_PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. Note that those options which are not mandatory
   must not be listed here. */
void
checkifset(struct imgarithparams *p)
{
  int intro=0;
  char comment[100];
  size_t numhdus=gal_linkedlist_num_in_stll(p->hdus);

  /* Operating modes: */
  if(p->cp.minmapsizeset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("minmapsize");

  /* Output parameters: */
  if(p->up.typeset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("type");

  /* Make sure the number of HDUs is not less than the total number of
     FITS images. If there are more HDUs than there are FITS images,
     there is no problem (since they can come from the configuration
     files). It is expected that the user can call their own desired
     number of HDUs, and not rely on the configuration files too much,
     however, if the configuration file does contain some HDUs, then
     it will be a real pain to first clean the configuration file and
     then re-run arithmetic. The best way is to simply ignore them. */
  if(numhdus<p->numfits)
    {
      sprintf(comment, "hdu (%zu FITS file(s), %zu HDUs)",
              p->numfits, numhdus);
      GAL_CONFIGFILES_REPORT_NOTSET(comment);
    }

  /* Report the (possibly) missing options. */
  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/

/* The dash of a negative number will cause problems for the users,
   so to work properly we will go over all the options/arguments and
   if any one starts with a dash and is followed by a number, then
   the dash is replaced by NEGDASHREPLACE. */
void
dashtonegchar(int argc, char *argv[])
{
  size_t i;
  for(i=0;i<argc;++i)
    if(argv[i][0]=='-' && isdigit(argv[i][1]))
      argv[i][0]=NEGDASHREPLACE;
}





/* Return the negative character back to the dash (to be read as a
   number in imgarith.c). When the token is not a FITS file name and
   since no operators or numbers begin with NEGDASHREPLACE, so if the
   token starts with NEGDASHREPLACE and its next character is a digit,
   it must be a negative number. If not, it is either an ordinary
   number or an operator.*/
void
negchartodashcountfits(struct imgarithparams *p)
{
  struct gal_linkedlist_stll *token;

  /* Initialize the numfits variable (just incase!) */
  p->numfits=0;

  /* Go through all the tokens and do the job(s). */
  for(token=p->tokens; token!=NULL; token=token->next)
    {
      if(gal_fits_name_is_fits(token->v))
        ++p->numfits;
      else if(token->v[0]==NEGDASHREPLACE && isdigit(token->v[1]) )
        token->v[0]='-';
    }
}





/* Standard sanity checks. */
void
sanitycheck(struct imgarithparams *p)
{
  struct gal_linkedlist_stll *token;

  /* Set the output file name (if any is needed). Note that since the
     lists are already reversed, the first FITS file encountered, is
     the first FITS file given by teh user. Also, notet that these
     file name operations are only necessary for the first FITS file
     in the token list. */
  for(token=p->tokens; token!=NULL; token=token->next)
    if(gal_fits_name_is_fits(token->v))
    {
      /* Set the p->up.maskname accordingly: */
      gal_fits_file_or_ext_name(token->v, p->cp.hdu,
                                p->up.masknameset, &p->up.maskname,
                                p->up.mhdu, p->up.mhduset, "mask");

      /* Set the name of the output file: */
      if(p->cp.outputset)
        gal_checkset_check_remove_file(p->cp.output, p->cp.dontdelete);
      else
        gal_checkset_automatic_output(token->v, "_arith.fits",
                                      p->cp.removedirinfo,
                                      p->cp.dontdelete, &p->cp.output);

      /* These were only necessary for the first FITS file in the
         tokens, so break out of the loop. */
      break;
    }
}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct imgarithparams *p)
{
  struct gal_commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = num_processors(NPROC_CURRENT);
  cp->removedirinfo = 1;

  p->hdus           = NULL;
  p->tokens         = NULL;
  p->up.maskname    = NULL;

  /* The hyphen of a negative number can be confused with a dash, so
     we will temporarily replace such hyphens with other
     characters. */
  dashtonegchar(argc, argv);

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");

  /* Revert the conversion of the hyphen above back to the original
     character. */
  negchartodashcountfits(p);

  /* Add the user default values and save them if asked. */
  GAL_CONFIGFILES_CHECK_SET_CONFIG;

  /* The inputs are put in a lastin-firstout (simple) linked list, so
     change them to the correct order so the order we pop a node is
     the same order that the user input a value. */
  gal_linkedlist_reverse_stll(&p->hdus);
  gal_linkedlist_reverse_stll(&p->tokens);

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
      GAL_CONFIGFILES_REPORT_PARAMETERS_SET;

  /* Do a sanity check. */
  sanitycheck(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct imgarithparams *p, struct timeval *t1)
{
  free(p->cp.output);

  /* If there are any remaining HDUs in the hdus linked list, then
     free them. */
  if(p->hdus)
    gal_linkedlist_free_stll(p->hdus, 1);

  /* Report the duration of the job */
  if(p->cp.verb)
    gal_timing_report(t1, SPACK_NAME" finished in", 0);
}
