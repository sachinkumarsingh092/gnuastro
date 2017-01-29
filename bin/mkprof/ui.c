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

#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>

#include <gnuastro/fits.h>
#include <gnuastro/table.h>
#include <gnuastro/linkedlist.h>

#include <timing.h>
#include <options.h>
#include <checkset.h>
#include <fixedstringmacros.h>

#include "main.h"

#include "ui.h"
#include "authors-cite.h"





/**************************************************************/
/*********      Argp necessary global entities     ************/
/**************************************************************/
/* Definition parameters for the argp: */
const char *
argp_program_version = PROGRAM_STRING "\n"
                       GAL_STRINGS_COPYRIGHT
                       "\n\nWritten/developed by "PROGRAM_AUTHORS;

const char *
argp_program_bug_address = PACKAGE_BUGREPORT;

static char
args_doc[] = "[BackgroundImage] Catalog";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will create a FITS image "
  "containing any number of mock astronomical profiles based on an input "
  "catalog. All the profiles will be built from the center outwards. First "
  "by Monte Carlo integration, then using the central pixel position. The "
  "tolerance level specifies when to switch to a latter.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




/* Option groups particular to this program. */
enum program_args_groups
{
  ARGS_GROUP_PROFILES = GAL_OPTIONS_GROUP_AFTER_COMMON,
  ARGS_GROUP_CATALOG,
  ARGS_GROUP_WCS,
};





/* Keys for each option.

   Available letters (-V which is used by GNU is also removed):

   a d f g j l u v
   A E G H I J L M O Q U W Z     */
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_BACKGROUND      = 'k',
  ARGS_OPTION_KEY_BACKHDU         = 'B',
  ARGS_OPTION_KEY_NAXIS1          = 'x',
  ARGS_OPTION_KEY_NAXIS2          = 'y',
  ARGS_OPTION_KEY_INPUTASCANVAS   = 'C',
  ARGS_OPTION_KEY_OVERSAMPLE      = 's',
  ARGS_OPTION_KEY_INDIVIDUAL      = 'i',
  ARGS_OPTION_KEY_NOMERGED        = 'm',
  ARGS_OPTION_KEY_TYPE            = 'T',
  ARGS_OPTION_KEY_NUMRANDOM       = 'r',
  ARGS_OPTION_KEY_TOLERANCE       = 't',
  ARGS_OPTION_KEY_TUNITINP        = 'p',
  ARGS_OPTION_KEY_XSHIFT          = 'X',
  ARGS_OPTION_KEY_YSHIFT          = 'Y',
  ARGS_OPTION_KEY_PREPFORCONV     = 'c',
  ARGS_OPTION_KEY_ZEROPOINT       = 'z',
  ARGS_OPTION_KEY_CIRCUMWIDTH     = 'w',
  ARGS_OPTION_KEY_REPLACE         = 'R',
  ARGS_OPTION_KEY_ENVSEED         = 'e',
  ARGS_OPTION_KEY_MFORFLATPIX     = 'F',

  /* Only with long version. */
  ARGS_OPTION_KEY_PSFINIMG        = 1000,
  ARGS_OPTION_KEY_MAGATPEAK,
  ARGS_OPTION_KEY_XCOL,
  ARGS_OPTION_KEY_YCOL,
  ARGS_OPTION_KEY_RACOL,
  ARGS_OPTION_KEY_DECCOL,
  ARGS_OPTION_KEY_FCOL,
  ARGS_OPTION_KEY_RCOL,
  ARGS_OPTION_KEY_NCOL,
  ARGS_OPTION_KEY_PCOL,
  ARGS_OPTION_KEY_QCOL,
  ARGS_OPTION_KEY_MCOL,
  ARGS_OPTION_KEY_TCOL,
  ARGS_OPTION_KEY_CRPIX1,
  ARGS_OPTION_KEY_CRPIX2,
  ARGS_OPTION_KEY_CRVAL1,
  ARGS_OPTION_KEY_CRVAL2,
  ARGS_OPTION_KEY_RESOLUTION,
};




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct mkprofparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;

  /* Set the necessary common parameters structure. */
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->coptions           = gal_commonopts_options;


  /* Modify the common options for this program. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    switch(cp->coptions[i].key)
      {
      case GAL_OPTIONS_KEY_HDU:
        cp->coptions[i].doc="Input catalog HDU name or number (if FITS).";
        break;

      case GAL_OPTIONS_KEY_TABLEFORMAT:
        cp->coptions[i].flags=OPTION_HIDDEN;
        break;

      case GAL_OPTIONS_KEY_SEARCHIN:
        cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
        break;

      case GAL_OPTIONS_KEY_MINMAPSIZE:
        cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
        break;
      }


  /* Read the number of threads available to the user, this should be done
     before reading command-line and configuration file options, since they
     can change it.  */
  gal_options_initialize_numthreads(cp);


  /* Set the non-zero initial values, the structure was initialized to have
     a zero/NULL value for all elements. So, this is necessary only when
     `0' is meaningful in the context of the variable. */
  p->type=GAL_DATA_TYPE_INVALID;
  p->crpix[0]=p->crpix[1]=p->crval[0]=p->crval[1]=NAN;
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct mkprofparams *p = state->input;

  /* Pass `gal_options_common_params' into the child parser.  */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then `arg`
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign (`=`). For short "
               "options, `=` should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  /* Set the key to this option. */
  switch(key)
    {

    /* Read the non-option tokens (arguments): */
    case ARGP_KEY_ARG:
      p->catname=arg;
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in `ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct mkprofparams *p)
{
  /* When a no-merged image is to be created, type is necessary. */
  if( p->type==GAL_DATA_TYPE_INVALID && p->nomerged==0)
    error(EXIT_FAILURE, 0, "an output type `--type' is necessary when a "
          "merged image is to be built.");

  /* Check if one of the coordinate columns has been given, the other is
     also given. To simplify the job, we use the fact that conditions in C
     return either a 0 (when failed) and 1 (when successful). Note that if
     neighter coordinates are specified there is no problem, the user might
     have input the other coordinate standard. We'll also check for that
     after this.*/
  if( ((p->xcol==NULL) + (p->ycol==NULL)) == 1 )
    error(EXIT_FAILURE, 0, "only `%s' has been given, please also specify "
          "a column for the position along the %s axis with the `%s' option",
          p->xcol?"xcol":"ycol", p->xcol?"Y":"X", p->xcol?"ycol":"xcol");

  if( ((p->racol==NULL) + (p->deccol==NULL)) == 1 )
    error(EXIT_FAILURE, 0, "only `%s' has been given, please also specify "
          "a column for the position along the %s axis with the `%s' option",
          p->racol?"racol":"deccol", p->racol?"Dec":"RA",
          p->xcol?"deccol":"racol");

}





/* Sanity check on options AND arguments. If only option values are to be
   checked, use `ui_read_check_only_options'. */
static void
ui_check_options_and_arguments(struct mkprofparams *p)
{
  /* Make sure an input table is given, and if it is FITS, that the HDU is
     also provided. */
  if(p->catname)
    {
      if( gal_fits_name_is_fits(p->catname) && p->cp.hdu==NULL)
        error(EXIT_FAILURE, 0, "no `hdu' specified for the input FITS table "
              "'%s', to ", p->catname);
    }
  else
    {
      error(EXIT_FAILURE, 0, "no input catalog provided. To build profiles, "
            "you need to give a catalog/table containing the information of "
            "the profiles");
    }
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_preparations(struct mkprofparams *p)
{
  struct gal_linkedlist_stll *cols=NULL;
  char *ax1col=p->racol?p->racol:p->xcol;
  char *ax2col=p->deccol?p->deccol:p->ycol;

  /* Correct/set based on the given oversampling. */
  p->naxes[0] *= p->oversample;
  p->naxes[1] *= p->oversample;
  p->halfpixel = 0.5f/p->oversample;

  /* Read the columns.
  gal_linkedlist_add_to_stll(cols, ax1col, 0);
  gal_table_read(char *filename, char *hdu, struct gal_linkedlist_stll *cols,
                 int searchin, int ignorecase, int minmapsize);
  */
}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
ui_read_check_inputs_setup(int argc, char *argv[], struct mkprofparams *p)
{
  struct gal_options_common_params *cp=&p->cp;


  /* Include the parameters necessary for argp from this program (`args.h')
     and for the common options to all Gnuastro (`commonopts.h'). We want
     to directly put the pointers to the fields in `p' and `cp', so we are
     simply including the header here to not have to use long macros in
     those headers which make them hard to read and modify. This also helps
     in having a clean environment: everything in those headers is only
     available within the scope of this function. */
#include <commonopts.h>
#include "args.h"


  /* Initialize the options and necessary information.  */
  ui_initialize_options(p, program_options, gal_commonopts_options);


  /* Read the command-line options and arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");


  /* Read the configuration files. */
  gal_options_read_config_set(&p->cp);


  /* Read the options into the program's structure, and check them and
     their relations prior to printing. */
  ui_read_check_only_options(p);


  /* Print the option values if asked. Note that this needs to be done
     after the sanity check so un-sane values are not printed in the output
     state. */
  gal_options_print_state(&p->cp);


  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments(p);


  /* Read/allocate all the necessary starting arrays. */
  ui_preparations(p);

  printf("\n... End of ui.c ...\n");
  exit(0);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct mkprofparams *p, struct timeval *t1)
{
  int status;

  /* Free all the allocated arrays. */
  free(p->cat);
  free(p->cp.hdu);
  free(p->outdir);
  free(p->basename);
  if(p->individual==0) free(p->log);

  /* p->cp.output might be equal to p->mergedimgname. In this case, if
     we simply free them after each other, there will be a double free
     error. So after freeing output, we set it to NULL since
     free(NULL) is ok.*/
  if(p->cp.output==p->mergedimgname)
    free(p->cp.output);
  else
    {
      free(p->cp.output);
      free(p->mergedimgname);
    }

  /* Free the WCS headers string that was defined for individual mode. */
  if(p->individual)
    free(p->wcsheader);

  /* Free the WCS structure. */
  if( (status=wcsvfree(&p->nwcs, &p->wcs)) )
    error(EXIT_FAILURE, 0, "wcsfree error %d: %s", status,
          wcs_errmsg[status]);

  /* Free the random number generator: */
  gsl_rng_free(p->rng);

  /* Report the duration of the job */
  if(!p->cp.quiet)
    gal_timing_report(t1,  PROGRAM_NAME" finished in", 0);
}
