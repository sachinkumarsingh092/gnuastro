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

#include <nproc.h>  /* from Gnulib, in Gnuastro's source */
#include <timing.h>
#include <options.h>
#include <checkset.h>

#include "main.h"

#include "ui.h"
#include "args.h"









/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* ONLY check options. When arguments are involved, do the check in
   `ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct mkprofparams *p)
{
  size_t i;
  struct argp_option *o=options;
  struct gal_linkedlist_stll *namell=NULL, *docll=NULL;

  /* Put the program's option values into the structure. */
  for(i=0; !gal_options_is_last(&o[i]); ++i)
    if( o[i].key && o[i].name )
      switch(o[i].key)
        {

          /* Input */
        case ARGS_OPTION_BACKHDU_KEY:
          gal_checkset_allocate_copy(o[i].value, &p->backhdu);
          break;


          /* Output */
        case ARGS_OPTION_NAXIS1_KEY:
          printf("\nJust before checking value to naxis1.\n"); exit(0);
          gal_options_check_set(&o[i], &p->naxes[0], GAL_OPTIONS_RANGE_GT_0);
          break;

        case ARGS_OPTION_NAXIS2_KEY:
          gal_options_check_set(&o[i], &p->naxes[1], GAL_OPTIONS_RANGE_GT_0);
          break;



          /* Operating mode */



        default:
          error(EXIT_FAILURE, 0, "option key %d not recognized in "
                "`ui_read_check_only_options'", o[i].key);
        }

  /* If any of the mandatory options were not given, then print an error
     listing them and abort. */
  if(namell)
    gal_options_mandatory_error(namell, docll);
}





/* Sanity check on options AND arguments. If only option values are to be
   checked, use `ui_read_check_only_options'. */
static void
ui_check_options_and_arguments(struct mkprofparams *p)
{

}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct mkprofparams *p)
{
  struct gal_options_common_params *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero/NULL value for all elements. */
  cp->poptions        = options;
  cp->program_name    = PROGRAM_NAME;
  cp->program_exec    = PROGRAM_EXEC;
  cp->program_bibtex  = PROGRAM_BIBTEX;
  cp->program_authors = PROGRAM_AUTHORS;
  cp->coptions        = gal_commonopts_options;
  cp->numthreads      = num_processors(NPROC_CURRENT);

  /* Read the command-line options and arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");

  /* Read the configuration files. */
  gal_options_read_config_files(cp);

  /* Read the options into the program's structure, and check them and
     their relations prior to printing. */
  ui_read_check_only_options(p);

  /* Print the option values if asked. Note that this needs to be done
     after the sanity check so un-sane values are not printed in the output
     state. */
  gal_options_print_state(cp);

  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments(p);

  /* Free all the allocated spaces in the option structures. */
  gal_options_free(options);
  gal_options_free(gal_commonopts_options);
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
