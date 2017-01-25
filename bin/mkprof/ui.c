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

#include "main.h"

#include "ui.h"
#include "args.h"









/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read the options into the main program structure. When an option wasn't
   given, it will not be given a value here, so it will have the
   initialized value of 0 (or NULL for pointers) after this function. If
   the value of `0' is meaningful in the context of the option, then it
   must be given the blank value for the type of the variable IN THIS
   FUNCTION.

   When the option is necessary for the program to run (independent of any
   arguments) but it wasn't given, call the `gal_options_add_to_not_given'
   function. The role of the final `gal_options_abort_if_mandatory_missing'
   function at the end of this function is to print the full list of
   mandatory options that were gathered by that function and abort with one
   error message to help the user. */
static void
ui_read_options(struct mkprofparams *p)
{
  size_t i;

  /* Put the program's option values into the structure. */
  for(i=0; !gal_options_is_last(&options[i]); ++i)
    if( options[i].key && options[i].name )
      switch(options[i].key)
        {

          /* Input */
        case ARGS_OPTION_BACKHDU_KEY:
          gal_checkset_allocate_copy(options[i].value, &p->backhdu);
          break;


          /* Output */

          /* naxis1 and naxis2 are only mandatory when a background is not
             given and `individual' is not called.  */
        case ARGS_OPTION_NAXIS1_KEY:
          printf("\n... just before checking `naxis1' ...\n");
          exit(0);
          break;

        case ARGS_OPTION_NAXIS2_KEY:

          break;



          /* Operating mode */



        default:
          error(EXIT_FAILURE, 0, "option key %d not recognized in "
                "`ui_read_check_only_options'", options[i].key);
        }


  /* If any of the mandatory options were not given, then print an error
     listing them and abort. */
  gal_options_abort_if_mandatory_missing(&p->cp);
}





/* Read and check ONLY the options. When arguments are involved, do the
   check in `ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct mkprofparams *p)
{

  /* Read all the options from the `argp_option' array into the main
     program structure to facilitate checks and running the program. */
  ui_read_options(p);

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
ui_read_check_inputs_setup(int argc, char *argv[], struct mkprofparams *p)
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


  /* Read the number of threads available to the user, this should be done
     before reading command-line and configuration file options, since they
     can change it.  */
  gal_options_initialize_numthreads(cp);


  /* If there are mandatory common options, add them to this list. The
     mandatory options for this program will be checked in
     `ui_read_check_only_options' and finally they will all be reported
     together if any of them are not given a value. */
  gal_linkedlist_add_to_ill(&cp->mand_common, GAL_OPTIONS_MINMAPSIZE_KEY);


  /* Read the command-line options and arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");


  /* Read the configuration files. */
  gal_options_read_config_set_common(cp);


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
  gal_linkedlist_free_ill(cp->mand_common);
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
