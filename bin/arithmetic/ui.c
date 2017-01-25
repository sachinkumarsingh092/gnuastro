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
ui_read_options(struct imgarithparams *p)
{
  size_t i;

  /* Put the program's option values into the structure. */
  for(i=0; !gal_options_is_last(&options[i]); ++i)
    if( options[i].key && options[i].name )
      switch(options[i].key)
        {
        /* Inputs */
        case ARGS_OPTION_HDU_KEY:
          gal_linked_list_copy_stll(options[i].value, &p->hdus);
          break;



        /* Output */



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
ui_read_check_only_options(struct imgarithparams *p)
{

  /* Read all the options from the `argp_option' array into the main
     program structure to facilitate checks and running the program. */
  ui_read_options(p);

}





/* Sanity check on options AND arguments. If only option values are to be
   checked, use `ui_read_check_only_options'. */
static void
ui_check_options_and_arguments(struct imgarithparams *p)
{
  int output_checked=0;
  size_t numfits=0, numhdus=0;
  struct gal_linkedlist_stll *token, *hdu;

  /* The inputs are put in a lastin-firstout (simple) linked list, so
     change them to the correct order so the order we pop a token is the
     same order that the user input a value. */
  gal_linkedlist_reverse_stll(&p->hdus);
  gal_linkedlist_reverse_stll(&p->tokens);

  /* Set the output file name (if any is needed). Note that since the
     lists are already reversed, the first FITS file encountered, is
     the first FITS file given by teh user. Also, notet that these
     file name operations are only necessary for the first FITS file
     in the token list. */
  for(token=p->tokens; token!=NULL; token=token->next)
    {
      /* This token is a FITS file, count them and use it to set the output
         filename if it has not been set. */
      if(gal_fits_name_is_fits(token->v))
        {
          /* Increment the counter for FITS files. */
          ++numfits;

          /* If the output filename isn't set yet, then set it. */
          if(output_checked==0)
            {
              if(p->cp.output)
                gal_checkset_check_remove_file(p->cp.output,
                                               p->cp.dontdelete);
              else
                p->cp.output=gal_checkset_automatic_output(&p->cp, token->v,
                                                           "_arith.fits");
              output_checked=1;
            }
        }

      /* This token is a number. Check if a negative dash was present that
         has been temporarily replaced with `NEG_DASH_REPLACE' before
         option parsing. */
      else if(token->v[0]==NEG_DASH_REPLACE && isdigit(token->v[1]) )
        token->v[0]='-';
    }

  /* Count the number of HDU values and check if its not less than the
     number of input FITS images. */
  for(hdu=p->hdus; hdu!=NULL; hdu=hdu->next) ++numhdus;
  if(numhdus<numfits)
    error(EXIT_FAILURE, 0, "not enough HDUs. There are %zu input FITS "
          "files, but only %zu HDUs. You can use the `--hdu' (`-h') option "
          "to specify the number or name of a HDU for each FITS file",
          numfits, numhdus);
}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
ui_read_check_inputs_setup(int argc, char *argv[], struct imgarithparams *p)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Set the non-zero initial values, the structure was initialized to
     have a zero/NULL value for all elements. */
  cp->poptions        = options;
  cp->program_name    = PROGRAM_NAME;
  cp->program_exec    = PROGRAM_EXEC;
  cp->program_bibtex  = PROGRAM_BIBTEX;
  cp->program_authors = PROGRAM_AUTHORS;
  cp->coptions        = gal_commonopts_options;


  /* The dash of a negative number will cause problems with the option
     readin. To work properly we will go over all the options/arguments and
     if any one starts with a dash and is followed by a number, then the
     dash is replaced by NEG_DASH_REPLACE. */
  for(i=0;i<argc;++i)
    if(argv[i][0]=='-' && isdigit(argv[i][1]))
      argv[i][0]=NEG_DASH_REPLACE;


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
freeandreport(struct imgarithparams *p, struct timeval *t1)
{
  free(p->cp.output);

  /* If there are any remaining HDUs in the hdus linked list, then
     free them. */
  if(p->hdus)
    gal_linkedlist_free_stll(p->hdus, 1);

  /* Report the duration of the job */
  if(!p->cp.quiet)
    gal_timing_report(t1,  PROGRAM_NAME" finished in", 0);
}
