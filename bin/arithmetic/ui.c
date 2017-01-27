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
args_doc[] = "ASTRdata or number [ASTRdata] OPERATOR ...";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will do arithmetic "
  "operations on one or multiple images and numbers. Simply put, the name "
  "of the image along with the arithmetic operators and possible numbers "
  "are given as arguments. The extensions of each input are specified with "
  "(possibly multiple) calls to the `--hdu' option."
  "\n\nCurrently "PROGRAM_NAME" only supports postfix or reverse polish "
  "notation. For example to get the result of `5+6', you should write "
  "`5 6 +', or to get the average of two images, you should write `a.fits "
  "b.fits + 2 /' (or more simply use the `average' operator with "
  "`a.fits b.fits average'). Please see the manual for more information. "
  "\n\nThe operators/functions recognized by "PROGRAM_NAME" are: +, -, *, "
  "/, abs, pow, sqrt, log, log10, minvalue, maxvalue, min, max, average, "
  "median, lt, le, gt, ge, eq, ne, and, or, not, isblank, and the full set "
  "of bitwise operators. Please run `info gnuastro \"Arithmetic "
  "operators\"' for detailed information on each operator. Note that "
  "multiplication should be quoted (like \"*\", or '*') to avoid shell "
  "expansion.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   a b c d e f g i j k l m n p r s t u v w x y z
   A B C E F G H I J L M O Q R T U W X Y Z                */
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_HDU        = 'h',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct imgarithparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  struct gal_options_common_params *cp=&p->cp;

  /* Set the necessary common parameters structure. */
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->coptions           = gal_commonopts_options;
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct imgarithparams *p = state->input;

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
      gal_linkedlist_add_to_stll(&p->tokens, arg, 0);
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
ui_read_check_only_options(struct imgarithparams *p)
{

}





/* Sanity check on options AND arguments. If only option values are to be
   checked, use `ui_read_check_only_options'. */
static void
ui_check_options_and_arguments(struct imgarithparams *p)
{
  int output_checked=0;
  size_t numfits=0, numhdus=0;
  struct gal_linkedlist_stll *token, *hdu;

  /* The input tokens are put in a lastin-firstout (simple) linked list, so
     change them to the correct order so the order we pop a token is the
     same order that the user input a value. Note that for the options this
     was done in `gal_options_read_config_set'. */
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


  /* Include the parameters necessary for argp from this program (`args.h')
     and for the common options to all Gnuastro (`commonopts.h'). We want
     to directly put the pointers to the fields in `p' and `cp', so we are
     simply including the header here to not have to use long macros in
     those headers which make them hard to read. */
#define NOT_COMMON_HDU_PARSER 1
#include <commonopts.h>
#include "args.h"


  /* Initialize the options and necessary information.  */
  ui_initialize_options(p, program_options, gal_commonopts_options);


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
  gal_options_read_config_set(cp);


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
