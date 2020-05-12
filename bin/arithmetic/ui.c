/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <string.h>

#include <gnuastro/wcs.h>
#include <gnuastro/list.h>
#include <gnuastro/fits.h>
#include <gnuastro/table.h>
#include <gnuastro/array.h>
#include <gnuastro/threads.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/fixedstringmacros.h>

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
  "(possibly multiple) calls to the '--hdu' option."
  "\n\nCurrently "PROGRAM_NAME" only supports postfix or reverse polish "
  "notation. For example to get the result of '5+6', you should write "
  "'5 6 +', or to get the average of two images, you should write 'a.fits "
  "b.fits + 2 /' (or more simply use the 'average' operator with "
  "'a.fits b.fits average'). Please see the manual for more information. "
  "\n\n"PROGRAM_NAME" recognizes a large collection of standard operators, "
  "including basic arithmetic (e.g., +, -, x, /), mathematical (e.g., abs, "
  "pow, sqrt, log), statistical (minvalue, min, max, average), comparison "
  "(e.g., lt, le, gt), logical (e.g., and, or, not), the full set of bitwise "
  "operators, and numeric type conversion operators to all known types. "
  "Please run the command below for a complete list describing all "
  "operators (press the 'SPACE' keyboard key, or arrow keys, to go down "
  "and 'q' to return to the command-line):\n\n"
  "     $ info gnuastro \"Arithmetic operators\"\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct arithmeticparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;

  /* Set the necessary common parameters structure. */
  cp->program_struct     = p;
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->numthreads         = gal_threads_number();
  cp->coptions           = gal_commonopts_options;

  /* Modify the common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_HDU:
          cp->coptions[i].value=&p->hdus;
          cp->coptions[i].type=GAL_TYPE_STRLL;
          cp->coptions[i].doc="Nth call, used for HDU of Nth input FITS.";
          break;

        case GAL_OPTIONS_KEY_TYPE:
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_IGNORECASE:
        case GAL_OPTIONS_KEY_STDINTIMEOUT:
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;

        case GAL_OPTIONS_KEY_MINMAPSIZE:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          break;
        }

      /* Select by group. */
      switch(cp->coptions[i].group)
        {
        case GAL_OPTIONS_GROUP_TESSELLATION:
          if(cp->coptions[i].key!=GAL_OPTIONS_KEY_INTERPMETRIC)
            cp->coptions[i].flags=OPTION_HIDDEN;
          break;
        }
    }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct arithmeticparams *p = state->input;

  /* Pass 'gal_options_common_params' into the child parser.  */
  state->child_inputs[0] = &p->cp;

  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then 'arg'
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign ('='). For short "
               "options, '=' should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  /* Set the key to this option. */
  switch(key)
    {

    /* Read the non-option tokens (arguments): */
    case ARGP_KEY_ARG:
      gal_list_str_add(&p->tokens, arg, 0);
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
static void
ui_read_check_only_options(struct arithmeticparams *p)
{
  if(p->wcsfile && strcmp(p->wcsfile,"none"))
    {
      if(gal_fits_name_is_fits(p->wcsfile)==0)
        error(EXIT_FAILURE, 0, "%s: file given to '--wcsfile' must be in "
              "FITS format with a recognizable FITS format suffix.",
              p->wcsfile);
      if(p->wcshdu==NULL)
        error(EXIT_FAILURE, 0, "%s: no HDU/extension specified (file given "
              "to '--wcsfile')! Please use '--wcshdu' to specify a "
              "HDU/extension to read from", p->wcsfile);
    }
}





/* Sanity check on options AND arguments. If only option values are to be
   checked, use 'ui_read_check_only_options'. */
static void
ui_check_options_and_arguments(struct arithmeticparams *p)
{
  char *filename;
  int output_checked=0;
  gal_list_str_t *token, *hdu;
  size_t nummultiext=0, numhdus=0;
  struct gal_options_common_params *cp=&p->cp;

  /* First, make sure that any tokens are actually given. */
  if(p->tokens==NULL)
    error(EXIT_FAILURE, 0, "no input tokens. Please specify a filename or "
          "number (as operands) along with operator(s) as input. Please run "
          "any of the following commands for more information.\n\n"
          "    $ astarithmetic --help           # Short info.\n"
          "    $ info astarithmetic             # Full invocation "
          "documentation.\n");

  /* The input tokens are put in a lastin-firstout (simple) linked list, so
     change them to the correct order so the order we pop a token is the
     same order that the user input a value. Note that for the options this
     was done in 'gal_options_read_config_set'. */
  gal_list_str_reverse(&p->tokens);

  /* To allow adding extensions to existing files, let the 'keep' flag be
     the same as the 'dontdelete'. */
  cp->keep=cp->dontdelete;

  /* Set the output file name (if any is needed). Note that since the lists
     are already reversed, the first FITS file encountered, is the first
     FITS file given by the user. Also, note that these file name
     operations are only necessary for the first FITS file in the token
     list. */
  for(token=p->tokens; token!=NULL; token=token->next)
    {
      /* Strings given to the 'tofile' operator are also considered as
         outputs and we should delete them before starting the parse. */
      if( strncmp(OPERATOR_PREFIX_TOFILE, token->v,
                  OPERATOR_PREFIX_LENGTH_TOFILE) )

        {
          /* This token is a file, count how many mult-extension files we
             have and use the first to set the output filename (if it has
             not been set). */
          if( gal_array_name_recognized(token->v) )
            {
              /* Increment the counter for FITS files (if they are
                 input). Recall that the 'tofile' operator can also have
                 '.fits' suffixes (they are the names of the output
                 files). */
              if( gal_array_name_recognized_multiext(token->v)  )
                ++nummultiext;

              /* If the output filename isn't set yet, then set it. */
              if(output_checked==0)
                {
                  if(cp->output)
                    gal_checkset_writable_remove(cp->output, cp->keep,
                                                 cp->dontdelete);
                  else
                    p->cp.output=gal_checkset_automatic_output(cp, token->v,
                                                               "_arith.fits");
                  output_checked=1;
                }
            }

          /* This token is a number. Check if a negative dash was present that
             has been temporarily replaced with 'NEG_DASH_REPLACE' before
             option parsing. */
          else if(token->v[0]==NEG_DASH_REPLACE && isdigit(token->v[1]) )
            token->v[0]='-';
        }

      /* We are on the 'tofile' operator. */
      else
        {
          filename=&token->v[ OPERATOR_PREFIX_LENGTH_TOFILE ];
          gal_checkset_writable_remove(filename, cp->keep, cp->dontdelete);
        }
    }

  /* Count the number of HDU values (if globalhdu isn't given) and check if
     its not less than the number of input FITS images. */
  if(p->globalhdu)
    { if(p->hdus) { gal_list_str_free(p->hdus, 1); p->hdus=NULL; }; }
  else
    {
      for(hdu=p->hdus; hdu!=NULL; hdu=hdu->next) ++numhdus;
      if(numhdus<nummultiext)
        error(EXIT_FAILURE, 0, "not enough HDUs. There are %zu input "
              "files in formats that may contain multiple extensions (for "
              "example FITS or TIFF). Therefore, the '--hdu' ('-h') option "
              "must be called atleaset %zu times (once for each "
              "multi-extension file). If the HDU value is the same for all "
              "the files, you may use '--globalhdu' ('-g') to specify a "
              "single HDU to be used for any number of input files",
              nummultiext, nummultiext);
    }
}





static void
ui_preparations(struct arithmeticparams *p)
{
  size_t ndim, *dsize;

  /* In case a file is specified to read the WCS from (and ignore input
     datasets), read the WCS prior to starting parsing of the arguments. */
  if(p->wcsfile && strcmp(p->wcsfile,"none"))
    {
      /* Read the number of dimensions and the size of each. */
      dsize=gal_fits_img_info_dim(p->wcsfile, p->wcshdu, &ndim);

      /* Read the WCS. */
      p->refdata.wcs=gal_wcs_read(p->wcsfile, p->wcshdu, 0, 0,
                                  &p->refdata.nwcs);
      if(p->refdata.wcs)
        {
          if(!p->cp.quiet)
            printf(" - WCS: %s (hdu %s).\n", p->wcsfile, p->wcshdu);
        }
      else
        fprintf(stderr, "WARNING: %s (hdu %s) didn't contain a "
                "(readable by WCSLIB) WCS.\n", p->wcsfile, p->wcshdu);

      /* Correct the WCS dimensions if necessary. Note that we don't need
         the 'ndim' or 'dsize' any more. */
      ndim=gal_dimension_remove_extra(ndim, dsize, p->refdata.wcs);

      /* Clean up. */
      free(dsize);
    }
}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
ui_read_check_inputs_setup(int argc, char *argv[], struct arithmeticparams *p)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Include the parameters necessary for argp from this program ('args.h')
     and for the common options to all Gnuastro ('commonopts.h'). We want
     to directly put the pointers to the fields in 'p' and 'cp', so we are
     simply including the header here to not have to use long macros in
     those headers which make them hard to read and modify. This also helps
     in having a clean environment: everything in those headers is only
     available within the scope of this function. */
#include <gnuastro-internal/commonopts.h>
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


  /* Print the option values if asked. Note that this needs to be done
     after the sanity check so un-sane values are not printed in the output
     state. */
  gal_options_print_state(cp);


  /* Sanity check only on options. */
  ui_read_check_only_options(p);


  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments(p);


  /* Initial preparations. */
  ui_preparations(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct arithmeticparams *p, struct timeval *t1)
{
  /* Free the simple strings. */
  free(p->cp.output);
  if(p->globalhdu) free(p->globalhdu);

  /* If there are any remaining HDUs in the hdus linked list, then
     free them. */
  if(p->hdus)
    gal_list_str_free(p->hdus, 1);

  /* Report the duration of the job */
  if(!p->cp.quiet)
    gal_timing_report(t1,  PROGRAM_NAME" finished in", 0);
}
