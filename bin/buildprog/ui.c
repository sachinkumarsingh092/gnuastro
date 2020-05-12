/*********************************************************************
BuildProgram: Compile and run programs using Gnuastro's library
BuildProgram is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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

#include <gnuastro/fits.h>

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
/* Definition parameters for the Argp: */
const char *
argp_program_version = PROGRAM_STRING "\n"
                       GAL_STRINGS_COPYRIGHT
                       "\n\nWritten/developed by "PROGRAM_AUTHORS;

const char *
argp_program_bug_address = PACKAGE_BUGREPORT;

static char
args_doc[] = "C-source [ARGUMENTS TO RUN]";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will compile and run a "
  "C program, while automatically linking with libraries that Gnuastro "
  "depends on. Hence you do not have to worry about explicitly linking "
  "with CFITSIO for example if you want to work on a FITS file, or with "
  "GSL if you want to use GNU Scientific Library's functions. The standard "
  "compiler options of '-I', '-L', and '-l' are also available for further "
  "customization of the build.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct buildprogparams *p,
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
  cp->coptions           = gal_commonopts_options;


  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_LOG:
        case GAL_OPTIONS_KEY_HDU:
        case GAL_OPTIONS_KEY_TYPE:
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_NUMTHREADS:
        case GAL_OPTIONS_KEY_TABLEFORMAT:
        case GAL_OPTIONS_KEY_STDINTIMEOUT:
          cp->coptions[i].flags=OPTION_HIDDEN;
          cp->coptions[i].mandatory=GAL_OPTIONS_NOT_MANDATORY;
          break;

        /* '--ignorecase's default short format is 'I', but here we want to
           follow the compiler format, hence we need 'I' for
           'include'. Therefore, here, we'll change the key for 'include'
           to some large number just to avoid confusion.*/
        case GAL_OPTIONS_KEY_IGNORECASE:
          cp->coptions[i].key=20000;
          cp->coptions[i].flags=OPTION_HIDDEN;
          cp->coptions[i].mandatory=GAL_OPTIONS_NOT_MANDATORY;
        }



      /* Select by group. */
      switch(cp->coptions[i].group)
        {
        case GAL_OPTIONS_GROUP_TESSELLATION:
          cp->coptions[i].doc=NULL; /* Necessary to remove title. */
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;
        }
    }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct buildprogparams *p = state->input;

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
      gal_list_str_add(&p->sourceargs, arg, 0);
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
   check in 'ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct buildprogparams *p)
{
  size_t len;

  /* If an '.la' file is given, make sure it has the correct suffix. */
  if(p->la)
    {
      len=strlen(p->la);
      if(len>=4)
        if(strcmp(&p->la[len-3], ".la"))
          error(EXIT_FAILURE, 0, "'%s' is not a Libtool control file name "
                "(with a '.la' suffix). The file name given to the '--la' "
                "('-a') option must be a Libtool control file", p->la);
    }
}





/* Check the options and arguments. */
static void
ui_check_options_and_arguments(struct buildprogparams *p)
{
  if(p->sourceargs==NULL)
    error(EXIT_FAILURE, 0, "no input (C source file) given");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/

void
ui_preparations(struct buildprogparams *p)
{
  /* Reverse the sourceargs list (note that the options are reversed in
     options.c). */
  gal_list_str_reverse(&p->sourceargs);

  /* Set the final output name. 'EXEEXT' comes from the configuration
     script (given by BuildProgram's 'Makefile.am'). */
  if(p->cp.output==NULL)
    p->cp.output=gal_checkset_automatic_output(&p->cp, p->sourceargs->v,
                                               EXEEXT);

  /* Set the C compiler. Later we can add a check to make sure that 'cc' is
     actually in the PATH. */
  if(p->cc==NULL)
    {                                        /* No C compiler chosen. */
      if(p->noenv==0)
        {
          p->cc=getenv("CC");                /* First check for 'CC'. */
          if(p->cc==NULL)
            p->cc=getenv("GCC");             /* Then check for 'GCC'. */
        }
      if(p->cc==NULL) p->cc="gcc";           /* Default: 'gcc'.       */
    }
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct buildprogparams *p)
{
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


  /* Read the command-line options and arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");


  /* Read the configuration files and set the common values. */
  gal_options_read_config_set(&p->cp);


  /* Read the options into the program's structure, and check them and
     their relations prior to printing. */
  ui_read_check_only_options(p);


  /* Print the option values if asked. Note that this needs to be done
     after the option checks so un-sane values are not printed in the
     output state. */
  gal_options_print_state(&p->cp);


  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments(p);


  /* Read/allocate all the necessary starting arrays. */
  ui_preparations(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct buildprogparams *p)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
  gal_list_str_free(p->include,    1);
  gal_list_str_free(p->linkdir,    1);
  gal_list_str_free(p->linklib,    1);
  gal_list_str_free(p->sourceargs, 0);
}
