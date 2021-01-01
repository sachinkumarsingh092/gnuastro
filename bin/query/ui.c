/*********************************************************************
Query - Retreive data from a remote data server.
Query is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2020-2021, Free Software Foundation, Inc.

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
#include <gnuastro/pointer.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/fixedstringmacros.h>

#include "main.h"

#include "ui.h"
#include "query.h"
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
args_doc[] = "ASTRdata";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" is just a place holder "
  "used as a minimal set of files and functions necessary for a program in "
  "Gnuastro. It can be used for learning or as a template to build new "
  "programs.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct queryparams *p,
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
        case GAL_OPTIONS_KEY_TYPE:
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_QUIETMMAP:
        case GAL_OPTIONS_KEY_IGNORECASE:
        case GAL_OPTIONS_KEY_NUMTHREADS:
        case GAL_OPTIONS_KEY_MINMAPSIZE:
        case GAL_OPTIONS_KEY_STDINTIMEOUT:
        case GAL_OPTIONS_KEY_KEEPINPUTDIR:
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;
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
  struct queryparams *p = state->input;

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
      argp_error(state, "no input arguments are needed");
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}





char *
ui_strlist_to_str(gal_list_str_t *input)
{
  char *out=NULL;
  gal_list_str_t *node;
  size_t n, nn, nnodes=0, alllen=0;

  /* First calculate the full length of all nodes. */
  for(node=input; node!=NULL; node=node->next)
    {
      /* We'll add two extra for each. One for the ',' that must come in
         between it and the next one. One just for a buffer, incase we
         haven't accounted for something. */
      alllen += strlen(node->v) + 2;
      ++nnodes;
    }

  /* Allocate the output string. */
  out=gal_pointer_allocate(GAL_TYPE_STRING, alllen, 1, "out", __func__);

  /* Write all the strings into the allocated space. */
  n=nn=0;
  for(node=input; node!=NULL; node=node->next)
    {
      if(nn++==nnodes-1)
        sprintf(out+n, "%s", node->v);
      else
        n += sprintf(out+n, "%s,", node->v);
    }

  /* Return the merged string. */
  return out;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in 'ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct queryparams *p)
{
  size_t i;
  gal_data_t *tmp;

  /* See if database has been specified. */
  if(p->databasestr==NULL)
    error(EXIT_FAILURE, 0, "no input dataset.\n\n"
          "Please use the '--database' ('-d') option to specify your "
          "desired database, see manual ('info gnuastro astquery' "
          "command) for the current databases");

  /* Convert the given string into a code. */
  if( !strcmp(p->databasestr, "gaia") ) p->database=QUERY_DATABASE_GAIA;
  else
    error(EXIT_FAILURE, 0, "'%s' is not a recognized database.\n\n"
          "For the full list of recognized databases, please see the "
          "documentation (with the command 'info astquery')", p->databasestr);


  /* Make sure that '--query' and '--center' are not called together. */
  if(p->query && (p->center || p->overlapwith) )
    error(EXIT_FAILURE, 0, "the '--query' option cannot be called together "
          "together with '--center' or '--overlapwith'");

  /* Overlapwith cannot be called with the manual query. */
  if( p->overlapwith && (p->center || p->width || p->radius) )
    error(EXIT_FAILURE, 0, "the '--overlapwith' option cannot be called "
          "with the manual region specifiers ('--center', '--width' or "
          "'--radius')");

  /* The radius and width cannot be called together. */
  if(p->radius && p->width)
    error(EXIT_FAILURE, 0, "the '--radius' and '--width' options cannot be "
          "called together");

  /* If radius is given, it should be one value and positive. */
  if(p->radius)
    {
      if(p->radius->size>1)
        error(EXIT_FAILURE, 0, "only one value can be given to '--radius' "
              "('-r') option");

      if( ((double *)(p->radius->array))[0]<0 )
        error(EXIT_FAILURE, 0, "the '--radius' option value cannot be negative");
    }

  /* If magnitude is given, it should only be two values. */
  i=0;
  if(p->range)
    for(tmp=p->range; tmp!=NULL; tmp=tmp->next)
      {
        ++i;
        if(tmp->size!=2)
          error(EXIT_FAILURE, 0, "two values (separated by ',' or ':') "
                "should be given to '--range'. But %zu values were given "
                "to the %zu%s call of this option (recall that the first "
                "value should be the column name in the given dataset)",
                tmp->size, i,
                i==1 ? "st" : i==2 ? "nd" : i==3 ? "rd" : "th");
      }

  /* Sanity checks on  width (if we are in the center-mode). */
  if(p->width && p->center)
    {
      /* Width should have the same number of elements as the center
         coordinates */
      if( p->width->size > 1 && p->width->size != p->center->size )
        error(EXIT_FAILURE, 0, "'--width' should either have a single "
              "value (used for all dimensions), or one value for each "
              "dimension. However, you have provided %zu coordinate "
              "values, and %zu width values", p->center->size,
              p->width->size);

      /* All values must be positive. */
      for(i=0;i<p->width->size;++i)
        if( ((double *)(p->width->array))[i]<0 )
          error(EXIT_FAILURE, 0, "the '--width' option value(s) cannot "
                "be negative");
    }

  /* Make sure that the output name is in a writable location and that it
     doesn't exist. If it exists, and the user hasn't called
     '--dontdelete', then delete the existing file. */
  gal_checkset_writable_remove(p->cp.output, p->cp.keep,
                               p->cp.dontdelete);

  /* Set the name for the downloaded and processed files. These are due to
     an internal low-level processing that will be done on the raw
     downloaded file. */
  p->processedname=gal_checkset_make_unique_suffix(p->cp.output
                                                   ? p->cp.output
                                                   : "query.fits",
                                                   ".fits");
  p->downloadname=gal_checkset_make_unique_suffix(p->processedname, NULL);
}





static void
ui_check_options_and_arguments(struct queryparams *p)
{

}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_preparations(struct queryparams *p)
{

}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
ui_read_check_inputs_setup(int argc, char *argv[], struct queryparams *p)
{
  struct gal_options_common_params *cp=&p->cp;


  /* Just to avoid warning on no minmapsize. It is irrelevant here. */
  p->cp.minmapsize=-1;


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


  /* Prepare all the options as FITS keywords to write in output later. */
  gal_options_as_fits_keywords(&p->cp);


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
ui_free_report(struct queryparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
}
