/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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
#include <gnuastro/linkedlist.h>

#include <options.h>
#include <checkset.h>
#include <fixedstringmacros.h>

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
args_doc[] = "ASTRdata";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" view and manipulate (add, "
  "delete, or modify) FITS extensions and header keywords. \n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Option groups particular to this program. */
enum program_args_groups
{
  ARGS_GROUP_EXTENSION = GAL_OPTIONS_GROUP_AFTER_COMMON,
  ARGS_GROUP_KEYWORD,
};

















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct fitsparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;

  /* Set the necessary common parameters structure. */
  cp->keep               = 1;
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->coptions           = gal_commonopts_options;

  /* For clarity and non-zero initializations. */
  p->mode                = FITS_MODE_INVALID;

  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_IGNORECASE:
        case GAL_OPTIONS_KEY_TYPE:
        case GAL_OPTIONS_KEY_TABLEFORMAT:
        case GAL_OPTIONS_KEY_DONTDELETE:
        case GAL_OPTIONS_KEY_LOG:
        case GAL_OPTIONS_KEY_MINMAPSIZE:
        case GAL_OPTIONS_KEY_NUMTHREADS:
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;

        case GAL_OPTIONS_KEY_OUTPUT:
          cp->coptions[i].doc="Output file name (only for writing HDUs).";
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
  struct fitsparams *p = state->input;

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
      /* Only FITS files are acceptable. */
      if( gal_fits_name_is_fits(arg) )
        {
          if(p->filename)
            argp_error(state, "only one input file should be given");
          else
            p->filename=arg;
        }
      else
        argp_error(state, "%s is not a recognized FITS file", arg);
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
ui_read_check_only_options(struct fitsparams *p)
{
  /* If any of the keyword manipulation options are requested, then set the
     mode flag to keyword-mode. */
  if( p->date || p->comment || p->history || p->asis || p->delete
      || p->rename || p->rename || p->update || p->write || p->printallkeys )
    {
      /* Set the mode. */
      p->mode=FITS_MODE_KEY;

      /* Check if a HDU is given. */
      if(p->cp.hdu==NULL)
        error(EXIT_FAILURE, 0, "a HDU (extension) is necessary for keywrod "
              "related options but none was defined. Please use the "
              "`--hdu' (or `-h') option to select one");
    }

  /* Same for the extension-related options */
  if( p->remove || p->copy || p->cut )
    {
      /* A small sanity check. */
      if(p->mode!=FITS_MODE_INVALID)
        error(EXIT_FAILURE, 0, "extension and keyword manipulation options "
              "cannot be called together");

      /* Set the mode and turn on the `needshdu' flag. */
      p->mode=FITS_MODE_HDU;

      /* Make sure the output name is set. */
      if(p->cp.output)
        gal_checkset_check_remove_file(p->filename, 1, p->cp.dontdelete);
      else
        p->cp.output=gal_checkset_automatic_output(&p->cp, p->filename,
                                                   "_ext.fits");
    }

  /* If no options are given, go into HDU mode, which will print the HDU
     information when nothing is asked. */
  if(p->mode==FITS_MODE_INVALID)
    p->mode=FITS_MODE_HDU;
}





static void
ui_check_options_and_arguments(struct fitsparams *p)
{
  /* Make sure an input file name was given and if it was a FITS file, that
     a HDU is also given. */
  if(p->filename==NULL)
    error(EXIT_FAILURE, 0, "no input file is specified");
}




















/**************************************************************/
/*****************       Preparations      ********************/
/**************************************************************/
/* The `--update' and `--write' options take multiple values for each
   keyword, so here, we tokenize them and put them into a `gal_fits_key_ll'
   list. */
static void
ui_fill_fits_headerll(struct gal_linkedlist_stll *input,
                          struct gal_fits_key_ll **output)
{
  long l, *lp;
  void *fvalue;
  double d, *dp;
  int i=0, type, vfree;
  char *c, *cf, *start, *tailptr;
  struct gal_linkedlist_stll *tmp;
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

      /* Find the of the value: */
      errno=0;
      l=strtol(value, &tailptr, 10);
      if(*tailptr=='\0' && errno==0)
        {
          vfree=1;
          type=GAL_DATA_TYPE_INT64;
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
              type=GAL_DATA_TYPE_FLOAT64;
              errno=0;
              fvalue=dp=malloc(sizeof *dp);
              if(dp==NULL)
                error(EXIT_FAILURE, errno, "%zu bytes for double",
                      sizeof *dp);
              *dp=d;
            }
          else
            { fvalue=value; type=GAL_DATA_TYPE_STRING; vfree=0; }
        }


      gal_fits_key_add_to_ll(output, type, keyname, 0,
                             fvalue, vfree, comment, 0, unit);
      free(original);
    }
}





static void
ui_preparations(struct fitsparams *p)
{
  /* Fill in the key linked lists. We want to do this here so if there is
     any error in parsing the user's input, the error is reported before
     any change is made in the input file. */
  if(p->write)  ui_fill_fits_headerll(p->write, &p->write_keys);
  if(p->update) ui_fill_fits_headerll(p->update, &p->update_keys);
}


















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct fitsparams *p)
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
ui_free_and_report(struct fitsparams *p)
{
  /* Free the allocated arrays: */
  free(p->cp.output);
}
