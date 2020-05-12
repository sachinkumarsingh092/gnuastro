/*********************************************************************
Fits - View and manipulate FITS extensions and/or headers.
Fits is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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
args_doc[] = "ASTRdata";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" allows you to view and "
  "manipulate (add, delete, or modify) FITS extensions (or HDUs) and FITS "
  "header keywords within one extension.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















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
  cp->program_struct     = p;
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
        case GAL_OPTIONS_KEY_NUMTHREADS:
        case GAL_OPTIONS_KEY_STDINTIMEOUT:
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
static void
ui_check_copykeys(struct fitsparams *p)
{
  long read;
  char *tailptr;
  /* size_t group=0; */
  char forl='f', *pt=p->copykeys;

  /* For copykeys, an output filename is mandatory. */
  if(p->cp.output==NULL || p->outhdu==NULL)
    error(EXIT_FAILURE, 0, "an output FITS extension (in an existing "
          "FITS file, specified with the '--output' and '--outhdu') are "
          "mandatory for running '--copykeys'");

  /* Initialize the values. */
  p->copykeysrange[0]=p->copykeysrange[1]=GAL_BLANK_LONG;

  /* Parse the string: 'forl': "first-or-last". */
  while(*pt!='\0')
    {
      switch(*pt)
        {
        case ':':
          forl='l';
          ++pt;
          break;
        case '.':
          error(EXIT_FAILURE, 0, "the numbers in the argument to "
                "'--section' ('-s') have to be integers. You input "
                "includes a float number: %s", p->copykeys);
          break;
        case ' ': case '\t':
          ++pt;
          break;

        /* Numerical characters signify the start of a number, so we don't
           need to increment the pointer and can just break out. */
        case '0': case '1': case '2': case '3': case '4': case '5':
        case '6': case '7': case '8': case '9': case '-':
          break;

        /* Identifier for next group of ranges. However, For the time
           being, we just support one group. So we are commenting the break
           here for it to follow onto default.
        case ',':
          ++group;
          forl='f';
          ++pt;
          break;
          */
        default:
          error(EXIT_FAILURE, 0, "value to '--copykeys' must only contain "
                "integer numbers and these special characters between them: "
                "':' when necessary. But it is '%s' (the first "
                "non-acceptable character is '%c').\n", p->copykeys, *pt);
          break;
        }

      /* Read the number: */
      read=strtol(pt, &tailptr, 0);

      /* Check the progress.
      printf("\n\n------\n%c: %ld (%s)\n", *pt, read, tailptr);
      */

      /* Make sure if a number was read at all? */
      if(tailptr==pt) continue;   /* No number was read!                 */

      /* Put it in the correct place. */
      p->copykeysrange[ forl=='f' ? 0 : 1 ]=read;
      pt=tailptr;
    }

  /* Basic sanity checks. */
  if( p->copykeysrange[1]==GAL_BLANK_LONG )
    error(EXIT_FAILURE, 0, "no ending keyword number given to '--copykeys'. "
          "If you want to copy all the keywords after a certain one "
          "(without worrying about how many there are), you can use '-1'.\n\n"
          "For example if you want to copy all the keywords after the 20th, "
          "you can use '--copykeys=20,-1'. Generally, you can use negative "
          "numbers for the last keyword number to count from the end.");
  if( p->copykeysrange[0]<=0 )
    error(EXIT_FAILURE, 0, "the first number given to '--copykeys' must be "
          "positive");
  if( p->copykeysrange[1]>=0 && p->copykeysrange[0]>=p->copykeysrange[1] )
    error(EXIT_FAILURE, 0, "the first number (%ld) given to '--copykeys' "
          "must be smaller than the second (%ld)", p->copykeysrange[0],
          p->copykeysrange[1]);

  /* For a check:
  printf("copykeys: %ld, %ld\n", p->copykeysrange[0], p->copykeysrange[1]);
  exit(0);
  */
}





/* Read and check ONLY the options. When arguments are involved, do the
   check in 'ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct fitsparams *p)
{
  uint8_t stdoutcheck;

  /* If any of the keyword manipulation options are requested, then set the
     mode flag to keyword-mode. */
  if( p->date || p->comment || p->history || p->asis || p->delete
      || p->rename || p->update || p->write || p->verify || p->printallkeys
      || p->copykeys || p->datetosec )
    {
      /* Check if a HDU is given. */
      if(p->cp.hdu==NULL)
        error(EXIT_FAILURE, 0, "a HDU (extension) is necessary for keyword "
              "related options but none was defined. Please use the "
              "'--hdu' (or '-h') option to select one");

      /* If Copy keys has been given, read it and make sure its setup. */
      if(p->copykeys)
        ui_check_copykeys(p);

      /* Currently 'datetosec' must be called alone. */
      if( p->datetosec
          && (p->date || p->comment || p->history || p->asis || p->delete
              || p->rename || p->update || p->write || p->verify
              || p->printallkeys || p->copykeys) )
        error(EXIT_FAILURE, 0, "'--datetosec' cannot currently be called "
              "with any other option");

      /* Set the operating mode. */
      p->mode=FITS_MODE_KEY;
    }

  /* Same for the extension-related options */
  if( p->remove || p->copy || p->cut || p->numhdus || p->datasum )
    {
      /* A small sanity check. */
      if(p->mode!=FITS_MODE_INVALID)
        error(EXIT_FAILURE, 0, "extension and keyword manipulation options "
              "cannot be called together");

      /* Some HDU options cannot be called with other options. */
      stdoutcheck = 0;
      stdoutcheck = p->numhdus + p->datasum;

      /* Make sure the output name is set. */
      if(stdoutcheck)
        {
          /* Make sure the other HDU-related options aren't called. */
          if(p->remove || p->copy || p->cut)
            error(EXIT_FAILURE, 0, "'--numhdus' or '--datasum' options "
                  "must be called alone");

          /* Make sure the HDU is given for the datasum option. */
          if( p->datasum && p->cp.hdu==NULL )
            error(EXIT_FAILURE, 0, "a HDU (extension) is necessary for the "
                  " '--checksum' or '--datasum' options. Please use the "
                  "'--hdu' (or '-h') option to select one");
        }
      else
        {
          if(p->cp.output)
            gal_checkset_writable_remove(p->cp.output, 1, p->cp.dontdelete);
          else
            p->cp.output=gal_checkset_automatic_output(&p->cp, p->filename,
                                                       "_ext.fits");
        }

      /* Set the operating mode. */
      p->mode=FITS_MODE_HDU;
    }

  /* If no options are given, go into HDU mode, which will print the HDU
     information when nothing is asked. */
  if(p->mode==FITS_MODE_INVALID)
    {
      if(p->hdu_in_commandline)
        {
          p->printallkeys=1;
          p->mode = FITS_MODE_KEY;
        }
      else
        p->mode = FITS_MODE_HDU;
    }
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
/* The '--update' and '--write' options take multiple values for each
   keyword, so here, we tokenize them and put them into a
   'gal_fits_list_key_t' list. */
static void
ui_fill_fits_headerll(gal_list_str_t *input, gal_fits_list_key_t **output,
                      char *option_name)
{
  long l, *lp;
  void *fvalue;
  double d, *dp;
  gal_list_str_t *tmp;
  char *c, *cf, *start, *tailptr;
  int i=0, type, vfree, needsvalue;
  char *original, *keyname, *value, *comment, *unit;

  for(tmp=input; tmp!=NULL; tmp=tmp->next)
    {
      i=0;
      tailptr=NULL;

      /* 'c' is created in case of an error, so the input value can be
         reported. */
      errno=0;
      original=malloc(strlen(tmp->v)+1);
      if(original==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'original'",
              __func__, strlen(tmp->v)+1);
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

      /* See if this is an option that needs a value or not.*/
      needsvalue=1;
      if(keyname)
        {
          if( !strcasecmp(keyname,"checksum")
              || !strcasecmp(keyname,"datasum") )
            needsvalue=0;
        }

      /* Make sure the keyname and value (if necessary) is given. */
      if( keyname==NULL || (needsvalue && value==NULL) )
        error(EXIT_FAILURE, 0, "'--%s' option string (%s) can't be parsed. "
              "The general expected format is (a comment string and unit "
              "are optional):\n\n"
              "    --%s=KEYWORD,value,\"a comment string\",unit\n\n"
              "Any space characters around the the comma (,) characters "
              "will be seen as part of the respective token.\n\n"
              "Note that there are some exceptions (where no value is need)"
              "please see the manual for more ('$ info astfits')",
              option_name, original, option_name);
      /*
      printf("\n\n-%s-\n-%s-\n-%s-\n-%s-\n", keyname, value, comment, unit);
      */


      /* Find the type of the value: */
      if(value)
        {
          errno=0;
          l=strtol(value, &tailptr, 10);
          if(*tailptr=='\0' && errno==0)
            {
              vfree=1;
              type=GAL_TYPE_INT64;
              errno=0;
              fvalue=lp=malloc(sizeof *lp);
              if(lp==NULL)
                error(EXIT_FAILURE, errno, "%s: %zu bytes for 'lp'",
                      __func__, sizeof *lp);
              *lp=l;
            }
          else
            {
              errno=0;
              d=strtod(value, &tailptr);
              if(*tailptr=='\0' && errno==0)
                {
                  vfree=1;
                  type=GAL_TYPE_FLOAT64;
                  errno=0;
                  fvalue=dp=malloc(sizeof *dp);
                  if(dp==NULL)
                    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes "
                          "for 'dp'", __func__, sizeof *dp);
                  *dp=d;
                }
              else
                { fvalue=value; type=GAL_TYPE_STRING; vfree=0; }
            }
        }
      else
        {
          fvalue=NULL; type=GAL_TYPE_UINT8; vfree=0;
        }


      /* Add it to the list of keywords. */
      gal_fits_key_list_add(output, type, keyname, 0, fvalue, vfree,
                            comment, 0, unit);
      free(original);
    }

  /* Reverse the list */
  gal_fits_key_list_reverse(output);
}





static void
ui_preparations(struct fitsparams *p)
{
  /* Fill in the key linked lists. We want to do this here so if there is
     any error in parsing the user's input, the error is reported before
     any change is made in the input file. */
  if(p->write)  ui_fill_fits_headerll(p->write, &p->write_keys, "write");
  if(p->update) ui_fill_fits_headerll(p->update, &p->update_keys, "update");
}


















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct fitsparams *p)
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


  /* Read the command-line options and arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");


  /* Check if the HDU is specified on the command-line. If so, then later,
     if no operation is requested, we will print the header of the given
     HDU.*/
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    if(cp->coptions[i].key==GAL_OPTIONS_KEY_HDU
       && cp->coptions[i].set)
      p->hdu_in_commandline=1;


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
