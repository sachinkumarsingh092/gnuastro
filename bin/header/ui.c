/*********************************************************************
Header - View and manipulate a data file header
Header is part of GNU Astronomy Utilities (Gnuastro) package.

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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" print the header "
  "information in any astronomical data file header. It can also "
  "manipulate (add, remove or modify) any of the existing keywords in "
  "a data header. \n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Available letters for short options:

   a b d e f g j k l m n p r s t u v w x y z
   A B C E F G H J L M O Q R T U W X Y Z  */
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_ASIS        = 'a',
  ARGS_OPTION_KEY_DELETE      = 'd',
  ARGS_OPTION_KEY_RENAME      = 'r',
  ARGS_OPTION_KEY_UPDATE      = 'u',
  ARGS_OPTION_KEY_WRITE       = 'w',
  ARGS_OPTION_KEY_COMMENT     = 'c',
  ARGS_OPTION_KEY_HISTORY     = 'h',
  ARGS_OPTION_KEY_DATE        = 't',
  ARGS_OPTION_KEY_QUITONERROR = 'Q',


  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};



















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct headerparams *p,
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
  struct headerparams *p = state->input;

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
ui_read_check_only_options(struct headerparams *p)
{
  /* See if the user just wants to view the header or actually do
     something. */
  if(p->delete || p->updatestr || p->writestr || p->asis || p->comment
     || p->history || p->date || p->rename)
    p->onlyview=0;
  else
    p->onlyview=1;
}





static void
ui_check_options_and_arguments(struct headerparams *p)
{
  /* Make sure an input file name was given and if it was a FITS file, that
     a HDU is also given. */
  if(p->filename)
    {
      if( gal_fits_name_is_fits(p->filename) && p->cp.hdu==NULL )
        error(EXIT_FAILURE, 0, "%s: no HDU specified. A FITS file can "
              "contain, multiple HDUs. You can use the `--hdu' (`-h') "
              "option and give it the HDU number (starting from zero), "
              "extension name, or anything acceptable by CFITSIO",
              p->filename);
    }
  else
    error(EXIT_FAILURE, 0, "no input file is specified");
}




















/**************************************************************/
/*****************       Preparations      ********************/
/**************************************************************/
static void
ui_setup_rename(struct headerparams *p)
{
  char *c;
  struct gal_linkedlist_stll *tmp;

  for(tmp=p->rename; tmp!=NULL; tmp=tmp->next)
    {
      /* `c' is created in case of an error, so the input value can be
         reported. */
      errno=0;
      c=malloc(strlen(tmp->v) + 1);
      if(c==NULL) error(EXIT_FAILURE, errno, "space for c in setuprename");
      strcpy(c, tmp->v);

      /* Tokenize the input. */
      gal_linkedlist_add_to_stll(&p->renamefrom, strtok(tmp->v, ", "), 1);
      gal_linkedlist_add_to_stll(&p->renameto, strtok(NULL, ", "), 1);
      if(p->renamefrom->v==NULL || p->renameto->v==NULL)
        error(EXIT_FAILURE, 0, "`%s' could not be tokenized in order to "
              "complete rename. There should be a space character "
              "or a comma (,) between the two keyword names. If you have "
              "used the space character, be sure to enclose the value to "
              "the `--rename' option in double quotation marks", c);
      free(c);
    }
  /*
  {
    struct gal_linkedlist_stll *tmp2=p->renameto;
    for(tmp=p->renamefrom; tmp!=NULL; tmp=tmp->next)
      {
        printf("%s to %s\n", tmp->v, tmp2->v);
        tmp2=tmp2->next;
      }
  }
  */
}





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
          type=GAL_DATA_TYPE_LONG;
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
              type=GAL_DATA_TYPE_DOUBLE;
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


      gal_fits_add_to_key_ll(output, type, keyname, 0,
                             fvalue, vfree, comment, 0, unit);
      free(original);
    }
}





static void
ui_preparations(struct headerparams *p)
{
  char *ffname;
  int status=0, iomode;

  /* Add hdu to filename: */
  asprintf(&ffname, "%s[%s#]", p->filename, p->cp.hdu);

  /* Open the FITS file: */
  iomode = p->onlyview ? READONLY : READWRITE;
  if( fits_open_file(&p->fptr, ffname, iomode, &status) )
    gal_fits_io_error(status, "reading file");
  free(ffname);

  /* Separate the comma-separated values:  */
  if(p->rename)
    ui_setup_rename(p);
  if(p->updatestr)
    ui_fill_fits_headerll(p->updatestr, &p->update);
  if(p->writestr)
    ui_fill_fits_headerll(p->writestr, &p->write);
}


















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct headerparams *p)
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
ui_free_and_report(struct headerparams *p)
{
  int status=0;

  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);

  /* Close the FITS file: */
  if(fits_close_file(p->fptr, &status))
    gal_fits_io_error(status, NULL);

  if(p->wcs)
    wcsvfree(&p->nwcs, &p->wcs);
}
