/*********************************************************************
Table - View and manipulate a FITS table structures.
Table is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <regex.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

#include <gnuastro/fits.h>
#include <gnuastro/table.h>
#include <gnuastro/linkedlist.h>

#include <nproc.h>               /* From Gnulib.                   */
#include <timing.h>              /* Includes time.h and sys/time.h */
#include <checkset.h>
#include <fixedstringmacros.h>

#include "main.h"

#include "ui.h"
#include "args.h"
#include "authors-cite.h"










/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/

/* This function is in charge of reading the option values from the arrays
   of the `argp_option' structure and put them into each program's own data
   structure.

   IMPORTANT: DO NOT USE THE SAME POINTERS. When the types are pointers
   (like stirings or linked lists), do not use the same pointers, allocate
   and keep a new copy. This is because the printing of options will take
   place after this sanity check, so users can be sure that the values that
   are printed and used as configuration files later have no non-sane
   values. Thus the poitners in `argp_option' will be freed after printing,
   but the program's pointers must stay to the end.

   After setting all the values, do other forms of sanity checks that
   involve more than one option.
*/
static void
ui_read_check_only_options(struct tableparams *p)
{
  size_t i;
  struct gal_linkedlist_stll *namell=NULL, *docll=NULL;

  /* Put the program's option values into the structure. */
  for(i=0; !gal_options_is_last(&options[i]); ++i)
    if( options[i].key && options[i].name )
      switch(options[i].key)
        {
        /* Inputs */
        case ARGS_OPTION_COLUMN_KEY:
          gal_linked_list_copy_stll(options[i].value, &p->columns);
          break;


        case ARGS_OPTION_SEARCHIN_KEY:
          if(options[i].value)
            p->searchin=gal_table_string_to_searchin(options[i].value);
          else
            gal_options_add_to_not_given(&options[i], &namell, &docll);
          break;


        case ARGS_OPTION_IGNORECASE_KEY:
          p->ignorecase = *(unsigned char *)options[i].value;
          break;



        /* Output */
        case ARGS_OPTION_TABLETYPE_KEY:

          /* If `tabletype' has a value, put it in. */
          if(options[i].value)
            {
              /* Set the value into the structure. */
              p->tabletype=gal_table_string_to_type(options[i].value);

              /* If the output name was set and is a FITS file, make sure
                 that the type of the table is not a `txt'. */
              if( p->cp.output && gal_fits_name_is_fits(p->cp.output)
                  && ( p->tabletype !=GAL_TABLE_TYPE_AFITS
                       && p->tabletype !=GAL_TABLE_TYPE_BFITS ) )
                error(EXIT_FAILURE, 0, "desired output file `%s' is a FITS "
                      "file, but `tabletype' is not a FITS table type. "
                      "Please set it to `fits-ascii', or `fits-binary'",
                      p->cp.output);
            }

          /* In Table, `tabletype' is only mandatory when the output is a
             FITS file, A text table output only has one possible type. */
          else if( gal_fits_name_is_fits(p->cp.output) )
            gal_options_add_to_not_given(&options[i], &namell, &docll);
          break;



        /* Operating mode */
        case ARGS_OPTION_INFORMATION_KEY:
          if(options[i].value)
            p->information = *(unsigned char *)options[i].value;
          break;


        default:
          error(EXIT_FAILURE, 0, "option key %d not recognized in "
                "`fill_params_from_options'", options[i].key);
        }

  /* If any of the mandatory options were not given, then print an error
     listing them and abort. */
  if(namell)
    gal_options_mandatory_error(namell, docll);
}





static void
ui_check_options_and_arguments(struct tableparams *p)
{
  /* Make sure an input file name was given and if it was a FITS file, that
     a HDU is also given. */
  if(p->up.filename)
    {
      if( gal_fits_name_is_fits(p->up.filename) && p->cp.hdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU specified. When the input is a FITS "
              "file, a HDU must also be specified, you can use the `--hdu' "
              "(`-h') option and give it the HDU number (starting from "
              "zero), extension name, or anything acceptable by CFITSIO");

    }
  else
    error(EXIT_FAILURE, 0, "no input file is specified");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
void
ui_preparations(struct tableparams *p)
{
  char *numstr;
  int tabletype;
  gal_data_t *allcols;
  size_t i, numcols, numrows;

  /* If there were no columns specified, we want the full set of
     columns. */
  if(p->columns==NULL)
    {
      /* Read the table information for the number of columns and rows. */
      allcols=gal_table_info(p->up.filename, p->cp.hdu, &numcols,
                             &numrows, &tabletype);

      /* If there was no actual data in the file, then inform the user */
      if(allcols==NULL)
        error(EXIT_FAILURE, 0, "%s: no usable data rows", p->up.filename);


      /* If the user just wanted information, then print it. */
      if(p->information)
        {
          /* Print the file information. */
          printf("--------\n");
          printf("%s", p->up.filename);
          if(gal_fits_name_is_fits(p->up.filename))
            printf(" (hdu: %s)\n", p->cp.hdu);
          else
            printf("\n");

          /* Print each column's information. */
          gal_table_print_info(allcols, numcols, numrows);
        }


      /* Free the information from all the columns. */
      for(i=0;i<numcols;++i)
        gal_data_free(&allcols[i], 1);
      free(allcols);


      /* If the user just wanted information, then free the allocated
         spaces and exit. Otherwise, add the number of columns to the list
         if the user wanted to print the columns (didn't just want their
         information. */
      if(p->information)
        {
          freeandreport(p);
          exit(EXIT_SUCCESS);
        }
      else
        for(i=1;i<=numcols;++i)
          {
            asprintf(&numstr, "%zu", i);
            gal_linkedlist_add_to_stll(&p->columns, numstr, 0);
          }
    }

  /* Reverse the list of column search criteria that we are looking for
     (since this is a last-in-first-out linked list, the order that
     elements were added to the list is the reverse of the order that they
     will be popped). */
  gal_linkedlist_reverse_stll(&p->columns);
  p->table=gal_table_read(p->up.filename, p->cp.hdu, p->columns,
                          p->searchin, p->ignorecase, p->cp.minmapsize);

  /* If there was no actual data in the file, then inform the user and
     abort. */
  if(p->table==NULL)
    error(EXIT_FAILURE, 0, "%s: no usable data rows (non-commented and "
          "non-blank lines)", p->up.filename);

  /* Now that the data columns are ready, we can free the string linked
     list. */
  gal_linkedlist_free_stll(p->columns, 1);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
setparams(int argc, char *argv[], struct tableparams *p)
{
  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  p->cp.poptions        = options;
  p->cp.program_name    = PROGRAM_NAME;
  p->cp.program_exec    = PROGRAM_EXEC;
  p->cp.program_bibtex  = PROGRAM_BIBTEX;
  p->cp.program_authors = PROGRAM_AUTHORS;
  p->cp.coptions        = gal_commonopts_options;
  p->cp.numthreads      = num_processors(NPROC_CURRENT);

  /* Initialize this utility's pointers to NULL. */
  p->columns=NULL;
  p->up.filename=NULL;

  /* Read the command-line arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");

  /* Read the configuration files. */
  gal_options_read_config_files(&p->cp);

  /* Read the options into the program's structure, and check them and
     their relations prior to printing. */
  ui_read_check_only_options(p);

  /* Print the option values if asked. Note that this needs to be done
     after the sanity check so un-sane values are not printed in the output
     state. */
  gal_options_print_state(&p->cp);

  /* Check that the options and arguments fit well with each other. Note
     that arguments don't go in a configuration file. So this test should
     be done after (possibly) printing the option values. */
  ui_check_options_and_arguments(p);

  /* Read/allocate all the necessary starting arrays. */
  ui_preparations(p);

  /* Free all the allocated spaces in the option structures. */
  gal_options_free(options);
  gal_options_free(gal_commonopts_options);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct tableparams *p)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
  gal_data_free_ll(p->table);
}
