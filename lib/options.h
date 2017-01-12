/*********************************************************************
Function to parse options and configuration file values.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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
#ifndef __GAL_OPTIONS_H__
#define __GAL_OPTIONS_H__

#include <argp.h>

#include <gnuastro/data.h>






/**********************************************************************/
/************                Common options             ***************/
/**********************************************************************/

/* Type for options that don't accept an argument/value. This macro is
   defined to help make the definition and processing of these options
   easier and less buggy. Please use this macro for such options. */
#define GAL_OPTIONS_NO_ARG_TYPE GAL_DATA_TYPE_UCHAR


/* External global variables that must be specified by the program using
   this header. */
extern char program_name[];     /* Defined in program's `main.h' */
extern char program_exec[];     /* Defined in program's `main.h' */
extern char program_bibtex[];   /* Defined in program's `cite.h' */
extern struct argp_option gal_commonopts_options[];


/* Key values for the common options, the free alphabetical keys are listed
   below. You can use any of the free letters for an option in a
   program. Note that `-V', which is used by GNU and implemented by Argp,
   is also removed.

   a b c d e f g i j k l m n p r s t u v w x y z
   A B C E F G H I J L M O Q R T U W X Y Z
*/
enum options_option_keys
{
  /* With short-option version */
  GAL_OPTIONS_HDU_KEY          = 'h',
  GAL_OPTIONS_OUTPUT_KEY       = 'o',
  GAL_OPTIONS_DONTDELETE_KEY   = 'D',
  GAL_OPTIONS_KEEPINPUTDIR_KEY = 'K',
  GAL_OPTIONS_QUIET_KEY        = 'q',
  GAL_OPTIONS_NUMTHREADS_KEY   = 'N',
  GAL_OPTIONS_PRINTPARAMS_KEY  = 'P',
  GAL_OPTIONS_SETDIRCONF_KEY   = 'S',
  GAL_OPTIONS_SETUSRCONF_KEY   = 'U',

  /* Only long option (integers for keywords). */
  GAL_OPTIONS_MINMAPSIZE_KEY   = 500,
  GAL_OPTIONS_LOG_KEY,
  GAL_OPTIONS_CITE_KEY,
  GAL_OPTIONS_CONFIG_KEY,
  GAL_OPTIONS_LASTCONFIG_KEY,
  GAL_OPTIONS_ONLYVERSION_KEY,
};


/* The structure keeping all the values of the common options in Gnuastro's
   programs. */
struct gal_options_common_params
{
  /* Input/Output: */
  char          *hdu;     /* Image extension.                            */
  char       *output;     /* Directory containg output.                  */
  int     dontdelete;     /* ==1: Don't delete existing file.            */
  int   keepinputdir;     /* Keep input directory for automatic output.  */

  /* Operating modes: */
  int          quiet;     /* ==1: don't print anything but errors.       */
  size_t  numthreads;     /* Number of threads to use.                   */
  size_t  minmapsize;     /* The minimum bytes necessary to use mmap.    */
  int            log;     /* Make a log file.                            */
};





/**********************************************************************/
/************         Main user functions/macros        ***************/
/**********************************************************************/

int
gal_options_is_last(struct argp_option *option);

int
gal_options_is_category_title(struct argp_option *option);

void
gal_options_free(struct argp_option *options);




/**********************************************************************/
/************            Command-line options           ***************/
/**********************************************************************/
error_t
gal_options_set_from_key(int key, char *arg, struct argp_option *options);

error_t
gal_options_common_argp_parse(int key, char *arg, struct argp_state *state);





/**********************************************************************/
/************            Configuration files            ***************/
/**********************************************************************/
void
gal_options_config_files(struct argp_option *poptions,
                         struct gal_options_common_params *cp);

void
gal_options_print_state(struct argp_option *poptions);


#endif
