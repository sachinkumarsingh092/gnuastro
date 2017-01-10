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


extern struct argp_option gal_commonopts_options[];


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

  /* Internal (before control goes back to the program). */
  int           cite;     /* BibTeX citation for program and abort.      */
  int    printparams;     /* Print all option values and abort.          */
  int     setdirconf;     /* Set default values for this dir and abort.  */
  int     setusrconf;     /* Set default values for this user and abort. */
  int    onlydirconf;     /* Only read current directory config file.    */
  int    onlyversion;     /* Only run program with this version.         */
};





/**********************************************************************/
/************         Main user functions/macros        ***************/
/**********************************************************************/

int
gal_options_is_last(struct argp_option *option);

int
gal_options_is_category_title(struct argp_option *option);





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
gal_options_parse_configs(char *progname, struct argp_option *progopts,
                          struct argp_option *commopts);


#endif
