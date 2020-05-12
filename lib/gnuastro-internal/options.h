/*********************************************************************
Function to parse options and configuration file values.

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
#ifndef __GAL_OPTIONS_H__
#define __GAL_OPTIONS_H__

#include <argp.h>

#include <gnuastro/tile.h>
#include <gnuastro/list.h>
#include <gnuastro/fits.h>




/**********************************************************************/
/************                Common options             ***************/
/**********************************************************************/

/* Type for options that don't accept an argument/value. This macro is
   defined to help make the definition and processing of these options
   easier and less buggy. Please use this macro for such options. */
#define GAL_OPTIONS_NO_ARG_TYPE   GAL_TYPE_UINT8

/* When printing the option names, values and comments, we want things to
   be clean and readable (all the comments starting on one line for most,
   ideally all, lines). But in some cases, option values can become too
   long (for example the '--polygon' option in Crop, which takes many
   coordinates). So simply using the maximum option length is going to make
   the whole thing unreadable and we need to have a maximum so this rule
   only applies to them. */
#define GAL_OPTIONS_MAX_VALUE_LEN 10

/* Statically allocated space for printing option values. */
#define GAL_OPTIONS_STATIC_MEM_FOR_VALUES 2000





/* Standard groups of options. From the Argp manual (in GNU C Library): "In
   a long help message, options are sorted alphabetically within each
   group, and the groups presented in the order 0, 1, 2, ..., N, -M, ...,
   -2, -1."

   We want the Operating mode group of options to be the last and the input
   and output groups to the be the first. So first we set the operating
   mdoe group code to '-1', and benefit from the definition of the
   Enumerator type in C, so each field afterwards will be one integer
   larger than the previous. The largest common option group is then used
   to build the codes of program-specific option groups. */
enum options_standard_groups
{
  GAL_OPTIONS_GROUP_OPERATING_MODE = -1,
  GAL_OPTIONS_GROUP_INPUT=1,
  GAL_OPTIONS_GROUP_TESSELLATION,
  GAL_OPTIONS_GROUP_OUTPUT,

  GAL_OPTIONS_GROUP_AFTER_COMMON,
};





/* Key values for the common options, the free alphabetical keys are listed
   below. You can use any of the free letters for an option in a
   program. Note that '-V', which is used by GNU and implemented by Argp,
   is also removed from this list.

   a b c d e f g i j k l m n p r s t u v w x y z
   A B C E G H J L O Q R W X Y
*/
enum options_common_keys
{
  /* With short-option version */
  GAL_OPTIONS_KEY_HDU           = 'h',
  GAL_OPTIONS_KEY_OUTPUT        = 'o',
  GAL_OPTIONS_KEY_TYPE          = 'T',
  GAL_OPTIONS_KEY_DONTDELETE    = 'D',
  GAL_OPTIONS_KEY_KEEPINPUTDIR  = 'K',
  GAL_OPTIONS_KEY_QUIET         = 'q',
  GAL_OPTIONS_KEY_NUMTHREADS    = 'N',
  GAL_OPTIONS_KEY_PRINTPARAMS   = 'P',
  GAL_OPTIONS_KEY_SETDIRCONF    = 'S',
  GAL_OPTIONS_KEY_SETUSRCONF    = 'U',
  GAL_OPTIONS_KEY_IGNORECASE    = 'I',
  GAL_OPTIONS_KEY_TILESIZE      = 'Z',
  GAL_OPTIONS_KEY_NUMCHANNELS   = 'M',
  GAL_OPTIONS_KEY_REMAINDERFRAC = 'F',

  /* Only long option (integers for keywords). */
  GAL_OPTIONS_KEY_STDINTIMEOUT = 500,
  GAL_OPTIONS_KEY_MINMAPSIZE,
  GAL_OPTIONS_KEY_QUIETMMAP,
  GAL_OPTIONS_KEY_LOG,
  GAL_OPTIONS_KEY_CITE,
  GAL_OPTIONS_KEY_CONFIG,
  GAL_OPTIONS_KEY_SEARCHIN,
  GAL_OPTIONS_KEY_LASTCONFIG,
  GAL_OPTIONS_KEY_CHECKCONFIG,
  GAL_OPTIONS_KEY_TABLEFORMAT,
  GAL_OPTIONS_KEY_ONLYVERSION,
  GAL_OPTIONS_KEY_WORKOVERCH,
  GAL_OPTIONS_KEY_CHECKTILES,
  GAL_OPTIONS_KEY_ONEELEMPERTILE,
  GAL_OPTIONS_KEY_INTERPONLYBLANK,
  GAL_OPTIONS_KEY_INTERPMETRIC,
  GAL_OPTIONS_KEY_INTERPNUMNGB,
};





/* Conditions to check */
enum gal_options_range_values
{
  GAL_OPTIONS_RANGE_ANY,

  GAL_OPTIONS_RANGE_GT_0,
  GAL_OPTIONS_RANGE_GE_0,
  GAL_OPTIONS_RANGE_0_OR_1,
  GAL_OPTIONS_RANGE_GE_0_LE_1,
  GAL_OPTIONS_RANGE_GE_0_LT_1,
  GAL_OPTIONS_RANGE_GT_0_LT_1,

  GAL_OPTIONS_RANGE_GT_0_ODD,
  GAL_OPTIONS_RANGE_0_OR_ODD,
};





/* What to do if option isn't given. Note that in each program's 'main.c'
   the main program structure is initialized to zero (or NULL for
   pointers). */
enum gal_options_mandatory_values
{
  GAL_OPTIONS_NOT_MANDATORY,     /* =0 in C standard. */
  GAL_OPTIONS_MANDATORY,
};





/* If the option has already been given a value or not. */
enum gal_options_set_values
{
  GAL_OPTIONS_NOT_SET,           /* =0 in C standard. */
  GAL_OPTIONS_SET,
};





/* The structure keeping all the values of the common options in Gnuastro's
   programs. */
struct gal_options_common_params
{
  /* Tessellation. */
  struct gal_tile_two_layer_params tl; /* Two layer tessellation params.  */
  uint8_t      interponlyblank; /* Only interpolate over blank values.    */
  uint8_t         interpmetric; /* Metric to use for nearest-ngb interp.  */
  size_t          interpnumngb; /* Number of neighbors for interpolation. */

  /* Input. */
  char                    *hdu; /* Image extension.                       */
  uint8_t             searchin; /* Column meta-data to match/search.      */
  uint8_t           ignorecase; /* Ignore case when matching col info.    */
  long            stdintimeout; /* Timeout (micro-seconds) for stdin.     */

  /* Output. */
  char                 *output; /* Directory containg output.             */
  uint8_t                 type; /* Data type of output.                   */
  uint8_t           dontdelete; /* ==1: Don't delete existing file.       */
  uint8_t         keepinputdir; /* Keep input directory for auto output.  */
  uint8_t          tableformat; /* Internal code for output table format. */

  /* Operating modes. */
  uint8_t                quiet; /* Only print errors.                     */
  size_t            numthreads; /* Number of threads to use.              */
  size_t            minmapsize; /* Minimum bytes necessary to use mmap.   */
  uint8_t            quietmmap; /* ==0: print mmap'd file name and size.  */
  uint8_t                  log; /* Make a log file.                       */
  char            *onlyversion; /* Redundant, kept/set for generality.    */

  /* Configuration files. */
  uint8_t          printparams; /* To print the full list of parameters.  */
  uint8_t           setdirconf; /* To write the directory config file.    */
  uint8_t           setusrconf; /* To write teh user config config file.  */
  uint8_t           lastconfig; /* This is the last configuration file.   */
  uint8_t          checkconfig; /* Check config files and values.         */

  /* Output files. */
  gal_fits_list_key_t  *okeys;  /* Configuration as FITS keys in output.  */

  /* For internal (to option processing) purposes. */
  uint8_t                 keep; /* Output file can exist.                 */
  void         *program_struct; /* Host program's main variable struct.   */
  char           *program_name; /* Official name to be used in text.      */
  char           *program_exec; /* Program's executable name.             */
  char         *program_bibtex; /* BibTeX record for this program.        */
  char        *program_authors; /* List of the program authors.           */
  struct argp_option *coptions; /* Common options to all programs.        */
  struct argp_option *poptions; /* Program specific options.              */
  gal_list_i32_t  *mand_common; /* Common mandatory options.  */
  gal_list_str_t  *novalue_doc; /* Mandatory opts, no value   */
  gal_list_str_t *novalue_name; /* Mandatory opts, no value   */
};





/**********************************************************************/
/************              Option utilities             ***************/
/**********************************************************************/
int
gal_options_is_last(struct argp_option *option);

int
gal_options_is_category_title(struct argp_option *option);

void
gal_options_add_to_not_given(struct gal_options_common_params *cp,
                             struct argp_option *option);

void
gal_options_abort_if_mandatory_missing(struct gal_options_common_params *cp);



/**********************************************************************/
/************     Parser functions for common options   ***************/
/**********************************************************************/
void *
gal_options_check_version(struct argp_option *option, char *arg,
                          char *filename, size_t lineno, void *junk);

void *
gal_options_print_citation(struct argp_option *option, char *arg,
                           char *filename, size_t lineno, void *pa);

void *
gal_options_check_config(struct argp_option *option, char *arg,
                         char *filename, size_t lineno, void *junk);

void *
gal_options_read_type(struct argp_option *option, char *arg,
                      char *filename, size_t lineno, void *junk);

void *
gal_options_read_searchin(struct argp_option *option, char *arg,
                          char *filename, size_t lineno, void *junk);

void *
gal_options_read_tableformat(struct argp_option *option, char *arg,
                             char *filename, size_t lineno, void *junk);

void *
gal_options_read_interpmetric(struct argp_option *option, char *arg,
                              char *filename, size_t lineno, void *junk);

gal_data_t *
gal_options_parse_list_of_numbers(char *string, char *filename,
                                  size_t lineno);

gal_data_t *
gal_options_parse_csv_strings_raw(char *string, char *filename,
                                  size_t lineno);

void *
gal_options_parse_csv_strings(struct argp_option *option, char *arg,
                              char *filename, size_t lineno, void *junk);

void *
gal_options_parse_sizes_reverse(struct argp_option *option, char *arg,
                                char *filename, size_t lineno, void *params);

void *
gal_options_parse_csv_float64(struct argp_option *option, char *arg,
                              char *filename, size_t lineno, void *junk);

void *
gal_options_read_sigma_clip(struct argp_option *option, char *arg,
                            char *filename, size_t lineno, void *junk);

void *
gal_options_parse_name_and_strings(struct argp_option *option, char *arg,
                                   char *filename, size_t lineno, void *junk);

void *
gal_options_parse_name_and_float64s(struct argp_option *option, char *arg,
                                    char *filename, size_t lineno, void *junk);

void *
gal_options_parse_colon_sep_csv(struct argp_option *option, char *arg,
                                char *filename, size_t lineno, void *junk);


/**********************************************************************/
/************            Command-line options           ***************/
/**********************************************************************/
error_t
gal_options_set_from_key(int key, char *arg, struct argp_option *options,
                         struct gal_options_common_params *cp);

error_t
gal_options_common_argp_parse(int key, char *arg, struct argp_state *state);

char *
gal_options_stdin_error(long stdintimeout, int precedence, char *name);

gal_list_str_t *
gal_options_check_stdin(char *inputname, long stdintimeout, char *name);




/**********************************************************************/
/************            Configuration files            ***************/
/**********************************************************************/
void *
gal_options_call_parse_config_file(struct argp_option *option, char *arg,
                                   char *filename, size_t lineno, void *cp);

void
gal_options_read_config_set(struct gal_options_common_params *cp);




/**********************************************************************/
/************              Printing/Writing             ***************/
/**********************************************************************/
void
gal_options_print_state(struct gal_options_common_params *cp);

void
gal_options_as_fits_keywords(struct gal_options_common_params *cp);

#endif
