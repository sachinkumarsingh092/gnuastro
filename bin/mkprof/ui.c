/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2015, Free Software Foundation, Inc.

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
#include <gnuastro/box.h>
#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/options.h>
#include <gnuastro-internal/checkset.h>
#include <gnuastro-internal/tableintern.h>
#include <gnuastro-internal/fixedstringmacros.h>

#include "main.h"

#include "ui.h"
#include "oneprofile.h"
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
args_doc[] = "[Options] [Catalog]";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will create a FITS image "
  "containing any number of mock astronomical profiles based on an input "
  "catalog. All the profiles will be built from the center outwards. First "
  "by Monte Carlo integration, then using the central pixel position. The "
  "tolerance level specifies when the switch will occur.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




/* Option groups particular to this program. */
enum program_args_groups
{
  ARGS_GROUP_PROFILES = GAL_OPTIONS_GROUP_AFTER_COMMON,
  ARGS_GROUP_CATALOG,
  ARGS_GROUP_WCS,
};




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static int
ui_profile_name_read(char *string, size_t row)
{
  if( !strcmp("sersic", string) )
    return PROFILE_SERSIC;

  else if ( !strcmp("moffat", string) )
    return PROFILE_MOFFAT;

  else if ( !strcmp("gaussian", string) )
    return PROFILE_GAUSSIAN;

  else if ( !strcmp("point", string) )
    return PROFILE_POINT;

  else if ( !strcmp("flat", string) )
    return PROFILE_FLAT;

  else if ( !strcmp("circum", string) )
    return PROFILE_CIRCUMFERENCE;

  else if ( !strcmp(GAL_BLANK_STRING, string) )
    error(EXIT_FAILURE, 0, "atleast one profile function is blank");

  else
    {
      if(row)
        error(EXIT_FAILURE, 0, "`%s' not recognized as a profile function "
              "name in row %zu", string, row);
      else
        error(EXIT_FAILURE, 0, "`%s' not recognized as a profile function "
              "name in values to `--kernel' option", string);
    }

  return PROFILE_INVALID;
}





static char *
ui_profile_name_write(int profile_code)
{
  switch(profile_code)
    {
    case PROFILE_SERSIC:         return "sersic";
    case PROFILE_MOFFAT:         return "moffat";
    case PROFILE_GAUSSIAN:       return "gaussian";
    case PROFILE_POINT:          return "point";
    case PROFILE_FLAT:           return "flat";
    case PROFILE_CIRCUMFERENCE:  return "circum";
    default:
      error(EXIT_FAILURE, 0, "%s: %d not recognized as a profile code",
            __func__, profile_code);
    }

  return NULL;
}






static void
ui_initialize_options(struct mkprofparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;

  /* Set the necessary common parameters structure. */
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->poptions           = program_options;
  cp->numthreads         = gal_threads_number();
  cp->coptions           = gal_commonopts_options;

  /* Default program parameters. */
  p->cp.type=GAL_TYPE_FLOAT32;


  /* Modify the common options for this program. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_HDU:
          cp->coptions[i].doc="Input catalog HDU name or number (if FITS).";
          break;

        case GAL_OPTIONS_KEY_TABLEFORMAT:
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;

        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_MINMAPSIZE:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
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
  struct mkprofparams *p = state->input;

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
      if(p->catname)
        argp_error(state, "only one argument (input catalog) may be given");
      else
        p->catname=arg;
      break;


    /* This is an option, set its value. */
    default:
      return gal_options_set_from_key(key, arg, p->cp.poptions, &p->cp);
    }

  return 0;
}





/* Parse the kernel properties, the format is like this:

     PROFILE_NAME,PARAM_1,PARAM_2,PARAM_3,...,PARAM_N       */
void *
ui_parse_kernel(struct argp_option *option, char *arg,
                char *filename, size_t lineno, void *junk)
{
  long profcode;
  double *darray;
  gal_data_t *kernel;
  size_t i, nc, numneeded=0;
  char *c, *profile, *tailptr;
  char *str, sstr[GAL_OPTIONS_STATIC_MEM_FOR_VALUES];

  /* We want to print the stored values. */
  if(lineno==-1)
    {
      /* Set the value pointer to kernel. */
      kernel=*(gal_data_t **)(option->value);
      darray = kernel->array;

      /* First write the profile function code into the output string. */
      nc=0;
      nc += sprintf(sstr+nc, "%s,", ui_profile_name_write(kernel->status));

      /* Write the values into a string. */
      for(i=0;i<kernel->size;++i)
        {
          if( nc > GAL_OPTIONS_STATIC_MEM_FOR_VALUES-100 )
            error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s so we "
                  "can address the problem. The number of necessary "
                  "characters in the statically allocated string has become "
                  "too close to %d", __func__, PACKAGE_BUGREPORT,
                  GAL_OPTIONS_STATIC_MEM_FOR_VALUES);
          nc += sprintf(sstr+nc, "%g,", darray[i]);
        }
      sstr[nc-1]='\0';

      /* Copy the string into a dynamically allocated space, because it
         will be freed later.*/
      gal_checkset_allocate_copy(sstr, &str);
      return str;
    }
  else
    {
      /* The first part of `arg' (before the first comma) is not
         necessarily a number. So we need to separate the first part from
         the rest.*/
      c=arg;while(*c!='\0' && *c!=',') ++c;
      profile=arg;
      arg = (*c=='\0') ? NULL : c+1;  /* `point' doesn't need any numbers. */
      *c='\0';

      /* Read the parameters. */
      kernel=gal_options_parse_list_of_numbers(arg, filename, lineno);
      *(gal_data_t **)(option->value) = kernel;

      /* Write the profile type code into `kernel->status'. If it starts
         with a digit, then the user might have given the code of the
         profile directly. In that case, parse the number. Otherwise,
         let `ui_profile_name_read' find the value. */
      if( isdigit(*profile) )
        {
          profcode=strtol(profile, &tailptr, 0);
          if(*tailptr!='\0')
            error_at_line(EXIT_FAILURE, 0, filename, lineno, "`%s' "
                          "couldn't be read as a profile code", profile);
          if(profcode<=0 || profcode>=PROFILE_MAXIMUM_CODE)
            error_at_line(EXIT_FAILURE, 0, filename, lineno, "`%s' "
                          "isn't a valid profile code. Please run with "
                          "`--help' and see the acceptable codes in "
                          "explanation of the `--fcol' option", profile);
          kernel->status=profcode;
        }
      else
        kernel->status=ui_profile_name_read(profile, 0);

      /* Make sure the number of parameters conforms with the profile. */
      switch(kernel->status)
        {
        case PROFILE_SERSIC:        numneeded=3;     break;
        case PROFILE_MOFFAT:        numneeded=3;     break;
        case PROFILE_GAUSSIAN:      numneeded=2;     break;
        case PROFILE_POINT:         numneeded=0;     break;
        case PROFILE_FLAT:          numneeded=1;     break;
        case PROFILE_CIRCUMFERENCE: numneeded=1;     break;
        default:
          error_at_line(EXIT_FAILURE, 0, filename, lineno, "%s: a bug! "
                        "Please contact us at %s to correct the issue. "
                        "Profile code %d is not recognized", __func__,
                        PACKAGE_BUGREPORT, kernel->status);
        }

      /* Make sure the number of parameters given are the same number that
         are needed. */
      if( kernel->size != numneeded )
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "as a kernel, a "
                      "`%s' profile needs %zu parameters, but %zu is given",
                      ui_profile_name_write(kernel->status), numneeded,
                      kernel->size);

      /* Our job is done, return NULL. */
      return NULL;
    }
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in `ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct mkprofparams *p)
{
  /* When a no-merged image is to be created, type is necessary. */
  if( p->cp.type==GAL_TYPE_INVALID && p->nomerged==0)
    error(EXIT_FAILURE, 0, "an output type `--type' is necessary when a "
          "merged image is to be built.");

  /* Check if one of the coordinate columns has been given, the other is
     also given. To simplify the job, we use the fact that conditions in C
     return either a 0 (when failed) and 1 (when successful). Note that if
     neighter coordinates are specified there is no problem, the user might
     have input the other coordinate standard. We'll also check for that
     after this.*/
  if(p->kernel==NULL)
    {
      if( ((p->xcol==NULL) + (p->ycol==NULL)) == 1 )
        error(EXIT_FAILURE, 0, "only `%s' has been given, please also "
              "specify a column for the position along the %s axis with "
              "the `%s' option", p->xcol?"xcol":"ycol", p->xcol?"Y":"X",
              p->xcol?"ycol":"xcol");

      if( ((p->racol==NULL) + (p->deccol==NULL)) == 1 )
        error(EXIT_FAILURE, 0, "only `%s' has been given, please also "
              "specify a column for the position along the %s axis with "
              "the `%s' option", p->racol?"racol":"deccol",
              p->racol?"Dec":"RA", p->xcol?"deccol":"racol");
    }
}





/* Sanity check on options AND arguments. If only option values are to be
   checked, use `ui_read_check_only_options'. */
static void
ui_check_options_and_arguments(struct mkprofparams *p)
{
  int d0f1;
  char *tmpname;

  /* If no kernel is given, make sure an input catalog is given, and if it
     is FITS, that the HDU is also provided. When a kernel option, we will
     set a fiducial catalog name called `kernel.txt' to automatic output
     filename generation. */
  if(p->kernel)
    {
      if(p->catname)
        error(EXIT_FAILURE, 0, "`--kernel' cannot be called with an input "
              "catalog (`%s'). The parameters necessary to build a single "
              "kernel output should be given to `--kernel', not in a "
              "catalog", p->catname);
      p->catname="kernel.option";
    }
  else
    {
      if(p->catname)
        {
          if( gal_fits_name_is_fits(p->catname) && p->cp.hdu==NULL)
            error(EXIT_FAILURE, 0, "no `hdu' specified for the input FITS "
                  "table '%s', to ", p->catname);
        }
      else
        error(EXIT_FAILURE, 0, "no input catalog provided. To build "
              "profiles, you need to give a catalog/table containing "
              "the information of the profiles");
    }


  /* If cp->output was not specified on the command line or in any of
     the configuration files, then automatic output should be used, in
     which case, cp->output should be the current directory. */
  if(p->cp.output==NULL)
      gal_checkset_allocate_copy("./", &p->cp.output);


  /* Set the necessary output names. */
  d0f1=gal_checkset_dir_0_file_1(p->cp.output, p->cp.dontdelete);
  if(d0f1)                        /* --output is a file name. */
    {
      p->mergedimgname=p->cp.output;
      p->outdir=gal_checkset_dir_part(p->mergedimgname);
    }
  else                            /* --output is a directory name. */
    {
      gal_checkset_allocate_copy(p->cp.output, &p->outdir);
      gal_checkset_check_dir_write_add_slash(&p->outdir);
      tmpname=gal_checkset_automatic_output(&p->cp, p->catname, ".fits");
      p->mergedimgname=gal_checkset_malloc_cat(p->outdir, tmpname);
      free(tmpname);
    }
  p->basename=gal_checkset_not_dir_part(p->mergedimgname);


  /* If a merged image is requested (or `--kernel' the option is called),
     then delete the final filename if it exists. */
  if(p->nomerged==0 && p->kernel)
    gal_checkset_check_remove_file(p->mergedimgname, p->cp.keep,
                                   p->cp.dontdelete);
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_read_cols(struct mkprofparams *p)
{
  int checkblank;
  size_t counter=0, i;
  char *colname=NULL, **strarr;
  gal_list_str_t *colstrs=NULL;
  gal_data_t *cols, *tmp, *corrtype=NULL;
  char *ax1col=p->racol?p->racol:p->xcol;
  char *ax2col=p->deccol?p->deccol:p->ycol;

  /* Specify the order of columns. */
  gal_list_str_add(&colstrs, ax1col, 0);
  gal_list_str_add(&colstrs, ax2col, 0);
  gal_list_str_add(&colstrs, p->fcol, 0);
  gal_list_str_add(&colstrs, p->rcol, 0);
  gal_list_str_add(&colstrs, p->ncol, 0);
  gal_list_str_add(&colstrs, p->pcol, 0);
  gal_list_str_add(&colstrs, p->qcol, 0);
  gal_list_str_add(&colstrs, p->mcol, 0);
  gal_list_str_add(&colstrs, p->tcol, 0);

  /* Read the desired columns from the file. */
  cols=gal_table_read(p->catname, p->cp.hdu, colstrs, p->cp.searchin,
                      p->cp.ignorecase, p->cp.minmapsize);

  /* Set the number of objects. */
  p->num=cols->size;

  /* For a sanity check, make sure that the total number of columns read is
     the same as those that were wanted (it might be more). */
  while(cols!=NULL)
    {
      /* Pop out the top column. */
      tmp=gal_list_data_pop(&cols);

      /* By default check if the column has blank values, but it can be
         turned off for some columns. */
      checkblank=1;

      /* Note that the input was a linked list, so the output order is the
         inverse of the input order. For the position, we will store the
         values into the `x' and `y' arrays even if they are RA/Dec. */
      switch(++counter)
        {
        case 9:
          colname="first axis position";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT64);
          p->x=corrtype->array;
          break;

        case 8:
          colname="second axis position";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT64);
          p->y=corrtype->array;
          break;

        case 7:
          if(tmp->type==GAL_TYPE_STRING)
            {
              p->f=gal_data_malloc_array(GAL_TYPE_INT32, p->num,
                                         __func__, "p->f");
              strarr=tmp->array;
              for(i=0;i<p->num;++i)
                p->f[i]=ui_profile_name_read(strarr[i], i+1);
              gal_data_free(tmp);
              corrtype=NULL;
            }
          else
            {
              /* Read the user's profile codes. */
              colname="profile function code (`fcol')";
              corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_INT32);
              p->f=corrtype->array;

              /* Check if they are in the correct range. */
              for(i=0;i<p->num;++i)
                if(p->f[i]<=PROFILE_INVALID || p->f[i]>=PROFILE_MAXIMUM_CODE)
                  error(EXIT_FAILURE, 0, "%s: table row %zu, the function "
                        "code is %d. It should be >%d and <%d. Please run "
                        "again with `--help' and check the acceptable "
                        "codes.\n\nAlternatively, you can use alphabetic "
                        "strings to specify the profile functions, see the "
                        "explanations under `fcol' from the command "
                        "below (press the `SPACE' key to go down, and the "
                        "`q' to return back to the command-line):\n\n"
                        "    $ info %s\n", p->catname, i+1, p->f[i],
                        PROFILE_INVALID, PROFILE_MAXIMUM_CODE, PROGRAM_EXEC);
            }
          break;

        case 6:
          colname="radius (`rcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->r=corrtype->array;
          break;

        case 5:
          colname="index (`ncol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->n=corrtype->array;
          break;

        case 4:
          colname="position angle (`pcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->p=corrtype->array;
          break;

        case 3:
          colname="axis ratio (`qcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->q=corrtype->array;
          break;

        case 2:
          colname="magnitude (`mcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->m=corrtype->array;
          checkblank=0;          /* Magnitude can be NaN: to mask regions. */
          break;

        case 1:
          colname="truncation (`tcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->t=corrtype->array;
          break;

        /* If the index isn't recognized, then it is larger, showing that
           there was more than one match for the given criteria */
        default:
          gal_tableintern_error_col_selection(p->catname, p->cp.hdu, "too "
                                              "many columns were selected "
                                              "by the given values to the "
                                              "options ending in `col'.");
        }

      /* Sanity check and clean up.  Note that it might happen that the
         input structure is already freed. In that case, `corrtype' will be
         NULL. */
      if(corrtype)
        {
          /* Make sure there are no blank values in this column. */
          if( checkblank && gal_blank_present(corrtype, 1) )
            error(EXIT_FAILURE, 0, "%s column has blank values. "
                  "Input columns cannot contain blank values", colname);

          /* Free the unnecessary sturcture information. The correct-type
             (`corrtype') data structure's array is necessary for later
             steps, so its pointer has been copied in the main program's
             structure. Hence, we should set the structure's pointer to
             NULL so the important data isn't freed.*/
          corrtype->array=NULL;
          gal_data_free(corrtype);
        }
    }
}





/* It is possible to define the internal catalog through a catalog or the
   `--kernel' option. This function will do the job. */
static void
ui_prepare_columns(struct mkprofparams *p)
{
  float r, n, t;
  double *karr;

  /* If the kernel option was called, then we need to build a series of
     single element columns to create an internal catalog. */
  if(p->kernel)
    {
      /* Number of profiles to be built. */
      p->num=1;

      /* Allocate the necessary columns. */
      p->x=gal_data_malloc_array(GAL_TYPE_FLOAT64, 1, __func__, "p->x");
      p->y=gal_data_malloc_array(GAL_TYPE_FLOAT64, 1, __func__, "p->y");
      p->f=gal_data_malloc_array(GAL_TYPE_UINT8,   1, __func__, "p->f");
      p->r=gal_data_malloc_array(GAL_TYPE_FLOAT32, 1, __func__, "p->r");
      p->n=gal_data_malloc_array(GAL_TYPE_FLOAT32, 1, __func__, "p->n");
      p->p=gal_data_malloc_array(GAL_TYPE_FLOAT32, 1, __func__, "p->p");
      p->q=gal_data_malloc_array(GAL_TYPE_FLOAT32, 1, __func__, "p->q");
      p->m=gal_data_malloc_array(GAL_TYPE_FLOAT32, 1, __func__, "p->m");
      p->t=gal_data_malloc_array(GAL_TYPE_FLOAT32, 1, __func__, "p->t");

      /* Set the values that need special consideration. */
      if(p->kernel->size)
        {
          karr=p->kernel->array;
          r = karr[0];
          n = p->kernel->size==2 ? 0.0f : karr[1];
          t = p->kernel->size==1 ? 1.0f : karr[ p->kernel->size - 1 ];
        }
      else r=n=t=0.0f;

      /* Fill the allocated spaces. */
      p->x[0] = 0.0f;
      p->y[0] = 0.0f;
      p->f[0] = p->kernel->status;
      p->r[0] = r;
      p->n[0] = n;
      p->p[0] = 0.0f;
      p->q[0] = 1.0f;
      p->m[0] = 0.0f;
      p->t[0] = t;
    }
  else
    ui_read_cols(p);
}





static void
ui_prepare_wcs(struct mkprofparams *p)
{
  int status;
  struct wcsprm *wcs;
  long os=p->oversample;

  /* If the WCS structure is already set, then return. */
  if(p->out->wcs) return;

  /* Allocate the memory necessary for the wcsprm structure. */
  errno=0;
  wcs=p->out->wcs=malloc(sizeof *wcs);
  if(wcs==NULL)
    error(EXIT_FAILURE, errno, "%zu for wcs in preparewcs", sizeof *wcs);

  /* Initialize the structure (allocate all its internal arrays). */
  wcs->flag=-1;
  if( (status=wcsini(1, 2, wcs)) )
    error(EXIT_FAILURE, 0, "wcsini error %d: %s",
          status, wcs_errmsg[status]);

  /* Correct the CRPIX values based on oversampling and shifting. */
  p->crpix[0] = p->crpix[0]*os + p->shift[0] - os/2;
  p->crpix[1] = p->crpix[1]*os + p->shift[1] - os/2;

  /* Fill in all the important WCS structure parameters. */
  wcs->altlin   = 0x1;
  wcs->equinox  = 2000.0f;
  wcs->crpix[0] = p->crpix[0];
  wcs->crpix[1] = p->crpix[1];
  wcs->crval[0] = p->crval[0];
  wcs->crval[1] = p->crval[1];
  wcs->pc[0]    = -1.0f;
  wcs->pc[3]    = 1.0f;
  wcs->pc[1]    = wcs->pc[2]=0.0f;
  wcs->cdelt[0] = wcs->cdelt[1]=p->resolution/3600;
  strcpy(wcs->cunit[0], "deg");
  strcpy(wcs->cunit[1], "deg");
  strcpy(wcs->ctype[0], "RA---TAN");
  strcpy(wcs->ctype[1], "DEC--TAN");

  /* Set up the wcs structure with the constants defined above. */
  status=wcsset(wcs);
  if(status)
    error(EXIT_FAILURE, 0, "wcsset error %d: %s", status,
          wcs_errmsg[status]);
}





static void
ui_prepare_canvas(struct mkprofparams *p)
{
  int status=0;
  float *f, *ff;
  double truncr;
  long width[2]={1,1};
  size_t i, ndim, dsize[2];

  /* If a background image is specified, then use that as the output
     image to build the profiles over. */
  if(p->backname)
    {
      /* Make sure the kernel option is not called. */
      if(p->kernel)
        error(EXIT_FAILURE, 0, "the `--kernel' and `--background' options "
              "cannot be called together");

      /* Small sanity check. */
      if(p->backhdu==NULL)
        error(EXIT_FAILURE, 0, "no hdu specified for the background image "
              "%s. Please run again `--backhdu' option", p->backname);

      /* Read in the background image and its coordinates, note that when
         no merged image is desired, we just need the WCS information of
         the background image. */
      if(p->nomerged)
        p->out=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 0, dsize, NULL,
                              1, p->cp.minmapsize, NULL, NULL, NULL);
      else
        {
          /* Read the image. */
          p->out=gal_fits_img_read_to_type(p->backname, p->backhdu,
                                           GAL_TYPE_FLOAT32,
                                           p->cp.minmapsize);
          p->naxes[0]=p->out->dsize[1];
          p->naxes[1]=p->out->dsize[0];

          /* Set all pixels to zero if the user wanted a clear canvas. */
          if(p->clearcanvas)
            {ff=(f=p->out->array)+p->out->size; do *f++=0.0f; while(f<ff);}
        }

      /* When a background image is specified, oversample must be 1 and
         there is no shifts. */
      p->oversample=1;
      p->shift[0]=p->shift[1]=0;

      /* Read the WCS structure of the background image. */
      p->out->wcs=gal_wcs_read(p->backname, p->backhdu, 0, 0, &p->out->nwcs);

    }
  else
    {
      /* If any of xshift or yshift is non-zero, the other should be too!
         Note that conditional operators return 1 if true and 0 if false,
         so if one is non-zero while the other is zero, then sum will be
         1. Otherwise the sum will either be 0 or 2.*/
      switch ( (p->shift[0]!=0) + (p->shift[1]!=0) )
        {
        case 0:
          /* `prepforconv' is only valid when xshift and yshift are both
             zero. Also, a PSF profile should exist in the image. */
          if(p->prepforconv)
            {
              /* Check if there is at least one Moffat or Gaussian profile. */
              for(i=0;i<p->num;++i)
                if( oneprofile_ispsf(p->f[i]) )
                  {
                    /* Calculate the size of the box holding the PSF. Note:

                       - For the Moffat and Gaussian profiles, the radius
                       columns is actually the FWHM which is actually the
                       diameter, not radius. So we have to divide it by
                       half.

                       - encloseellipse outputs the total width, we only want
                       half of it for the shift. */
                    truncr = p->tunitinp ? p->t[i] : p->t[i] * p->r[i]/2;
                    gal_box_ellipse_in_box(truncr, p->q[i]*truncr,
                                           p->p[i]*DEGREESTORADIANS, width);
                    p->shift[0]  = (width[0]/2)*p->oversample;
                    p->shift[1]  = (width[1]/2)*p->oversample;
                  }
            }
          break;

        case 1:
          error(EXIT_FAILURE, 0, "at least one of `--xshift` (`-X`) or "
                "`--yshift` (`-Y`) are zero, they must either both be zero "
                "or both given a positive value");
          break;

        case 2:
          p->shift[0] *= p->oversample;
          p->shift[1] *= p->oversample;
          break;

        default:
          error(EXIT_FAILURE, 0, "a bug in ui_prepare_canvas! In checks "
                "for shifts. Please contact us at %s so we can fix it",
                PACKAGE_BUGREPORT);
        }

      /* Prepare the sizes of the final merged image (if it is to be
         made). */
      if(p->nomerged)
        ndim=0;
      else
        {
          /* Sanity check */
          if(p->naxes[0]==0 || p->naxes[1]==0)
            error(EXIT_FAILURE, 0, "No final merged image size is specified, "
                  "please use the `--naxis1', and `--naxis2' options to "
                  "specify the respective sizes");

          /* Set the final merged image size. Note that even if we don't
             want a merged image, we still need its WCS structure.*/
          p->naxes[0] = (p->naxes[0] * p->oversample) + (2 * p->shift[0]);
          p->naxes[1] = (p->naxes[1] * p->oversample) + (2 * p->shift[1]);
          dsize[0]    = p->naxes[1];
          dsize[1]    = p->naxes[0];
          ndim        = 2;
        }

      /* Make the output structure. */
      p->out=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, ndim, dsize, NULL, 1,
                            p->cp.minmapsize, NULL, NULL, NULL);
    }


  /* Make the WCS structure of the output data structure if it has not
     been set yet. */
  ui_prepare_wcs(p);


  /* Set the name, comments and units of the output structure/file. Note
     that when no merged image is to be created, the array in `p->out' will
     be NULL.*/
  if(p->out->array)
    {
      if(p->out->name) free(p->out->name);
      gal_checkset_allocate_copy("Mock profiles", &p->out->name);
      if(p->out->unit==NULL)
        gal_checkset_allocate_copy("Brightness", &p->out->unit);
    }



  /* When individual mode is requested, write the WCS structure to a header
     string to speed up the process: if we don't do it here, this process
     will be necessary on every individual profile's output. So it is much
     more efficient done once here. */
  if(p->individual && p->out->wcs)
    {
      status=wcshdo(WCSHDO_safe, p->out->wcs, &p->wcsnkeyrec, &p->wcsheader);
      if(status)
        error(EXIT_FAILURE, 0, "wcshdo error %d: %s", status,
              wcs_errmsg[status]);
    }
}





static void
ui_finalize_coordinates(struct mkprofparams *p)
{
  size_t i;
  double *x=NULL, *y=NULL;

  /* When the user specified RA and Dec columns, the respective values
     where stored in the `p->x' and `p->y' arrays. So before proceeding, we
     need to change them into actual image coordinates. */
  if(p->racol)
    {
      /* Note that we read the RA and Dec columns into the `p->x' and `p->y'
         arrays temporarily before. Here, we will convert them, free the old
         ones and replace them with the proper X and Y values. */
      gal_wcs_world_to_img(p->out->wcs, p->x, p->y, &x, &y, p->num);

      /* If any conversions created a WCSLIB error, both the outputs will be
         set to NaN. */
      for(i=0;i<p->num;++i)
        if( isnan(x[i]) )
          error(EXIT_FAILURE, 0, "catalog row %zu: WCSLIB could not convert "
                "(%f, %f) coordinates into image coordinates", i, p->x[i],
                p->y[i]);

      /* Free the RA and Dec arrays and put in the new image values. */
      free(p->x);
      free(p->y);
      p->x=x;
      p->y=y;
    }


  /* Correct the WCS scale. Note that when the WCS is read from a
     background image, oversample is set to 1. This is done here because
     the conversion of WCS to pixel coordinates needs to be done with the
     non-over-sampled image. */
  p->out->wcs->cdelt[0] /= p->oversample;
  p->out->wcs->cdelt[1] /= p->oversample;



  /* For a sanity check:
  printf("\nui_finalize_coordinates sanity check:\n");
  for(i=0;i<p->num;++i)
    printf("%f, %f\n", p->x[i], p->y[i]);
  */
}





/* Add all the columns of the log file. Just note that since this is a
   linked list, we have to add them in the opposite order. */
static void
ui_make_log(struct mkprofparams *p)
{
  char *name, *comment;

  /* Return if no long file is to be created. */
  if(p->cp.log==0) return;

  /* Individual created. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_UINT8, 1, &p->num, NULL,
                          1, p->cp.minmapsize, "INDIV_CREATED", "bool",
                          "If an individual image was made (1) or not (0).");

  /* Fraction of monte-carlo. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_FLOAT32, 1, &p->num, NULL,
                          1, p->cp.minmapsize, "FRAC_MONTECARLO", "frac",
                          "Fraction of brightness in Monte-carlo integrated "
                          "pixels.");

  /* Number of monte-carlo. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_UINT64, 1, &p->num, NULL,
                          1, p->cp.minmapsize, "NUM_MONTECARLO", "count",
                          "Number of Monte Carlo integrated pixels.");

  /* Magnitude of profile overlap. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_FLOAT32, 1, &p->num, NULL,
                          1, p->cp.minmapsize, "MAG_OVERLAP", "mag",
                          "Magnitude of profile's overlap with merged "
                          "image.");

  /* Row number in input catalog. */
  name=gal_fits_name_save_as_string(p->catname, p->cp.hdu);
  asprintf(&comment, "Row number of profile in %s.", name);
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_UINT64, 1, &p->num, NULL,
                          1, p->cp.minmapsize, "INPUT_ROW_NO", "count",
                          comment);
  free(comment);
  free(name);
}





static void
ui_preparations(struct mkprofparams *p)
{
  /* If the kernel option was given, some parameters need to be
     over-written: */
  if(p->kernel)
    {
      p->nomerged=1;
      p->psfinimg=0;
      p->individual=1;
    }

  /* Read in all the columns. */
  ui_prepare_columns(p);

  /* Prepare the output canvas. */
  ui_prepare_canvas(p);

  /* Read the (possible) RA/Dec inputs into X and Y for the builder.*/
  ui_finalize_coordinates(p);

  /* Allocate the random number generator: */
  gsl_rng_env_setup();
  p->rng=gsl_rng_alloc(gsl_rng_default);

  /* Make the log linked list. */
  ui_make_log(p);
}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
static void
ui_print_intro(struct mkprofparams *p)
{
  char *jobname;

  if(p->cp.quiet) return;

  printf(PROGRAM_NAME" started on %s", ctime(&p->rawtime));

  if(p->kernel)
    asprintf(&jobname, "Building one %s kernel",
             ui_profile_name_write(p->kernel->status));
  else
    asprintf(&jobname, "%zu profile%sread from %s", p->num,
             p->num>1?"s ":" ", p->catname);
  gal_timing_report(NULL, jobname, 1);
  free(jobname);

  if(p->backname)
    {
      if(p->nomerged)
        asprintf(&jobname, "WCS information read from %s", p->backname);
      else
        asprintf(&jobname, "%s is read and will be used as canvas",
                 p->backname);
      gal_timing_report(NULL, jobname, 1);
      free(jobname);
    }

  asprintf(&jobname, "Random number generator (RNG) type: %s",
           gsl_rng_name(p->rng));
  gal_timing_report(NULL, jobname, 1);
  free(jobname);
  if(p->envseed)
    {
      asprintf(&jobname, "RNG seed for all profiles: %lu",
               gsl_rng_default_seed);
      gal_timing_report(NULL, jobname, 1);
      free(jobname);
    }

  if(p->kernel==NULL)
    {
      asprintf(&jobname, "Using %zu threads.", p->cp.numthreads);
      gal_timing_report(NULL, jobname, 1);
      free(jobname);
    }
}





void
ui_read_check_inputs_setup(int argc, char *argv[], struct mkprofparams *p)
{
  struct gal_options_common_params *cp=&p->cp;


  /* Include the parameters necessary for argp from this program (`args.h')
     and for the common options to all Gnuastro (`commonopts.h'). We want
     to directly put the pointers to the fields in `p' and `cp', so we are
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


  /* Read the configuration files. */
  gal_options_read_config_set(&p->cp);


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

  /* Print introductory information. */
  ui_print_intro(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct mkprofparams *p, struct timeval *t1)
{
  /* Free all the allocated arrays. */
  free(p->cat);
  free(p->cp.hdu);
  free(p->outdir);
  free(p->basename);

  /* p->cp.output might be equal to p->mergedimgname. In this case, if
     we simply free them after each other, there will be a double free
     error. So after freeing output, we set it to NULL since
     free(NULL) is ok.*/
  if(p->cp.output==p->mergedimgname)
    free(p->cp.output);
  else
    {
      free(p->cp.output);
      free(p->mergedimgname);
    }

  /* Free the WCS headers string that was defined for individual mode. */
  if(p->individual)
    free(p->wcsheader);

  /* Free the random number generator: */
  gsl_rng_free(p->rng);

  /* Free the log file information. */
  if(p->cp.log)
    gal_list_data_free(p->log);

  /* Report the duration of the job */
  if(!p->cp.quiet)
    gal_timing_report(t1,  PROGRAM_NAME" finished in", 0);
}
