/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <gnuastro/array.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>
#include <gnuastro/pointer.h>

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




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static uint8_t
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

  else if ( !strcmp("distance", string) )
    return PROFILE_DISTANCE;

  else if ( !strcmp(GAL_BLANK_STRING, string) )
    error(EXIT_FAILURE, 0, "atleast one profile function is blank");

  else
    {
      if(row)
        error(EXIT_FAILURE, 0, "'%s' not recognized as a profile function "
              "name in row %zu", string, row);
      else
        error(EXIT_FAILURE, 0, "'%s' not recognized as a profile function "
              "name in values to '--kernel' option", string);
    }

  return PROFILE_INVALID;
}





char *
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
    case PROFILE_DISTANCE:       return "distance";
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
  cp->program_struct     = p;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->poptions           = program_options;
  cp->numthreads         = gal_threads_number();
  cp->coptions           = gal_commonopts_options;

  /* Default program parameters. */
  p->zeropoint           = NAN;
  p->cp.type             = GAL_TYPE_FLOAT32;

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
  size_t i, nc, need=0;
  char *c, *dstr, *profile, *tailptr;
  char *str, sstr[GAL_OPTIONS_STATIC_MEM_FOR_VALUES];

  /* We want to print the stored values. */
  if(lineno==-1)
    {
      /* Set the value pointer to kernel. */
      kernel=*(gal_data_t **)(option->value);
      darray = kernel->array;

      /* First write the profile function code into the output string. */
      nc=0;
      profile=ui_profile_name_write(kernel->status);
      switch(kernel->flag)
        {
        case 2: nc += sprintf(sstr+nc, "%s,",    profile); break;
        case 3: nc += sprintf(sstr+nc, "%s-3d,", profile); break;
        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                "fix the problem. %u is not a recognized kernel "
                "dimensionality", __func__, PACKAGE_BUGREPORT, kernel->flag);
        }

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
      /* The first part of 'arg' (before the first comma) is not
         necessarily a number. So we need to separate the first part from
         the rest.*/
      c=arg;while(*c!='\0' && *c!=',') ++c;
      profile=arg;
      arg = (*c=='\0') ? NULL : c+1;  /* the 'point' profile doesn't need */
      *c='\0';                        /* any numbers.                     */


      /* Read the parameters. */
      kernel=gal_options_parse_list_of_numbers(arg, filename, lineno);
      *(gal_data_t **)(option->value) = kernel;


      /* All parameters must be positive. */
      darray=kernel->array;
      for(i=0;i<kernel->size;++i)
        if(darray[i]<=0)
          error(EXIT_FAILURE, 0, "value number %zu (%g) in the given list "
                "of kernel parameters ('%s') is not acceptable. All "
                "parameters to the '--kernel' option must be non-zero and "
                "positive", i+1, darray[i], arg);


      /* See if a 2D kernel is requested or a 3D kernel and keep the value
         in 'kernel->flag'. If no dimensionality is defined, then by
         default, we'll assume it is 2D.*/
      c=profile;while(*c!='\0' && *c!='-') ++c;
      if(*c=='\0')
        kernel->flag=2;
      else
        {
          *c='\0';
          dstr=c+1;
          if( (dstr[1]!='d' && dstr[1]!='D') || dstr[2]!='\0')
            error(EXIT_FAILURE, 0, "bad formatting in '--kernel' "
                  "dimensionality. The dimensionality suffix must be either "
                  "2d, 3d (not case sensitive). You have given '%s'", dstr);
          switch(dstr[0])
            {
            case '2': kernel->flag=2; break;
            case '3': kernel->flag=3; break;
            default:
              error(EXIT_FAILURE, 0, "only 2 or 3 dimensional kernels can "
                    "currently be built, you have asked for a %c "
                    "dimensional kernel", dstr[0]);
            }
        }


      /* Write the profile type code into 'kernel->status'. If it starts
         with a digit, then the user might have given the code of the
         profile directly. In that case, parse the number. Otherwise,
         let 'ui_profile_name_read' find the value. */
      if( isdigit(*profile) )
        {
          profcode=strtol(profile, &tailptr, 0);
          if(*tailptr!='\0')
            error_at_line(EXIT_FAILURE, 0, filename, lineno, "'%s' "
                          "couldn't be read as a profile code", profile);
          if(profcode<=0 || profcode>=PROFILE_MAXIMUM_CODE)
            error_at_line(EXIT_FAILURE, 0, filename, lineno, "'%s' "
                          "isn't a valid profile code. Please run with "
                          "'--help' and see the acceptable codes in "
                          "explanation of the '--fcol' option", profile);
          kernel->status=profcode;
        }
      else
        kernel->status=ui_profile_name_read(profile, 0);


      /* Make sure the number of parameters conforms with the profile. */
      switch(kernel->status)
        {

        case PROFILE_SERSIC:        need = kernel->flag==2 ? 3 : 4;  break;
        case PROFILE_MOFFAT:        need = kernel->flag==2 ? 3 : 4;  break;
        case PROFILE_GAUSSIAN:      need = kernel->flag==2 ? 2 : 3;  break;
        case PROFILE_POINT:         need = 0;                        break;
        case PROFILE_FLAT:          need = kernel->flag==2 ? 1 : 2;  break;
        case PROFILE_CIRCUMFERENCE: need = kernel->flag==2 ? 1 : 2;  break;
        case PROFILE_DISTANCE:      need = kernel->flag==2 ? 1 : 2;  break;
        default:
          error_at_line(EXIT_FAILURE, 0, filename, lineno, "%s: a bug! "
                        "Please contact us at %s to correct the issue. "
                        "Profile code %d is not recognized", __func__,
                        PACKAGE_BUGREPORT, kernel->status);
        }


      /* Make sure the number of parameters given are the same number that
         are needed. */
      if( kernel->size != need )
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "as a %uD kernel, "
                      "a '%s' profile needs %zu parameters, but %zu "
                      "parameter%s given to '--kernel'", kernel->flag,
                      ui_profile_name_write(kernel->status), need,
                      kernel->size, kernel->size>1?"s are":" is");


      /* Our job is done, return NULL. */
      return NULL;
    }
}





/* Parse the mode to interpret the given coordinates. */
void *
ui_parse_coordinate_mode(struct argp_option *option, char *arg,
                         char *filename, size_t lineno, void *junk)
{
  char *outstr;

  /* We want to print the stored values. */
  if(lineno==-1)
    {
      gal_checkset_allocate_copy( *(uint8_t *)(option->value)==MKPROF_MODE_IMG
                                  ? "img" : "wcs", &outstr );
      return outstr;
    }
  else
    {
      if(!strcmp(arg, "img"))
        *(uint8_t *)(option->value)=MKPROF_MODE_IMG;
      else if (!strcmp(arg, "wcs"))
        *(uint8_t *)(option->value)=MKPROF_MODE_WCS;
      else
        error_at_line(EXIT_FAILURE, 0, filename, lineno, "'%s' (value to "
                      "'--mode') not recognized as a coordinate standard "
                      "mode. Recognized values are 'img' and 'wcs'. This "
                      "option is necessary to identify the nature of your "
                      "input coordinates", arg);
      return NULL;
    }
}



















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
/* Read and check ONLY the options. When arguments are involved, do the
   check in 'ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct mkprofparams *p)
{
  size_t i;

  /* When a no-merged image is to be created, type is necessary. */
  if( p->cp.type==GAL_TYPE_INVALID && p->nomerged==0)
    error(EXIT_FAILURE, 0, "an output type '--type' is necessary when a "
          "merged image is to be built.");

  /* Check if one of the coordinate columns has been given, the other is
     also given. To simplify the job, we use the fact that conditions in C
     return either a 0 (when failed) and 1 (when successful). Note that if
     neighter coordinates are specified there is no problem, the user might
     have input the other coordinate standard. We'll also check for that
     after this.*/
  if(p->kernel==NULL)
    {
      if(p->mode==0)
        error(EXIT_FAILURE, 0, "the '--mode' option is necessary when "
              "building profiles from a catalog. It can take two values: "
              "'img' or 'wcs' which specify how to interpret the "
              "coordinate columns");
    }

  /* The zeropoint magnitude is only necessary when 'mcolisbrightness' is
     not called.  */
  if( p->mcolisbrightness==0 && isnan(p->zeropoint) )
    error(EXIT_FAILURE, 0, "no zeropoint magnitude given. A zeropoint "
          "magnitude is necessary when '--mcolisbrightness' is not called "
          "(i.e., when the contents of '--mcol' must be interpretted as "
          "a magnitude, not brightness).");

  /* Make sure no zero value is given for the '--mergedsize' option (only
     when it is necessary). */
  if(p->dsize && p->backname==NULL)
    for(i=0;p->dsize[i]!=GAL_BLANK_SIZE_T;++i)
      if(p->dsize[i]==0)
        error(EXIT_FAILURE, 0, "values to '--mergedsize' option must not "
              "be zero");
}





/* Sanity check on options AND arguments. If only option values are to be
   checked, use 'ui_read_check_only_options'. */
static void
ui_check_options_and_arguments(struct mkprofparams *p)
{
  int d0f1;
  char *tmpname;

  /* If no kernel is given, make sure an input catalog is given, and if it
     is FITS, that the HDU is also provided. When a kernel option, we will
     set a fiducial catalog name called 'kernel.txt' to automatic output
     filename generation. */
  if(p->kernel)
    {
      if(p->catname)
        error(EXIT_FAILURE, 0, "'--kernel' cannot be called with an input "
              "catalog ('%s'). The parameters necessary to build a single "
              "kernel output should be given to '--kernel', not in a "
              "catalog", p->catname);
      p->catname="kernel.option";
    }
  else
    {
      if(p->catname)
        {
          if( gal_fits_name_is_fits(p->catname) && p->cp.hdu==NULL)
            error(EXIT_FAILURE, 0, "no 'hdu' specified for the input FITS "
                  "table '%s', to ", p->catname);
        }
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
      tmpname=gal_checkset_automatic_output(&p->cp,
                                            ( p->catname
                                              ? p->catname
                                              : "makeprofiles" ), ".fits");
      p->mergedimgname=gal_checkset_malloc_cat(p->outdir, tmpname);
      free(tmpname);
    }
  p->basename=gal_checkset_not_dir_part(p->mergedimgname);


  /* If a merged image is requested (or '--kernel' the option is called),
     then delete the final filename if it exists. */
  if(p->nomerged==0 && p->kernel)
    gal_checkset_writable_remove(p->mergedimgname, p->cp.keep,
                                 p->cp.dontdelete);
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_read_cols_2d(struct mkprofparams *p)
{
  int checkblank;
  size_t i, counter=0;
  char *colname=NULL, **strarr;
  gal_data_t *cols, *tmp, *corrtype=NULL;
  gal_list_str_t *ccol, *lines, *colstrs=NULL;

  /* The coordinate columns are a linked list of strings. */
  ccol=p->ccol;
  for(i=0; i<p->ndim; ++i)
    {
      gal_list_str_add(&colstrs, ccol->v, 0);
      ccol=ccol->next;
    }

  /* Add the rest of the columns in a specific order. */
  gal_list_str_add(&colstrs, p->fcol, 0);
  gal_list_str_add(&colstrs, p->rcol, 0);
  gal_list_str_add(&colstrs, p->ncol, 0);
  gal_list_str_add(&colstrs, p->pcol, 0);
  gal_list_str_add(&colstrs, p->qcol, 0);
  gal_list_str_add(&colstrs, p->mcol, 0);
  gal_list_str_add(&colstrs, p->tcol, 0);

  /* Reverse the order to make the column orders correspond to how we added
     them here and avoid possible bugs. */
  gal_list_str_reverse(&colstrs);

  /* Read the desired columns from the file. */
  lines=gal_options_check_stdin(p->catname, p->cp.stdintimeout, "input");
  cols=gal_table_read(p->catname, p->cp.hdu, lines, colstrs, p->cp.searchin,
                      p->cp.ignorecase, p->cp.minmapsize, p->cp.quietmmap,
                      NULL);
  gal_list_str_free(lines, 1);

  /* The name of the input catalog is only for informative steps from now
     on (we won't be dealing with the actual file any more). So if the
     standard input was used (therefore 'catname==NULL', set it to
     'stdin'). */
  if(p->catname==NULL)
    gal_checkset_allocate_copy("standard-input", &p->catname);

  /* Set the number of objects. */
  p->num = cols ? cols->size : 0;

  /* Put each column's data in the respective internal array. */
  while(cols!=NULL)
    {
      /* Pop out the top column. */
      tmp=gal_list_data_pop(&cols);

      /* By default check if the column has blank values, but it can be
         turned off for some columns. */
      checkblank=1;

      /* See which column we are currently reading. */
      switch(++counter)
        {
        case 1:
        case 2:
          colname = ( counter==1
                      ? "first coordinate column ('--coordcol')"
                      : "second coordinate column ('--coordcol')" );
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT64);
          switch(counter)
            {
            case 1: p->x=corrtype->array; break;
            case 2: p->y=corrtype->array; break;
            }
          break;


        case 3:
          if(tmp->type==GAL_TYPE_STRING)
            {
              p->f=gal_pointer_allocate(GAL_TYPE_UINT8, p->num, 0,
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
              colname="profile function code ('fcol')";
              corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_UINT8);
              p->f=corrtype->array;

              /* Check if they are in the correct range. */
              for(i=0;i<p->num;++i)
                if(p->f[i]<=PROFILE_INVALID || p->f[i]>=PROFILE_MAXIMUM_CODE)
                  error(EXIT_FAILURE, 0, "%s: row %zu, the function "
                        "code is %u. It should be >%d and <%d. Please run "
                        "again with '--help' and check the acceptable "
                        "codes.\n\nAlternatively, you can use alphabetic "
                        "strings to specify the profile functions, see the "
                        "explanations under 'fcol' from the command "
                        "below (press the 'SPACE' key to go down, and the "
                        "'q' to return back to the command-line):\n\n"
                        "    $ info %s\n", p->catname, i+1, p->f[i],
                        PROFILE_INVALID, PROFILE_MAXIMUM_CODE, PROGRAM_EXEC);
            }
          break;


        case 4:
          colname="radius ('rcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->r=corrtype->array;

          /* Check if there is no negative or zero-radius profile. */
          for(i=0;i<p->num;++i)
            if(p->f[i]!=PROFILE_POINT && p->r[i]<=0.0f)
              error(EXIT_FAILURE, 0, "%s: row %zu, the radius value %g is "
                    "not acceptable for a '%s' profile. It has to be larger "
                    "than 0", p->catname, i+1, p->r[i],
                    ui_profile_name_write(p->f[i]));
          break;


        case 5:
          colname="index ('ncol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->n=corrtype->array;
          break;


        case 6:
          colname="position angle ('pcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->p1=corrtype->array;
          break;


        case 7:
          colname="axis ratio ('qcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->q1=corrtype->array;

          /* Check if there is no negative or >1.0f axis ratio. */
          for(i=0;i<p->num;++i)
            if( p->f[i]!=PROFILE_POINT && (p->q1[i]<=0.0f || p->q1[i]>1.0f) )
              error(EXIT_FAILURE, 0, "%s: row %zu, the axis ratio value %g "
                    "is not acceptable for a '%s' profile. It has to be >0 "
                    "and <=1", p->catname, i+1, p->q1[i],
                    ui_profile_name_write(p->f[i]));
          break;


        case 8:
          colname="magnitude ('mcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->m=corrtype->array;
          checkblank=0;          /* Magnitude can be NaN: to mask regions. */
          break;


        case 9:
          colname="truncation ('tcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->t=corrtype->array;

          /* Check if there is no negative or zero truncation radius. */
          for(i=0;i<p->num;++i)
            if(p->f[i]!=PROFILE_POINT && p->t[i]<=0.0f)
              error(EXIT_FAILURE, 0, "%s: row %zu, the truncation radius "
                    "value %g is not acceptable for a '%s' profile. It has "
                    "to be larger than 0", p->catname, i+1, p->t[i],
                    ui_profile_name_write(p->f[i]));
          break;


        /* If the index isn't recognized, then it is larger, showing that
           there was more than one match for the given criteria */
        default:
          gal_tableintern_error_col_selection(p->catname, p->cp.hdu, "too "
                                              "many columns were selected "
                                              "by the given values to the "
                                              "options ending in 'col'.");
        }

      /* Sanity check and clean up.  Note that it might happen that the
         input structure is already freed. In that case, 'corrtype' will be
         NULL. */
      if(corrtype)
        {
          /* Make sure there are no blank values in this column. */
          if( checkblank && gal_blank_present(corrtype, 1) )
            error(EXIT_FAILURE, 0, "%s column has blank values. "
                  "Input columns cannot contain blank values", colname);

          /* Free the unnecessary sturcture information. The correct-type
             ('corrtype') data structure's array is necessary for later
             steps, so its pointer has been copied in the main program's
             structure. Hence, we should set the structure's pointer to
             NULL so the important data isn't freed.*/
          corrtype->array=NULL;
          gal_data_free(corrtype);
        }
    }

  /* Multi-column sanity checks. */
  counter=0;
  if( !p->cp.quiet && (p->mforflatpix || p->mcolisbrightness) )
    for(i=0;i<p->num;++i)
      if( p->m[i]==0.0 && ( p->f[i]==PROFILE_POINT
                            || p->f[i]==PROFILE_FLAT
                            || p->f[i]==PROFILE_CIRCUMFERENCE ) )
        {
          error(0, 0, "WARNING: atleast one single-valued profile (point, "
                "flat, or circumference profiles) has a magnitude column "
                "value of 0.0 while '--mforflatpix' or "
                "'--mcolforbrightness' have also been given. In such cases "
                "the profile's pixels will have a value of zero and thus "
                "they will not be identifiable from the zero-valued "
                "background. If this behavior is intended, this warning "
                "can be suppressed with the '--quiet' (or '-q') option.\n");
          break;
        }
}





/* Read the columns for a 3D profile. */
static void
ui_read_cols_3d(struct mkprofparams *p)
{
  int checkblank;
  size_t i, counter=0;
  char *colname=NULL, **strarr;
  gal_data_t *cols, *tmp, *corrtype=NULL;
  gal_list_str_t *lines, *ccol, *colstrs=NULL;

  /* The 3D-specific columns are not mandatory in 'args.h', so we need to
     check here if they are given or not before starting to read them. */
  if(p->p2col==NULL || p->p3col==NULL || p->q2col==NULL)
    error(EXIT_FAILURE, 0, "at least one of '--p2col', '--p3col', or "
          "'--q2col' have not been identified. When building a 3D profile, "
          "these three columns are also mandatory");

  /* The coordinate columns are a linked list of strings. */
  ccol=p->ccol;
  for(i=0; i<p->ndim; ++i)
    {
      gal_list_str_add(&colstrs, ccol->v, 0);
      ccol=ccol->next;
    }

  /* Put the columns a specific order to read. */
  gal_list_str_add(&colstrs, p->fcol,  0);
  gal_list_str_add(&colstrs, p->rcol,  0);
  gal_list_str_add(&colstrs, p->ncol,  0);
  gal_list_str_add(&colstrs, p->pcol,  0);
  gal_list_str_add(&colstrs, p->p2col, 0);
  gal_list_str_add(&colstrs, p->p3col, 0);
  gal_list_str_add(&colstrs, p->qcol,  0);
  gal_list_str_add(&colstrs, p->q2col, 0);
  gal_list_str_add(&colstrs, p->mcol,  0);
  gal_list_str_add(&colstrs, p->tcol,  0);

  /* Reverse the order to make the column orders correspond to how we added
     them here and avoid possible bugs. */
  gal_list_str_reverse(&colstrs);

  /* Read the desired columns from the file. */
  lines=gal_options_check_stdin(p->catname, p->cp.stdintimeout, "input");
  cols=gal_table_read(p->catname, p->cp.hdu, lines, colstrs, p->cp.searchin,
                      p->cp.ignorecase, p->cp.minmapsize, p->cp.quietmmap,
                      NULL);
  gal_list_str_free(lines, 1);

  /* Set the number of objects. */
  p->num=cols->size;

  /* Put each column's data in the respective internal array. */
  while(cols!=NULL)
    {
      /* Pop out the top column. */
      tmp=gal_list_data_pop(&cols);

      /* By default check if the column has blank values, but it can be
         turned off for some columns. */
      checkblank=1;

      /* See which column we are currently reading. */
      switch(++counter)
        {
        case 1:
        case 2:
        case 3:
          colname = ( counter==1
                      ? "first coordinate column ('--coordcol')"
                      : ( counter==2
                          ? "second coordinate column ('--coordcol')"
                          : "third coordinate column ('--coordcol')" ) );
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT64);
          switch(counter)
            {
            case 1: p->x=corrtype->array; break;
            case 2: p->y=corrtype->array; break;
            case 3: p->z=corrtype->array; break;
            }
          break;

        case 4:
          if(tmp->type==GAL_TYPE_STRING)
            {
              p->f=gal_pointer_allocate(GAL_TYPE_UINT8, p->num, 0,
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
              colname="profile function code ('fcol')";
              corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_UINT8);
              p->f=corrtype->array;

              /* Check if they are in the correct range. */
              for(i=0;i<p->num;++i)
                if(p->f[i]<=PROFILE_INVALID || p->f[i]>=PROFILE_MAXIMUM_CODE)
                  error(EXIT_FAILURE, 0, "%s: row %zu, the function "
                        "code is %u. It should be >%d and <%d. Please run "
                        "again with '--help' and check the acceptable "
                        "codes.\n\nAlternatively, you can use alphabetic "
                        "strings to specify the profile functions, see the "
                        "explanations under 'fcol' from the command "
                        "below (press the 'SPACE' key to go down, and the "
                        "'q' to return back to the command-line):\n\n"
                        "    $ info %s\n", p->catname, i+1, p->f[i],
                        PROFILE_INVALID, PROFILE_MAXIMUM_CODE, PROGRAM_EXEC);
            }
          break;

        case 5:
          colname="radius ('rcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->r=corrtype->array;

          /* Check if there is no negative or zero-radius profile. */
          for(i=0;i<p->num;++i)
            if(p->f[i]!=PROFILE_POINT && p->r[i]<=0.0f)
              error(EXIT_FAILURE, 0, "%s: row %zu, the radius value %g is "
                    "not acceptable for a '%s' profile. It has to be larger "
                    "than 0", p->catname, i+1, p->r[i],
                    ui_profile_name_write(p->f[i]));
          break;

        case 6:
          colname="index ('ncol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->n=corrtype->array;
          break;

        case 7:
          colname="first euler angle ('pcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->p1=corrtype->array;
          break;

        case 8:
          colname="second euler angle ('p2col')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->p2=corrtype->array;
          break;

        case 9:
          colname="third euler angle ('p3col')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->p3=corrtype->array;
          break;

        case 10:
          colname="axis ratio 1 ('qcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->q1=corrtype->array;

          /* Check if there is no negative or >1.0f axis ratio. */
          for(i=0;i<p->num;++i)
            if( p->f[i]!=PROFILE_POINT && (p->q1[i]<=0.0f || p->q1[i]>1.0f) )
              error(EXIT_FAILURE, 0, "%s: row %zu, the first axis ratio "
                    "value %g is not acceptable for a '%s' profile. It has "
                    "to be >0 and <=1", p->catname, i+1, p->q1[i],
                    ui_profile_name_write(p->f[i]));
          break;

        case 11:
          colname="axis ratio 2 ('q2col')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->q2=corrtype->array;

          /* Check if there is no negative or >1.0f axis ratio. */
          for(i=0;i<p->num;++i)
            if( p->f[i]!=PROFILE_POINT && (p->q2[i]<=0.0f || p->q2[i]>1.0f) )
              error(EXIT_FAILURE, 0, "%s: row %zu, the second axis ratio "
                    "value %g is not acceptable for a '%s' profile. It has "
                    "to be >0 and <=1", p->catname, i+1, p->q2[i],
                    ui_profile_name_write(p->f[i]));
          break;

        case 12:
          colname="magnitude ('mcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->m=corrtype->array;
          checkblank=0;          /* Magnitude can be NaN: to mask regions. */
          break;

        case 13:
          colname="truncation ('tcol')";
          corrtype=gal_data_copy_to_new_type_free(tmp, GAL_TYPE_FLOAT32);
          p->t=corrtype->array;

          /* Check if there is no negative or zero truncation radius. */
          for(i=0;i<p->num;++i)
            if(p->f[i]!=PROFILE_POINT && p->t[i]<=0.0f)
              error(EXIT_FAILURE, 0, "%s: row %zu, the truncation radius "
                    "value %g is not acceptable for a '%s' profile. It has "
                    "to be larger than 0", p->catname, i+1, p->t[i],
                    ui_profile_name_write(p->f[i]));
          break;

        /* If the index isn't recognized, then it is larger, showing that
           there was more than one match for the given criteria */
        default:
          gal_tableintern_error_col_selection(p->catname, p->cp.hdu, "too "
                                              "many columns were selected "
                                              "by the given values to the "
                                              "options ending in 'col'.");
        }

      /* Sanity check and clean up.  Note that it might happen that the
         input structure is already freed. In that case, 'corrtype' will be
         NULL. */
      if(corrtype)
        {
          /* Make sure there are no blank values in this column. */
          if( checkblank && gal_blank_present(corrtype, 1) )
            error(EXIT_FAILURE, 0, "%s column has blank values. "
                  "Input columns cannot contain blank values", colname);

          /* Free the unnecessary sturcture information. The correct-type
             ('corrtype') data structure's array is necessary for later
             steps, so its pointer has been copied in the main program's
             structure. Hence, we should set the structure's pointer to
             NULL so the important data isn't freed.*/
          corrtype->array=NULL;
          gal_data_free(corrtype);
        }
    }
}





/* It is possible to define the internal catalog through a catalog or the
   '--kernel' option. This function will do the job. */
static void
ui_prepare_columns(struct mkprofparams *p)
{
  double *karr;
  float r, n, t, q2;

  /* If the kernel option was called, then we need to build a series of
     single element columns to create an internal catalog. */
  if(p->kernel)
    {
      /* Number of profiles to be built. */
      p->num=1;

      /* Allocate the necessary columns. */
      p->x  = gal_pointer_allocate(GAL_TYPE_FLOAT64, 1, 1, __func__, "p->x");
      p->y  = gal_pointer_allocate(GAL_TYPE_FLOAT64, 1, 1, __func__, "p->y");
      p->f  = gal_pointer_allocate(GAL_TYPE_UINT8,   1, 1, __func__, "p->f");
      p->r  = gal_pointer_allocate(GAL_TYPE_FLOAT32, 1, 1, __func__, "p->r");
      p->n  = gal_pointer_allocate(GAL_TYPE_FLOAT32, 1, 1, __func__, "p->n");
      p->p1 = gal_pointer_allocate(GAL_TYPE_FLOAT32, 1, 1, __func__, "p->p1");
      p->q1 = gal_pointer_allocate(GAL_TYPE_FLOAT32, 1, 1, __func__, "p->q1");
      p->m  = gal_pointer_allocate(GAL_TYPE_FLOAT32, 1, 1, __func__, "p->m");
      p->t  = gal_pointer_allocate(GAL_TYPE_FLOAT32, 1, 1, __func__, "p->t");
      if(p->ndim==3)
        {
          p->z =gal_pointer_allocate(GAL_TYPE_FLOAT64, 1, 1, __func__, "p->z");
          p->p2=gal_pointer_allocate(GAL_TYPE_FLOAT32, 1, 1, __func__, "p->p2");
          p->p3=gal_pointer_allocate(GAL_TYPE_FLOAT32, 1, 1, __func__, "p->p3");
          p->q2=gal_pointer_allocate(GAL_TYPE_FLOAT32, 1, 1, __func__, "p->q2");
        }

      /* For profiles that need a different number of input values. Note
         that when a profile doesn't need a value, it will be ignored. */
      karr=p->kernel->array;
      if(p->kernel->size)
        {
          r = karr[0];
          n = p->kernel->size==2 ? 0.0f : karr[1];
          t = ( p->ndim==2
                ? p->kernel->size==1 ? 1.0f : karr[ p->kernel->size - 1 ]
                : p->kernel->size==1 ? 1.0f : karr[ p->kernel->size - 2 ] );
        }
      else r=n=t=0.0f;

      /* Fill the allocated spaces. */
      p->x[0]  = 0.0f;
      p->y[0]  = 0.0f;
      p->f[0]  = p->kernel->status;
      p->r[0]  = r;
      p->n[0]  = n;
      p->p1[0] = 0.0f;
      p->q1[0] = 1.0f;
      p->m[0]  = 0.0f;
      p->t[0]  = t;
      if(p->ndim==3)
        {
          /* Parameters for any case. */
          p->z[0] = 0.0f;
          q2      = p->kernel->size ? karr[ p->kernel->size - 1 ] : 0.0f;

          /* 3rd-dim axis ratio > 1: Set the major axis in the direction of
             the 3rd dimension (90 degree rotation for all three
             rotations). Also set the two axis ratios to the inverse of the
             requested value. */
          if(q2>1.0)
            {
              p->q1[0] = p->q2[0] = 1/q2;
              p->p1[0] = p->p2[0] = p->p3[0] = 90.0;
            }

          /* 3rd-dim axis ratio <=1: No extra rotation is necessary and
             'q2'can simply be put in the respective column. */
          else
            {
              p->q2[0] = q2;
              p->p2[0] = p->p3[0] = 0.0;
            }
        }
    }
  else
    {
      /* Make sure the number of coordinate columns and number of
         dimensions in outputs are the same. There is no problem if it is
         more than 'ndim'. In that case, the last values (possibly in
         configuration files) will be ignored. */
      if( gal_list_str_number(p->ccol) < p->ndim )
        error(EXIT_FAILURE, 0, "%zu coordinate columns (calls to "
              "'--coordcol') given but output has %zu dimensions",
              gal_list_str_number(p->ccol), p->ndim);

      /* Call the respective function. */
      switch(p->ndim)
        {
        case 2: ui_read_cols_2d(p);   break;
        case 3: ui_read_cols_3d(p);   break;
        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                "resolve the issue. %zu not recognized for 'p->ndim'",
                __func__, PACKAGE_BUGREPORT, p->ndim);
        }
    }
}





/* To keep things clean, we'll do the WCS sanity checks in this small
   function. If everything is ok, this function will return 0 (so an if
   condition won't be executed). If any of the necessary inputs aren't
   given, it will return 1. */
static int
ui_wcs_sanity_check(struct mkprofparams *p)
{
  size_t ndim=p->ndim;

  if(p->crpix)
    {
      if(p->crpix->size!=ndim)
        error(EXIT_FAILURE, 0, "%zu values given to '--crpix'. This must be "
              "the same as the output dimension (%zu)", p->crpix->size, ndim);
      return 0;
    }
  else return 1;

  if(p->crval)
    {
      if(p->crval->size!=ndim)
        error(EXIT_FAILURE, 0, "%zu values given to '--crval'. This must be "
              "the same as the output dimension (%zu)", p->crval->size, ndim);
      return 0;
    }
  else return 1;

  if(p->cdelt)
    {
      if(p->cdelt->size!=ndim)
        error(EXIT_FAILURE, 0, "%zu values given to '--cdelt'. This must be "
              "the same as the output dimension (%zu)", p->cdelt->size, ndim);
      return 0;
    }
  else return 1;

  if(p->pc)
    {
      if(p->pc->size!=ndim*ndim)
        error(EXIT_FAILURE, 0, "%zu values given to '--pc'. This must be "
              "the square as the output dimension (%zu)", p->pc->size,
              ndim*ndim);
      return 0;
    }
  else return 1;

  if(p->cunit)
    {
      if(p->cunit->size!=ndim)
        error(EXIT_FAILURE, 0, "%zu values given to '--cunit'. This must be "
              "the same as the output dimension (%zu)", p->cunit->size, ndim);
      return 0;
    }
  else return 1;

  if(p->ctype)
    {
      if(p->ctype->size!=ndim)
        error(EXIT_FAILURE, 0, "%zu values given to '--ctype'. This must be "
              "the same as the output dimension (%zu)", p->ctype->size, ndim);
      return 0;
    }
  else return 1;
}





static void
ui_prepare_wcs(struct mkprofparams *p)
{
  int status;
  struct wcsprm *wcs;
  char **cunit, **ctype;
  size_t i, ndim=p->ndim;
  double *crpix, *crval, *cdelt, *pc;


  /* Check and initialize the WCS information. If any of the necessary WCS
     parameters are missing, then don't build any WCS. */
  if( ui_wcs_sanity_check(p) ) return;
  crpix = p->crpix->array;
  crval = p->crval->array;
  cdelt = p->cdelt->array;
  pc    = p->pc->array;
  cunit = p->cunit->array;
  ctype = p->ctype->array;


  /* Allocate the memory necessary for the wcsprm structure. */
  errno=0;
  wcs=p->wcs=malloc(sizeof *wcs);
  if(wcs==NULL)
    error(EXIT_FAILURE, errno, "%zu for wcs in preparewcs", sizeof *wcs);


  /* Initialize the structure (allocate all its internal arrays). */
  wcs->flag=-1;
  if( (status=wcsini(1, ndim, wcs)) )
    error(EXIT_FAILURE, 0, "wcsini error %d: %s",
          status, wcs_errmsg[status]);


  /* Fill in all the important WCS structure parameters. */
  wcs->altlin   = 0x1;
  wcs->equinox  = 2000.0f;
  for(i=0;i<ndim;++i)
    {
      /* IMPORTANT: At this point, we don't want the WCS to be over-sampled
         because if the user has given RA and Dec for the profiles, they
         need to be converted to non-oversampled and shifted image
         coordinates. After the conversion (in 'ui_finalize_coordinates')
         we are going to correct for the oversampling in the WCS.*/
      wcs->crpix[i] = crpix[i];
      wcs->crval[i] = crval[i];
      wcs->cdelt[i] = cdelt[i];
      strcpy(wcs->cunit[i], cunit[i]);
      strcpy(wcs->ctype[i], ctype[i]);
    }
  for(i=0;i<ndim*ndim;++i) wcs->pc[i]=pc[i];

  /* Set up the wcs structure with the constants defined above. */
  status=wcsset(wcs);
  if(status)
    error(EXIT_FAILURE, 0, "wcsset error %d: %s", status,
          wcs_errmsg[status]);
}





static void
ui_prepare_canvas(struct mkprofparams *p)
{
  float *f, *ff;
  long width[3]={1,1,1};
  size_t tndim, *tdsize;
  int status=0, setshift=0;
  double truncr, semiaxes[3], euler_deg[3];
  size_t i, nshift=0, *dsize=NULL, ndim_counter;

  /* If a background image is specified, then use that as the output
     image to build the profiles over. */
  if(p->backname)
    {
      /* Read in the background image and its coordinates, note that when
         no merged image is desired, we just need the WCS information of
         the background image and the number of its dimensions. So
         'ndim==0' and what 'dsize' points to is irrelevant. */
      tdsize=gal_fits_img_info_dim(p->backname, p->backhdu, &tndim);
      p->wcs=gal_wcs_read(p->backname, p->backhdu, 0, 0, &p->nwcs);
      tndim=gal_dimension_remove_extra(tndim, tdsize, p->wcs);
      free(tdsize);
      if(p->nomerged==0)
        {
          /* If p->dsize was given as an option, free it. */
          if( p->dsize ) free(p->dsize);

          /* Write the size of the background image into 'dsize'. */
          p->dsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->ndim, 0,
                                        __func__, "p->dsize");
          for(i=0;i<p->ndim;++i) p->dsize[i] = p->out->dsize[i];

          /* Set all pixels to zero if the user wanted a clear canvas. */
          if(p->clearcanvas)
            {ff=(f=p->out->array)+p->out->size; do *f++=0.0f; while(f<ff);}
        }

      /* When a background image is specified, oversample must be 1 and
         there is no shifts. */
      p->oversample=1;
      if(p->shift) free(p->shift);
      p->shift=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->ndim, 1, __func__,
                                    "p->shift (1)");
    }
  else
    {
      /* If any of the shift elements are zero, the others should be too!*/
      if(p->shift && p->shift[0] && p->shift[1])
        {
          /* Multiply the shift by the over-sample. */
          for(i=0;p->shift[i]!=GAL_BLANK_SIZE_T;++i)
            {
              ++nshift;
              p->shift[i] *= p->oversample;
            }

          /* Make sure it has the same number of elements as naxis. */
          if(p->ndim!=nshift)
            error(EXIT_FAILURE, 0, "%zu and %zu elements given to '--ndim' "
                  "and '--shift' respectively. These two numbers must be the "
                  "same", p->ndim, nshift);
        }
      else
        {
          /* 'prepforconv' is only valid when xshift and yshift are both
             zero. Also, a PSF profile should exist in the image. */
          if(p->prepforconv)
            {
              /* Check if there is at least one Moffat or Gaussian profile. */
              for(i=0;i<p->num;++i)
                if( oneprofile_ispsf(p->f[i]) )
                  {
                    /* Calculate the size of the box holding the PSF. Note:

                       - For the Moffat and Gaussian profiles, the radius
                       column is actually the FWHM which is actually the
                       diameter, not radius. So we have to divide it by
                       half.

                       - encloseellipse outputs the total width, we only want
                       half of it for the shift. */
                    setshift=1;
                    truncr = p->tunitinp ? p->t[i] : p->t[i] * p->r[i]/2;
                    if(p->ndim==2)
                      gal_box_bound_ellipse(truncr, p->q1[i]*truncr,
                                            p->p1[i], width);
                    else
                      {
                        euler_deg[0] = p->p1[i];
                        euler_deg[1] = p->p2[i];
                        euler_deg[2] = p->p3[i];
                        semiaxes[0]  = truncr;
                        semiaxes[1]  = truncr * p->q1[i];
                        semiaxes[2]  = truncr * p->q2[i];
                        gal_box_bound_ellipsoid( semiaxes, euler_deg, width);
                      }
                  }

              /* Either set the shifts to zero or to the values set from
                 the PSF. Note that the user might have given any number of
                 shifts (from zero). So, we'll just free it and reset
                 it. */
              if(p->shift) free(p->shift);
              p->shift=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->ndim, 1,
                                            __func__, "p->shift (2)");
              if(setshift)
                {
                  p->shift[0]  = (width[0]/2)*p->oversample;
                  p->shift[1]  = (width[1]/2)*p->oversample;
                  if(p->ndim==3) p->shift[2] = (width[2]/2)*p->oversample;
                }
            }
        }

      /* If shift has not been set until now, set it. */
      if(p->shift==NULL)
        p->shift=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->ndim, 1,
                                      __func__, "p->shift (3)");

      /* Prepare the sizes of the final merged image (if it is to be
         made). Note that even if we don't want a merged image, we still
         need its WCS structure. */
      if(p->nomerged==0)
        {
          ndim_counter=0;
          for(i=0;p->dsize[i]!=GAL_BLANK_SIZE_T;++i)
            {
              /* Count the number of dimensions. */
              ++ndim_counter;

              /* Correct dsize. */
              p->dsize[i] = (p->dsize[i]*p->oversample) + (2*p->shift[i]);
            }
          dsize = p->dsize;

          /* Make the output structure. */
          p->out=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, ndim_counter, dsize,
                                NULL, 1, p->cp.minmapsize, p->cp.quietmmap,
                                NULL, NULL, NULL);
        }
    }


  /* Make the WCS structure of the output data structure if it has not
     been set yet. */
  if(p->wcs==NULL)
    ui_prepare_wcs(p);


  /* Set the name, comments and units of the final merged output. */
  if(p->out)
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
  if(p->individual && p->wcs)
    {
      status=wcshdo(WCSHDO_safe, p->wcs, &p->wcsnkeyrec, &p->wcsheader);
      if(status)
        error(EXIT_FAILURE, 0, "wcshdo error %d: %s", status,
              wcs_errmsg[status]);
    }
}





static void
ui_finalize_coordinates(struct mkprofparams *p)
{
  void *arr=NULL;
  size_t i=0, ndim=p->ndim;
  uint8_t os=p->oversample;
  gal_data_t *tmp, *coords=NULL;
  double *cdelt=p->wcs->cdelt, *crpix=p->wcs->crpix;

  /* When the user specified RA and Dec columns, the respective values
     where stored in the 'p->x' and 'p->y' arrays. So before proceeding, we
     need to change them into actual image coordinates. */
  if(p->mode==MKPROF_MODE_WCS)
    {
      /* Make list of coordinates for input of 'gal_wcs_world_to_img'. */
      for(i=0;i<ndim;++i)
        {
          /* Set the array pointer. Note that we read the WCS columns into
         the 'p->x', 'p->y' and 'p->z' arrays temporarily before. Here, we
         will convert them to image coordinates in place. */
          switch(i)
            {
            /* Note that the linked list gets filled in a first-in-last-out
               order, so the last column added should be the first WCS
               dimension. */
            case 0: arr = ndim==2 ? p->y : p->z;   break;
            case 1: arr = ndim==2 ? p->x : p->y;   break;
            case 2: arr = p->x;                    break;
            default:
              error(EXIT_FAILURE, 0, "conversion from WCS to image "
                    "coordinates is not supported for %zu-dimensional "
                    "datasets", ndim);
            }

          /* Allocate the list of coordinates. */
          gal_list_data_add_alloc(&coords, arr, GAL_TYPE_FLOAT64, 1, &p->num,
                                  NULL, 0, -1, 1, NULL, NULL, NULL);
        }

      /* Convert the world coordinates to image coordinates (inplace). */
      gal_wcs_world_to_img(coords, p->wcs, 1);


      /* If any conversions created a WCSLIB error, both the outputs will be
         set to NaN. */
      for(i=0;i<p->num;++i)
        if( isnan(p->x[i]) )
          error(EXIT_FAILURE, 0, "catalog row %zu: WCSLIB could not convert "
                "(%f, %f) coordinates into image coordinates", i, p->x[i],
                p->y[i]);

      /* We want the actual arrays of each 'coords' column. So, first we'll
         set all the array elements to NULL, then free it. */
      for(tmp=coords;tmp!=NULL;tmp=tmp->next) tmp->array=NULL;
      gal_list_data_free(coords);
    }

  /* Correct the WCS scale. Note that when the WCS is read from a
     background image, oversample is set to 1. This is done here because
     the conversion of WCS to pixel coordinates needs to be done with the
     non-over-sampled image.*/
  for(i=0;i<p->ndim;++i)
    {
      /* Oversampling has already been applied in 'p->shift'. Also note
         that shift is in the C dimension ordring, while crpix is in FITS
         ordering. */
      crpix[i]  = crpix[i]*os + p->shift[ndim-i-1] - os/2;
      cdelt[i] /= os;
    }

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
                          1, p->cp.minmapsize, p->cp.quietmmap,
                          "INDIV_CREATED", "bool",
                          "If an individual image was made (1) or not (0).");

  /* Fraction of monte-carlo. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_FLOAT32, 1, &p->num, NULL,
                          1, p->cp.minmapsize, p->cp.quietmmap,
                          "FRAC_MONTECARLO", "frac",
                          "Fraction of brightness in Monte-carlo integrated "
                          "pixels.");

  /* Number of monte-carlo. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_UINT64, 1, &p->num, NULL,
                          1, p->cp.minmapsize, p->cp.quietmmap,
                          "NUM_MONTECARLO", "count",
                          "Number of Monte Carlo integrated pixels.");

  /* Magnitude of profile overlap. */
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_FLOAT32, 1, &p->num, NULL,
                          1, p->cp.minmapsize, p->cp.quietmmap, "MAG_OVERLAP",
                          "mag", "Magnitude of profile's overlap with merged "
                          "image.");

  /* Row number in input catalog. */
  name=gal_fits_name_save_as_string(p->catname, p->cp.hdu);
  if( asprintf(&comment, "Row number of profile in %s.", name)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_list_data_add_alloc(&p->log, NULL, GAL_TYPE_UINT64, 1, &p->num, NULL,
                          1, p->cp.minmapsize, p->cp.quietmmap, "INPUT_ROW_NO",
                          "count", comment);
  free(comment);
  free(name);
}





static void
ui_read_ndim(struct mkprofparams *p)
{
  size_t i, *dsize, ndim_counter;

  if(p->kernel)
    {
      /* The kernel's dimensionality is fixed. */
      p->ndim=p->kernel->flag;

      /* Make sure the kernel and background are not given together. */
      if(p->backname)
        error(EXIT_FAILURE, 0, "the '--kernel' and '--background' options "
              "cannot be called together");
    }
  else
    {
      /* Packground image is given. */
      if(p->backname)
        {
          /* Small sanity check. */
          if(p->backhdu==NULL)
            error(EXIT_FAILURE, 0, "no hdu specified for the background "
                  "image %s. Please run again '--backhdu' option",
                  p->backname);

          /* If '--nomerged' is given, we don't actually need to load the
             image, we just need its WCS later. */
          if(p->nomerged)
            {
              /* Get the number of the background image's dimensions. */
              dsize=gal_fits_img_info_dim(p->backname, p->backhdu, &p->ndim);
              p->ndim=gal_dimension_remove_extra(p->ndim, dsize, NULL);
              free(dsize);
            }
          else
            {
              /* Read the image. */
              p->out=gal_array_read_one_ch_to_type(p->backname, p->backhdu,
                                                   NULL, GAL_TYPE_FLOAT32,
                                                   p->cp.minmapsize,
                                                   p->cp.quietmmap);
              p->out->ndim=gal_dimension_remove_extra(p->out->ndim,
                                                      p->out->dsize, NULL);
              p->ndim=p->out->ndim;
            }

          /* Make sure the dimensionality is supported. */
          if(p->ndim!=2 && p->ndim!=3)
            error(EXIT_FAILURE, 0, "%s (hdu %s) has %zu dimensions. Currently "
                  "only 2 or 3 dimensional outputs can be produced",
                  p->backname, p->backhdu, p->ndim);
        }
      else
        {
          /* Get the number of dimensions from the user's options. */
          ndim_counter=0;
          for(i=0;p->dsize[i]!=GAL_BLANK_SIZE_T;++i) ++ndim_counter;
          p->ndim=ndim_counter;

          /* Make sure the dimensionality is supported. */
          if(p->ndim!=2 && p->ndim!=3)
            error(EXIT_FAILURE, 0, "%zu values given to '--mergedsize'. "
                  "Currently only 2 or 3 dimensional outputs can be produced",
                  p->ndim);
        }
    }
}





static void
ui_preparations(struct mkprofparams *p)
{
  /* Set the output dimensionality (necessary to know which columns to
     use). */
  ui_read_ndim(p);

  /* Read in all the columns (necessary for '--prepforconf' when we want to
     build the profiles). */
  ui_prepare_columns(p);

  /* If the kernel option was given, some parameters need to be
     over-written: */
  if(p->kernel)
    {
      /* Set the necessary constants. */
      p->nomerged=1;
      p->psfinimg=0;
      p->individual=1;
      p->ndim=p->kernel->flag;

      /* Set the shift array. */
      p->shift=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->ndim, 1,
                                    __func__, "p->shift");
    }
  else
    ui_prepare_canvas(p);

  /* Read the (possible) RA/Dec inputs into X and Y for the builder.  NOTE:
     It may happen that there are no input columns, in that case, just
     ignore this step.*/
  if(p->wcs && p->num)
    ui_finalize_coordinates(p);

  /* Prepare the random number generator. */
  p->rng=gal_checkset_gsl_rng(p->envseed, &p->rng_name, &p->rng_seed);

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

  printf(PROGRAM_NAME" "PACKAGE_VERSION" started on %s", ctime(&p->rawtime));

  if(p->kernel)
    {
      if( asprintf(&jobname, "Building one %s kernel",
                   ui_profile_name_write(p->kernel->status))<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
    }
  else
    {
      if( asprintf(&jobname, "%zu profile%sread from %s", p->num,
                   p->num>1?"s ":" ", p->catname)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
    }
  gal_timing_report(NULL, jobname, 1);
  free(jobname);

  if(p->backname)
    {
      if(p->nomerged)
        {
          if( asprintf(&jobname, "WCS information read from %s",
                       p->backname)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      else
        {
          if( asprintf(&jobname, "%s is read and will be used as canvas",
                       p->backname)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }
      gal_timing_report(NULL, jobname, 1);
      free(jobname);
    }

  if( asprintf(&jobname, "Random number generator (RNG) type: %s",
               gsl_rng_name(p->rng))<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_timing_report(NULL, jobname, 1);
  free(jobname);

  if( asprintf(&jobname, "Basic RNG seed: %lu", p->rng_seed)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  gal_timing_report(NULL, jobname, 1);
  free(jobname);

  if(p->kernel==NULL)
    {
      if( asprintf(&jobname, "Using %zu threads.", p->cp.numthreads)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(NULL, jobname, 1);
      free(jobname);
    }
}





void
ui_read_check_inputs_setup(int argc, char *argv[], struct mkprofparams *p)
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


  /* Read the configuration files. */
  gal_options_read_config_set(&p->cp);


  /* Read the options into the program's structure, and check them and
     their relations prior to printing. */
  ui_read_check_only_options(p);


  /* Print the option values if asked. Note that this needs to be done
     after the sanity check so un-sane values are not printed in the output
     state. */
  gal_options_print_state(&p->cp);


  /* Prepare all the options as FITS keywords to write in output later. */
  gal_options_as_fits_keywords(&p->cp);


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
