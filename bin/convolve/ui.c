/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <gnuastro/fits.h>
#include <gnuastro/tile.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>
#include <gnuastro/array.h>
#include <gnuastro/threads.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>
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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will convolve an input "
  "image with a given spatial kernel (image) in the spatial domain (no "
  "edge effects) or frequency domain. The latter suffers from edge effects, "
  "but can be much faster.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct convolveparams *p,
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

  /* Set the mandatory common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    switch(cp->coptions[i].key)
      {
      case GAL_OPTIONS_KEY_HDU:
      case GAL_OPTIONS_KEY_TYPE:
      case GAL_OPTIONS_KEY_MINMAPSIZE:
        cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
        break;

      case GAL_OPTIONS_KEY_LOG:
      case GAL_OPTIONS_KEY_IGNORECASE:
      case GAL_OPTIONS_KEY_INTERPNUMNGB:
      case GAL_OPTIONS_KEY_INTERPONLYBLANK:
        cp->coptions[i].flags=OPTION_HIDDEN;
        break;
      }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct convolveparams *p = state->input;

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
      if(p->filename)
        argp_error(state, "only one argument (input file) should be given");
      else
        p->filename=arg;
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
ui_read_check_only_options(struct convolveparams *p)
{
  struct gal_options_common_params *cp=&p->cp;


  /* Read the domain from a string into an integer. */
  if( !strcmp("spatial", p->domainstr) )
    p->domain=CONVOLVE_DOMAIN_SPATIAL;
  else if( !strcmp("frequency", p->domainstr) )
    p->domain=CONVOLVE_DOMAIN_FREQUENCY;
  else
    error(EXIT_FAILURE, 0, "domain value '%s' not recognized. Please use "
          "either 'spatial' or 'frequency'", p->domainstr);


  /* If we are in the spatial domain, make sure that the necessary
     parameters are set. */
  if( p->domain==CONVOLVE_DOMAIN_SPATIAL )
    if( cp->tl.tilesize==NULL || cp->tl.numchannels==NULL )
      {
        if( cp->tl.tilesize==NULL && cp->tl.numchannels==NULL )
          error(EXIT_FAILURE, 0, "in spatial convolution, '--numchannels' "
                "and '--tilesize' are mandatory");
        else
          error(EXIT_FAILURE, 0, "in spatial convolution, '--%s' is "
                "mandatory: you should use it to set the %s",
                cp->tl.tilesize ? "numchannels" : "tilesize",
                ( cp->tl.tilesize
                  ? "number of channels along each dimension of the input"
                  : "size of tiles to cover the input along each "
                  "dimension" ) );
      }
}





static void
ui_check_options_and_arguments(struct convolveparams *p)
{
  int kernel_type;

  if(p->filename)
    {
      /* If input is FITS. */
      if( (p->isfits=gal_fits_name_is_fits(p->filename)) )
        {
          /* Check if a HDU is given. */
          if( p->cp.hdu==NULL )
            error(EXIT_FAILURE, 0, "no HDU specified. When the input is a "
                  "FITS file, a HDU must also be specified, you can use "
                  "the '--hdu' ('-h') option and give it the HDU number "
                  "(starting from zero), extension name, or anything "
                  "acceptable by CFITSIO");

          /* If its an image, make sure column isn't given (in case the
             user confuses an image with a table). */
          p->hdu_type=gal_fits_hdu_format(p->filename, p->cp.hdu);
          if(p->hdu_type==IMAGE_HDU && p->column)
            error(EXIT_FAILURE, 0, "%s (hdu: %s): is a FITS image "
                  "extension. The '--column' option is only applicable "
                  "to tables.", p->filename, p->cp.hdu);
        }
    }

  if(p->kernelname)
    {
      /* If input is FITS. */
      if( gal_fits_name_is_fits(p->kernelname) )
        {
          /* Check if a HDU is given. */
          if( p->khdu==NULL )
            error(EXIT_FAILURE, 0, "no HDU specified. When the kernel is a "
                  "FITS file, a HDU must also be specified, you can use "
                  "the '--khdu' ('-u') option and give it the HDU number "
                  "(starting from zero), extension name, or anything "
                  "acceptable by CFITSIO");

          /* If its an image, make sure column isn't given (in case the
             user confuses an image with a table). */
          kernel_type=gal_fits_hdu_format(p->kernelname, p->khdu);
          if(kernel_type==IMAGE_HDU && p->kernelcolumn)
            error(EXIT_FAILURE, 0, "%s (hdu: %s): is a FITS image "
                  "extension. The '--kernelcolumn' option is only "
                  "applicable to tables.", p->kernelname, p->khdu);
        }
    }
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static gal_data_t *
ui_read_column(struct convolveparams *p, int i0k1)
{
  int tformat;
  char *source;
  gal_data_t *out, *cinfo;
  gal_list_str_t *column=NULL;
  size_t ncols, nrows, counter=0;
  char *hdu        = i0k1==0 ? p->cp.hdu   : p->khdu;
  char *name       = i0k1==0 ? "input"     : "kernel";
  char *filename   = i0k1==0 ? p->filename : p->kernelname;
  char *columnname = i0k1==0 ? p->column   : p->kernelcolumn;
  gal_list_str_t *lines = gal_options_check_stdin(filename,
                                                  p->cp.stdintimeout, name);

  /* If no column is specified, Convolve will abort and an error will be
     printed when the table has more than one column. If there is only one
     column, there is no need to specify any, so Convolve will use it. */
  if(columnname==NULL)
    {
      /* Get the basic table information. */
      cinfo=gal_table_info(filename, hdu, lines, &ncols, &nrows, &tformat);
      gal_data_array_free(cinfo, ncols, 1);

      /* See how many columns it has and take the proper action. */
      switch(ncols)
        {
        case 0:
          error(EXIT_FAILURE, 0, "%s contains no usable information",
                ( filename
                  ? gal_checkset_dataset_name(filename, hdu)
                  : "Standard input" ));
        case 1:
          gal_checkset_allocate_copy("1", &columnname);
          break;
        default:
          error(EXIT_FAILURE, 0, "%s is a table containing more than one "
                "column. However, the specific column to work on isn't "
                "specified.\n\n"
                "Please use the '--column' ('-c') or '--kernelcolumn' "
                "options (depending on which dataset it is) to specify a "
                "column. You can either give it the column number "
                "(couting from 1), or a match/search in its meta-data (e.g., "
                "column names).\n\n"
                "For more information, please run the following command "
                "(press the 'SPACE' key to go down and 'q' to return to the "
                "command-line):\n\n"
                "    $ info gnuastro \"Selecting table columns\"\n",
                ( filename
                  ? gal_checkset_dataset_name(filename, hdu)
                  : "Standard input" ));
        }

    }
  gal_list_str_add(&column, columnname, 0);

  /* Read the desired column(s). */
  out=gal_table_read(filename, hdu, lines, column, p->cp.searchin,
                     p->cp.ignorecase, p->cp.minmapsize, p->cp.quietmmap,
                     NULL);
  gal_list_str_free(lines, 1);

  /* Confirm if only one column was read (it is possible to match more than
     one column). */
  if(out->next!=NULL)
    {
      if(filename)
        source=gal_checkset_dataset_name(filename, hdu);
      else
        source="standard-input";
      error(EXIT_FAILURE, 0, "%s: more than one column in input table mached "
            "the search criteria. Please limit the match by specifying the "
            "exact name (if its unique) or column number", source);
    }

  /* Make sure it is a usable datatype. */
  switch(out->type)
    {
    case GAL_TYPE_BIT:
    case GAL_TYPE_STRLL:
    case GAL_TYPE_STRING:
    case GAL_TYPE_COMPLEX32:
    case GAL_TYPE_COMPLEX64:
      error(EXIT_FAILURE, 0, " read column number %zu has a %s type, "
            "which is not currently supported by %s", counter,
            gal_type_name(out->type, 1), PROGRAM_NAME);
    }
  out=gal_data_copy_to_new_type_free(out, INPUT_USE_TYPE);

  /* If the input was from standard input, we can actually write this into
     it (for future reporting). */
  if(filename==NULL)
    gal_checkset_allocate_copy("standard-input",
                               ( i0k1==0 ? &p->filename : &p->kernelname ));

  /* Clean up and return. */
  gal_list_str_free(column, 0);
  return out;
}





/* Read the input dataset. */
static void
ui_read_input(struct convolveparams *p)
{
  /* To see if we should read it as a table. */
  p->input=NULL;

  /* If the input is a FITS image or any recognized array file format, then
     read it as an array, otherwise, as a table. */
  if( p->filename && gal_array_name_recognized(p->filename) )
    if (p->isfits && p->hdu_type==IMAGE_HDU)
      {
        p->input=gal_array_read_one_ch_to_type(p->filename, p->cp.hdu, NULL,
                                               INPUT_USE_TYPE,
                                               p->cp.minmapsize,
                                               p->cp.quietmmap);
        p->input->wcs=gal_wcs_read(p->filename, p->cp.hdu, 0, 0,
                                   &p->input->nwcs);
        p->input->ndim=gal_dimension_remove_extra(p->input->ndim,
                                                  p->input->dsize,
                                                  p->input->wcs);
      }

  /* The input isn't an image (wasn't read yet), so we'll read it as a
     column. */
  if(p->input==NULL)
    p->input=ui_read_column(p, 0);
}





/* Read the kernel. VERY IMPORTANT: We can't use the 'fits_img_read_kernel'
   because the Convolve program also does de-convolution. */
static void
ui_read_kernel(struct convolveparams *p)
{
  /* Read the image into file. */
  if( p->kernelname
      && p->input->ndim>1
      && gal_array_name_recognized(p->kernelname)  )
    {
      p->kernel = gal_array_read_one_ch_to_type(p->kernelname, p->khdu,
                                                NULL, INPUT_USE_TYPE,
                                                p->cp.minmapsize,
                                                p->cp.quietmmap);
      p->kernel->ndim=gal_dimension_remove_extra(p->kernel->ndim,
                                                 p->kernel->dsize,
                                                 p->kernel->wcs);
    }
  else
    p->kernel=ui_read_column(p, 1);

  /* Make sure that the kernel and input have the same number of
     dimensions. */
  if(p->kernel->ndim!=p->input->ndim)
    error(EXIT_FAILURE, 0, "input datasets must have the same number of "
          "dimensions");
}





static void
ui_preparations(struct convolveparams *p)
{
  int check=0;
  double sumv=0;
  size_t i, size;
  gal_data_t *sum;
  float *f, *fp, tmp, *kernel;
  struct gal_options_common_params *cp=&p->cp;
  char *outsuffix = p->makekernel ? "_kernel.fits" : "_convolved.fits";


  /* Read the input dataset. */
  ui_read_input(p);


  /* Currently Convolve only works on 1D, 2D and 3D datasets. */
  if(p->input->ndim>3)
    error(EXIT_FAILURE, 0, "%s (hdu %s) has %zu dimensions. Currently "
          "Convolve only operates on 1D (table column, spectrum), 2D "
          "(image), and 3D (data cube) datasets", p->filename, cp->hdu,
          p->input->ndim);


  /* Domain-specific checks. */
  if(p->domain==CONVOLVE_DOMAIN_FREQUENCY)
    {
      /* Check the dimensionality. */
      if(p->input->ndim!=2)
        error(EXIT_FAILURE, 0, "%s (hdu %s) has %zu dimensions. Frequency "
              "domain convolution currently only operates on 2D images",
              p->filename, cp->hdu, p->input->ndim);

      /* Blank values. */
      if( gal_blank_present(p->input, 1) )
        fprintf(stderr, "\n----------------------------------------\n"
                "######## %s WARNING ########\n"
                "There are blank pixels in '%s' (hdu: '%s') and you have "
                "asked for frequency domain convolution. As a result, all "
                "the pixels in the output ('%s') will be blank. Only "
                "spatial domain convolution can account for blank pixels "
                "in the input data. You can run %s again with "
                "'--domain=spatial'\n"
                "----------------------------------------\n\n",
                PROGRAM_NAME, p->filename, cp->hdu, cp->output,
                PROGRAM_NAME);

      /* Frequency domain is only implemented in 2D. */
      if( p->input->ndim==1 )
        error(EXIT_FAILURE, 0, "Frequency domain convolution is currently "
              "not implemented on 1D datasets. Please use '--domain=spatial' "
              "to convolve this dataset");
    }
  else
    {
      if(p->input->ndim>1)
        gal_tile_full_sanity_check(p->filename, cp->hdu, p->input, &cp->tl);
    }


  /* Read the file specified by --kernel. If makekernel is specified, then
     this is actually the sharper image and the input image (given as an
     argument) is the blurry image. */
  if(p->makekernel)
    {
      /* Currently this is not implemented in 1D. */
      if(p->kernel->ndim==1)
        error(EXIT_FAILURE, 0, "'--makekernel' is currently not available "
              "on 1D datasets");
      else
        {
          /* Read in the kernel array. */
          ui_read_kernel(p);

          /* Make sure the size of the kernel is the same as the input */
          if( p->input->dsize[0]!=p->kernel->dsize[0]
              || p->input->dsize[1]!=p->kernel->dsize[1] )
            error(EXIT_FAILURE, 0, "with the '--makekernel' ('-m') option, "
                  "the input image and the image specified with the "
                  "'--kernel' ('-k') option should have the same size. The "
                  "lower resolution input image (%s) has %zux%zu pixels "
                  "while the sharper image (%s) specified with the kernel "
                  "option has %zux%zu pixels", p->filename,
                  p->input->dsize[1], p->input->dsize[0], p->kernelname,
                  p->kernel->dsize[1], p->kernel->dsize[0]);

          /* Divide both images by their sum so their lowest frequency becomes
             1 and their division (in the frequency domain) would be
             meaningful. */
          sum=gal_statistics_sum(p->input);
          sum=gal_data_copy_to_new_type_free(sum, GAL_TYPE_FLOAT32);
          p->input = gal_arithmetic(GAL_ARITHMETIC_OP_DIVIDE, 1,
                                    GAL_ARITHMETIC_FLAGS_ALL, p->input, sum);
          sum=gal_statistics_sum(p->kernel);
          sum=gal_data_copy_to_new_type_free(sum, GAL_TYPE_FLOAT32);
          p->kernel = gal_arithmetic(GAL_ARITHMETIC_OP_DIVIDE, 1,
                                     GAL_ARITHMETIC_FLAGS_ALL, p->kernel, sum);
        }
    }

  /* Read the kernel. If there is anything particular to Convolve, then
     don't use the standard kernel reading function in fits.c. Otherwise
     just use the same one that all programs use. The standard one is
     faster because it mixes the NaN conversion and also the normalization
     into one loop. */
  else
    {
      /* Read in the kernel array: */
      ui_read_kernel(p);

      /* Check if the size along each dimension of the kernel is an odd
         number. If they are all an odd number, then the for each dimension,
         check will be incremented once. */
      for(i=0;i<p->kernel->ndim;++i)
        check += p->kernel->dsize[i]%2;
      if(check!=p->kernel->ndim)
        error(EXIT_FAILURE, 0, "%s: the kernel has to have an odd number of "
              "elements in all dimensions (there has to be one element/pixel "
              "in the center). At least one of the dimensions of doesn't "
              "have an odd number of pixels",
              gal_checkset_dataset_name(p->kernelname, p->khdu));


      /* If there are any NaN pixels, set them to zero and normalize it. A
         blank pixel in a kernel is going to make a completely blank
         output.*/
      if( !p->nokernelnorm )
        {
          sumv=0;
          fp=(f=p->kernel->array)+p->kernel->size;
          do
            {
              if(isnan(*f)) *f=0.0f;
              else          sumv+=*f;
            }
          while(++f<fp);
          p->kernel->flag |= GAL_DATA_FLAG_BLANK_CH;
          p->kernel->flag &= ~GAL_DATA_FLAG_HASBLANK;
          f=p->kernel->array; do *f++ *= 1/sumv; while(f<fp);
        }

      /* Flip the kernel: */
      if( !p->nokernelflip )
        {
          size=p->kernel->size;
          kernel=p->kernel->array;
          for(i=0;i<p->kernel->size/2;++i)
            {
              tmp                     = kernel[ i            ];
              kernel[ i             ] = kernel[ size - i - 1 ];
              kernel[ size -  i - 1 ] = tmp;
            }
        }
    }


  /* Set the output name if the user hasn't set it. */
  if(cp->output==NULL)
    cp->output=gal_checkset_automatic_output(cp, p->filename, outsuffix);
  gal_checkset_writable_remove(cp->output, 0, cp->dontdelete);
  if(p->checkfreqsteps)
    {
      p->freqstepsname=gal_checkset_automatic_output(cp, p->filename,
                                                     "_freqsteps.fits");
      gal_checkset_writable_remove(p->freqstepsname, 0, cp->dontdelete);
    }
  if(cp->tl.checktiles)
    {
      cp->tl.tilecheckname=gal_checkset_automatic_output(cp, p->filename,
                                                         "_tiled.fits");
      gal_checkset_writable_remove(cp->tl.tilecheckname, 0,
                                   cp->dontdelete);
    }
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
static void
ui_print_intro(struct convolveparams *p)
{
  printf("%s started on %s", PROGRAM_NAME, ctime(&p->rawtime));
  printf("  - Using %zu CPU threads.\n", p->cp.numthreads);
  printf("  - Input: %s\n",
         gal_checkset_dataset_name(p->filename, p->cp.hdu));
  printf("  - Kernel: %s\n",
         gal_checkset_dataset_name(p->kernelname, p->khdu));
}





void
ui_read_check_inputs_setup(int argc, char *argv[], struct convolveparams *p)
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


  /* Read the configuration files and set the common values. */
  gal_options_read_config_set(&p->cp);


  /* Do a sanity check only on options. */
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


  /* Everything is ready, print the intro if not in quiet mode. */
  if(!p->cp.quiet)
    ui_print_intro(p);
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct convolveparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->khdu);
  free(p->cp.hdu);
  free(p->cp.output);
  gal_data_free(p->input);
  gal_data_free(p->kernel);

  /* Print the final message. */
  if(!p->cp.quiet)
    gal_timing_report(t1, PROGRAM_NAME" finished in: ", 0);
}
