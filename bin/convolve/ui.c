/*********************************************************************
Convolve - Convolve input data with a given kernel.
Convolve is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <gnuastro/fits.h>
#include <gnuastro/tile.h>
#include <gnuastro/blank.h>
#include <gnuastro/table.h>
#include <gnuastro/threads.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/statistics.h>
#include <gnuastro/linkedlist.h>

#include <timing.h>
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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" can be used to view the "
  "information, select columns, or convert tables. The inputs and outputs "
  "can be plain text (with whitespace or comma as delimiters), FITS ascii, "
  "or FITS binary tables. The output columns can either be selected by "
  "number (counting from 1), name or using regular expressions. For regular "
  "expressions, enclose the value to the `--column' (`-c') option in "
  "slashes (`\\', as in `-c\\^mag\\'). To print the selected columns on the "
  "command-line, don't specify an output file.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




/* Option groups particular to this program. */
enum program_args_groups
{
  ARGS_GROUP_MESH_GRID = GAL_OPTIONS_GROUP_AFTER_COMMON,
};




















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
      case GAL_OPTIONS_KEY_SEARCHIN:
      case GAL_OPTIONS_KEY_IGNORECASE:
      case GAL_OPTIONS_KEY_TABLEFORMAT:
        cp->coptions[i].flags=OPTION_HIDDEN;
        break;
      }
}





/* Parse a single option: */
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct convolveparams *p = state->input;

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
  /* Read the domain from a string into an integer. */
  if( !strcmp("spatial", p->domainstr) )
    p->domain=CONVOLVE_DOMAIN_SPATIAL;
  else if( !strcmp("frequency", p->domainstr) )
    p->domain=CONVOLVE_DOMAIN_FREQUENCY;
  else
    error(EXIT_FAILURE, 0, "domain value `%s' not recognized. Please use "
          "either `spatial' or `frequency'", p->domainstr);

  /* If we are in the spatial domain, make sure that the necessary
     parameters are set. */
  if( p->domain==CONVOLVE_DOMAIN_SPATIAL )
    if( p->tile==NULL || p->numchannels==NULL )
      {
        if( p->tile==NULL && p->numchannels==NULL )
          error(EXIT_FAILURE, 0, "in spatial convolution, `--numchannels' "
                "and `--tile' are mandatory");
        else
          error(EXIT_FAILURE, 0, "in spatial convolution, `--%s' is "
                "mandatory: you should use it to set the %s",
                p->tile ? "numchannels" : "tile",
                ( p->tile
                  ? "number of channels along each dimension of the input"
                  : "size of tiles to cover the input along each "
                  "dimension" ) );
      }
}





static void
ui_check_options_and_arguments(struct convolveparams *p)
{
  /* Make sure an input file name was given and if it was a FITS file, that
     a HDU is also given. */
  if(p->filename==NULL)
    error(EXIT_FAILURE, 0, "no input file is specified");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
/* Read the kernel. VERY IMPORTANT: We can't use the `fits_img_read_kernel'
   because the Convolve program also does de-convolution. */
static void
ui_read_kernel(struct convolveparams *p)
{
  float *f, *ff;

  /* Read the image into file. */
  p->kernel = gal_fits_img_read_to_type(p->kernelname, p->khdu,
                                        GAL_DATA_TYPE_FLOAT32,
                                        p->cp.minmapsize);

  /* Convert all the NaN pixels to zero if the kernel contains blank
     pixels. */
  if(gal_blank_present(p->kernel))
    {
      ff = (f=p->kernel->array) + p->kernel->size;
      do *f = isnan(*f) ? 0.0f : *f; while(++f<ff);
    }
}





static void
ui_preparations(struct convolveparams *p)
{
  size_t i, size;
  gal_data_t *sum;
  float *kernel, tmp;
  char *outsuffix = p->makekernel ? "_kernel.fits" : "_convolved.fits";


  /* Set the output name if the user hasn't set it. */
  if(p->cp.output==NULL)
    p->cp.output=gal_checkset_automatic_output(&p->cp, p->filename,
                                               outsuffix);
  gal_checkset_check_remove_file(p->cp.output, 0, p->cp.dontdelete);
  if(p->checkfreqsteps)
    {
      p->freqstepsname=gal_checkset_automatic_output(&p->cp, p->filename,
                                                     "_freqsteps.fits");
      gal_checkset_check_remove_file(p->freqstepsname, 0, p->cp.dontdelete);
    }
  if(p->checktiles)
    {
      p->tilesname=gal_checkset_automatic_output(&p->cp, p->filename,
                                                 "_tiled.fits");
      gal_checkset_check_remove_file(p->tilesname, 0, p->cp.dontdelete);
    }


  /* Read the input image as a float64 array and its WCS info. */
  p->input=gal_fits_img_read_to_type(p->filename, p->cp.hdu,
                                     GAL_DATA_TYPE_FLOAT64, p->cp.minmapsize);
  gal_wcs_read(p->filename, p->cp.hdu, 0, 0, &p->input->nwcs, &p->input->wcs);


  /* See if there are any blank values. */
  if(p->domain==CONVOLVE_DOMAIN_FREQUENCY)
    {
      if( gal_blank_present(p->input) )
        fprintf(stderr, "\n----------------------------------------\n"
                "######## %s WARNING ########\n"
                "There are blank pixels in `%s' (hdu: `%s') and you have "
                "asked for frequency domain convolution. As a result, all "
                "the pixels in the output (`%s') will be blank. Only "
                "spatial domain convolution can account for blank pixels "
                "in the input data. You can run %s again with "
                "`--domain=spatial'\n"
                "----------------------------------------\n\n",
                PROGRAM_NAME, p->filename, p->cp.hdu, p->cp.output,
                PROGRAM_NAME);
    }
  else
    p->channel=gal_tile_all_sanity_check(p->filename, p->cp.hdu, p->input,
                                         p->tile, p->numchannels);



  /* Read the file specified by --kernel. If makekernel is specified, then
     this is actually the sharper image and the input image (given as an
     argument) is the blurry image. */
  if(p->makekernel)
    {
      /* Read in the kernel array. */
      ui_read_kernel(p);

      /* Make sure the size of the kernel is the same as the input */
      if( p->input->dsize[0]!=p->kernel->dsize[0]
          || p->input->dsize[1]!=p->kernel->dsize[1] )
        error(EXIT_FAILURE, 0, "with the `--makekernel' (`-m') option, "
              "the input image and the image specified with the `--kernel' "
              "(`-k') option should have the same size. The lower resolution "
              "input image (%s) has %zux%zu pixels while the sharper image "
              "(%s) specified with the kernel option has %zux%zu pixels",
              p->filename, p->input->dsize[1], p->input->dsize[0],
              p->kernelname, p->kernel->dsize[1], p->kernel->dsize[0]);

      /* Divide both images by their sum so their lowest frequency becomes
         1 and their division (in the frequency domain) would be
         meaningful. */
      sum=gal_statistics_sum(p->input);
      p->input = gal_arithmetic(GAL_ARITHMETIC_OP_DIVIDE,
                                GAL_ARITHMETIC_FLAGS_ALL, p->input, sum);
      sum=gal_statistics_sum(p->kernel);
      p->kernel = gal_arithmetic(GAL_ARITHMETIC_OP_DIVIDE,
                                GAL_ARITHMETIC_FLAGS_ALL, p->kernel, sum);
    }

  /* Read the kernel. If there is anything particular to Convolve, then
     don't use the standard kernel reading function in fits.c. Otherwise
     just use the same one that all programs use. The standard one is
     faster because it mixes the NaN conversion and also the normalization
     into one loop. */
  else
    {
      if(p->nokernelnorm || p->nokernelflip)
        {
          /* Read in the kernel array: */
          ui_read_kernel(p);

          /* Check its size (must be odd). */
          if(p->kernel->dsize[0]%2==0 || p->kernel->dsize[1]%2==0)
            error(EXIT_FAILURE, 0, "the kernel image has to have an odd "
                  "number of pixels on both sides (there has to be on pixel "
                  "in the center). %s (hdu: %s) is %zu by %zu",
                  p->kernelname, p->khdu, p->kernel->dsize[1],
                  p->kernel->dsize[0]);

          /* Normalize the kernel: */
          if( !p->nokernelnorm )
            {
              sum=gal_statistics_sum(p->kernel);
              p->kernel = gal_arithmetic(GAL_ARITHMETIC_OP_DIVIDE,
                                         GAL_ARITHMETIC_FLAGS_ALL,
                                         p->kernel, sum);
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
      else
        p->kernel = gal_fits_img_read_kernel(p->kernelname, p->khdu,
                                             p->cp.minmapsize);
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
  printf("  - Input: %s (hdu: %s)\n", p->filename, p->cp.hdu);
  printf("  - Kernel: %s (hdu: %s)\n", p->kernelname, p->khdu);
}





void
ui_read_check_inputs_setup(int argc, char *argv[], struct convolveparams *p)
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


  /* Do a sanity check only on options. */
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
