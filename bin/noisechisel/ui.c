/*********************************************************************
NoiseChisel - Detect signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <gnuastro/array.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>

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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" Detects and segments signal "
  "that is deeply burried in noise. It employs a noise-based detection and "
  "segmentation method enabling it to be very resilient to the rich "
  "diversity of shapes in astronomical targets.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct noisechiselparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Set the necessary common parameters structure. */
  cp->program_struct     = p;
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->numthreads         = gal_threads_number();
  cp->coptions           = gal_commonopts_options;

  /* Program specific initialization. */
  p->snthresh=NAN;

  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_LOG:
        case GAL_OPTIONS_KEY_TYPE:
        case GAL_OPTIONS_KEY_SEARCHIN:
        case GAL_OPTIONS_KEY_IGNORECASE:
        case GAL_OPTIONS_KEY_STDINTIMEOUT:
          cp->coptions[i].flags=OPTION_HIDDEN;
          break;

        case GAL_OPTIONS_KEY_TILESIZE:
        case GAL_OPTIONS_KEY_MINMAPSIZE:
        case GAL_OPTIONS_KEY_NUMCHANNELS:
        case GAL_OPTIONS_KEY_INTERPMETRIC:
        case GAL_OPTIONS_KEY_INTERPNUMNGB:
        case GAL_OPTIONS_KEY_REMAINDERFRAC:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          break;

        case GAL_OPTIONS_KEY_TABLEFORMAT:
          cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
          cp->coptions[i].doc="'txt', 'fits-ascii', 'fits-binary'.";
          break;
        }
    }
}





/* Parse a single option: */
static error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct noisechiselparams *p = state->input;

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
      if(p->inputname)
        argp_error(state, "only one argument (input file) should be given");
      else
        p->inputname=arg;
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
   check in 'ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct noisechiselparams *p)
{
  /* If the convolved option is given, then the convolved HDU is also
     mandatory. */
  if(p->convolvedname && p->chdu==NULL)
    error(EXIT_FAILURE, 0, "no value given to '--chdu'. When the "
          "'--convolved' option is called (to specify a convolved image "
          "and avoid convolution) it is mandatory to also specify a HDU "
          "for it");

  /* Make sure that the no-erode-quantile is not smaller or equal to
     qthresh. */
  if( p->noerodequant <= p->qthresh)
    error(EXIT_FAILURE, 0, "the quantile for no erosion ('--noerodequant') "
          "must be larger than the base quantile threshold ('--qthresh', "
          "or '-t'). You have provided %.4f and %.4f for the former and "
          "latter, respectively", p->noerodequant, p->qthresh);

  /* For the options that make tables, the table formation option is
     mandatory. */
  if( p->checksn && p->cp.tableformat==0 )
    error(EXIT_FAILURE, 0, "'--tableformat' is necessary with the "
          "'--checksn' option.\n"
          "Please see description for '--tableformat' after running the "
          "following command for more information (use 'SPACE' to go down "
          "the page and 'q' to return to the command-line):\n\n"
          "    $ info gnuastro \"Input Output options\"");

  /* Kernel checks. */
  if(p->kernelname && strcmp(p->kernelname, UI_NO_CONV_KERNEL_NAME))
    {
      /* Check if it exists. */
      gal_checkset_check_file(p->kernelname);

      /* If its FITS, see if a HDU has been provided. */
      if( gal_fits_name_is_fits(p->kernelname) && p->khdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU specified for kernel. When the "
              "kernel is a FITS file, a HDU must also be specified. You "
              "can use the '--khdu' option and give it the HDU number "
              "(starting from zero), extension name, or anything "
              "acceptable by CFITSIO");
    }

  /* Wide kernel checks. */
  if(p->widekernelname)
    {
      /* Check if it exists. */
      gal_checkset_check_file(p->widekernelname);

      /* If its FITS, see if a HDU has been provided. */
      if( gal_fits_name_is_fits(p->widekernelname) && p->whdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU specified for the given wide kernel "
              "('%s'). When the wide kernel is a FITS file, a HDU must also "
              "be specified. You can use the '--whdu' option and give it the "
              "HDU number (starting from zero), extension name, or any "
              "HDU identifier acceptable by CFITSIO", p->widekernelname);
    }

  /* If the S/N quantile is less than 0.1 (an arbitrary small value), this
     is probably due to forgetting that this is the purity level
     (higher-is-better), not the contamination level
     (lower-is-better). This actually happened in a few cases: where we
     wanted a false detection rate of 0.0001 (a super-high value!), and
     instead of inputing 0.9999, we mistakenly gave '--snquant' a value of
     '0.0001'. We were thus fully confused with the output (an extremely
     low value) and thought its a bug, while it wasn't! */
  if(p->snquant<0.1)
    fprintf(stderr, "\nWARNING: Value of '--snquant' ('-c') is %g. Note "
            "that this is not a contamination rate (where lower is "
            "better), it is a purity rate (where higher is better). If you "
            "intentionally asked for such a low purity level, please "
            "ignore this warning\n\n", p->snquant);
}





static void
ui_check_options_and_arguments(struct noisechiselparams *p)
{
  /* Basic input file checks. */
  if(p->inputname)
    {
      /* Check if it exists. */
      gal_checkset_check_file(p->inputname);

      /* If its FITS, see if a HDU has been provided. */
      if( gal_fits_name_is_fits(p->inputname) && p->cp.hdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU specified for input. When the input "
              "is a FITS file, a HDU must also be specified, you can use "
              "the '--hdu' ('-h') option and give it the HDU number "
              "(starting from zero), extension name, or anything "
              "acceptable by CFITSIO");
    }
  else
    error(EXIT_FAILURE, 0, "no input file is specified");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_set_output_names(struct noisechiselparams *p)
{
  char *output=p->cp.output;
  char *basename = output ? output : p->inputname;

  /* Main program output. */
  if(output)
    {
      /* Delete the file if it already exists. */
      gal_checkset_writable_remove(p->cp.output, 0, p->cp.dontdelete);

      /* When the output name is given (possibly with directory
         information), the check images will also be put in that same
         directory.. */
      p->cp.keepinputdir=1;
    }
  else
    p->cp.output=gal_checkset_automatic_output(&p->cp, p->inputname,
                                               "_detected.fits");

  /* Tile check. */
  if(p->cp.tl.checktiles)
    p->cp.tl.tilecheckname=gal_checkset_automatic_output(&p->cp, basename,
                                                         "_tiles.fits");

  /* Quantile threshold. */
  if(p->checkqthresh)
    p->qthreshname=gal_checkset_automatic_output(&p->cp, basename,
                                                 "_qthresh.fits");

  /* Initial detection Sky values. */
  if(p->checkdetsky)
    p->detskyname=gal_checkset_automatic_output(&p->cp, basename,
                                                "_detsky.fits");

  /* Pseudo-detection S/N values. */
  if(p->checksn)
    {
      p->detsn_s_name=gal_checkset_automatic_output(&p->cp, basename,
                 ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                   ? "_detsn_sky.txt" : "_detsn.fits") );
      p->detsn_d_name=gal_checkset_automatic_output(&p->cp, basename,
                 ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                   ? "_detsn_det.txt" : "_detsn.fits") );
      p->detsn_D_name=gal_checkset_automatic_output(&p->cp, basename,
                 ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                   ? "_detsn_grown.txt" : "_detsn.fits") );
    }

  /* Detection steps. */
  if(p->checkdetection)
    p->detectionname=gal_checkset_automatic_output(&p->cp, basename,
                                               "_detcheck.fits");

  /* Sky checks. */
  if(p->checksky)
    p->skyname=gal_checkset_automatic_output(&p->cp, basename, "_sky.fits");
}





static void
ui_prepare_kernel(struct noisechiselparams *p)
{
  float *f, *ff, *k;

/* Impor the default kernel. */
#include "kernel-2d.h"

  /* If a kernel file is given, then use it. Otherwise, use the default
     kernel. */
  if(p->kernelname)
    {
      if( strcmp(p->kernelname, UI_NO_CONV_KERNEL_NAME) )
        p->kernel=gal_fits_img_read_kernel(p->kernelname, p->khdu,
                                           p->cp.minmapsize, p->cp.quietmmap);
      else
        p->kernel=NULL;
    }
  else
    {
      /* Allocate space for the kernel (we don't want to use the statically
         allocated array. */
      p->kernel=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 2,
                               kernel_2d_dsize, NULL, 0, p->cp.minmapsize,
                               p->cp.quietmmap, NULL, NULL, NULL);

      /* Now copy the staticly allocated array into it. */
      k=p->kernel->array;
      ff=(f=kernel_2d)+gal_dimension_total_size(2, p->kernel->dsize);
      do *k++=*f; while(++f<ff);
    }


  /* If a wide kernel is given, then read it into memory. Otherwise, just
     ignore it. */
  if(p->widekernelname)
    p->widekernel=gal_fits_img_read_kernel(p->widekernelname, p->whdu,
                                           p->cp.minmapsize, p->cp.quietmmap);
}





static void
ui_prepare_tiles(struct noisechiselparams *p)
{
  gal_data_t *check;
  struct gal_tile_two_layer_params *tl=&p->cp.tl, *ltl=&p->ltl;


  /* Check the tile parameters for the small tile sizes and make the tile
     structure. We will also need the dimensions of the tile with the
     maximum required memory. */
  p->maxtsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->input->ndim, 0,
                                   __func__, "p->maxtsize");
  gal_tile_full_sanity_check(p->inputname, p->cp.hdu, p->input, tl);
  gal_tile_full_two_layers(p->input, tl);
  gal_tile_full_permutation(tl);
  for(check=tl->tiles; check!=NULL; check=check->next)
    if( check->size > p->maxtcontig )/* p->maxtcontig initialized to 0. */
      {
        p->maxtcontig=check->size;
        memcpy(p->maxtsize, check->dsize, tl->ndim*sizeof *p->maxtsize);
      }


  /* Make the large tessellation, except for the size, the rest of the
     parameters are the same as the small tile sizes. */
  ltl->numchannels    = tl->numchannels;
  ltl->remainderfrac  = tl->remainderfrac;
  ltl->workoverch     = tl->workoverch;
  ltl->checktiles     = tl->checktiles;
  ltl->oneelempertile = tl->oneelempertile;
  p->maxltsize=gal_pointer_allocate(GAL_TYPE_SIZE_T, p->input->ndim, 0,
                                    __func__, "p->maxltsize");
  gal_tile_full_sanity_check(p->inputname, p->cp.hdu, p->input, ltl);
  gal_tile_full_two_layers(p->input, ltl);
  gal_tile_full_permutation(ltl);
  for(check=ltl->tiles; check!=NULL; check=check->next)
    if( check->size > p->maxltcontig )/* p->maxltcontig initialized to 0. */
      {
        p->maxltcontig=check->size;
        memcpy(p->maxltsize, check->dsize, ltl->ndim*sizeof *p->maxltsize);
      }


  /* If the input has blank elements, then set teh appropriate flag for
     each tile.*/
  if( p->input->flag & GAL_DATA_FLAG_HASBLANK )
    {
      gal_tile_block_blank_flag(tl->tiles,  p->cp.numthreads);
      gal_tile_block_blank_flag(ltl->tiles, p->cp.numthreads);
    }


  /* Make the tile check image if requested. */
  if(tl->checktiles)
    {
      /* Large tiles. */
      check=gal_tile_block_check_tiles(ltl->tiles);
      gal_fits_img_write(check, tl->tilecheckname, NULL, PROGRAM_NAME);
      gal_data_free(check);

      /* Small tiles. */
      check=gal_tile_block_check_tiles(tl->tiles);
      gal_fits_img_write(check, tl->tilecheckname, NULL, PROGRAM_NAME);
      gal_data_free(check);

      /* If 'continueaftercheck' hasn't been called, abort NoiseChisel. */
      if(!p->continueaftercheck)
        ui_abort_after_check(p, tl->tilecheckname, NULL,
                             "showing all tiles over the image");

      /* Free the name. */
      free(tl->tilecheckname);
    }
}





static void
ui_ngb_check(size_t value, char *optionname, size_t ndim)
{
  switch(ndim)
    {
    case 2:
      if(value!=4 && value!=8)
        error(EXIT_FAILURE, 0, "%zu is not an acceptable value for "
              "'--%s'. Acceptable values for 2D inputs are 4 or 8",
              value, optionname);
      break;
    case 3:
      error(EXIT_FAILURE, 0, "3D input data is not yet supported");
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
            "problem. Dimension value %zu is not recognized.", __func__,
            PACKAGE_BUGREPORT, ndim);
    }
}





/* Read the input image and do the basic checks */
static void
ui_preparations_read_input(struct noisechiselparams *p)
{
  float *f;
  size_t ndim;

  /* Read the input as a single precision floating point dataset, also load
     the WCS and finally remove any possibly existing extra dimensions
     (with a length of 1). */
  p->input = gal_array_read_one_ch_to_type(p->inputname, p->cp.hdu,
                                           NULL, GAL_TYPE_FLOAT32,
                                           p->cp.minmapsize,
                                           p->cp.quietmmap);
  p->input->wcs = gal_wcs_read(p->inputname, p->cp.hdu, 0, 0,
                               &p->input->nwcs);
  p->input->ndim=gal_dimension_remove_extra(p->input->ndim,
                                            p->input->dsize,
                                            p->input->wcs);

  /* When the input doesn't have a name, use 'INPUT'. */
  if(p->input->name==NULL)
    gal_checkset_allocate_copy("INPUT", &p->input->name);

  /* NoiseChisel currently only works on 2D datasets (images). */
  if(p->input->ndim!=2)
    error(EXIT_FAILURE, 0, "%s (hdu: %s) has %zu dimensions but NoiseChisel "
          "can only operate on 2D datasets (images)", p->inputname, p->cp.hdu,
          p->input->ndim);

  /* Check the values of dimension-related options. */
  ndim=p->input->ndim;
  ui_ngb_check(p->holengb, "holengb", ndim);
  ui_ngb_check(p->erodengb, "erodengb", ndim);
  ui_ngb_check(p->openingngb, "openingngb", ndim);
  ui_ngb_check(p->dopeningngb, "dopeningngb", ndim);
  ui_ngb_check(p->pseudoconcomp, "pseudoconcomp", ndim);

  /* A small check to see if the edges of the dataset aren't zero valued:
     they should be masked. */
  f=p->input->array;
  if( (f[0]==0.0 && f[1]==0.0)
      || (f[ p->input->size-1 ]==0.0 && f[ p->input->size-2 ]==0.0) )
    error(0, 0, "%s (hdu %s): [*** WARNING ***]: The first and/or last few "
          "pixels have a value of 0.0. As described below, the result of "
          "this run may thus not be reasonable/optimal.\n\n"
          "Some data reduction pipelines put 0.0 where there isn't data "
          "(most commonly on the edges). However, NoiseChisel's "
          "noise-based detection paradigm starts from the lower values of "
          "the dataset (not high S/N peaks): its initial threshold is "
          "mostly below the Sky value (0.0 in processed images). Therefore "
          "0.0 is meaningful for NoiseChisel and must not be used for a "
          "blank value.\n\n"
          "To ignore certain pixels, they must have a blank/NaN value. "
          "To mask (set to blank/NaN) the 0.0 valued elements, you can use "
          "Gnuastro's Arithmetic program with a command like this:\n\n"
          "    $ astarithmetic %s %s 0.0 eq nan where -g%s\n\n"
          "If the few 0.0 valued pixels on the edges are meaningful for "
          "your analysis, please ignore this warming message.\n"
          "--------------------------",
          p->inputname, p->cp.hdu, p->inputname, p->inputname, p->cp.hdu);
}





static void
ui_preparations(struct noisechiselparams *p)
{
  /* Prepare the names of the outputs. */
  ui_set_output_names(p);

  /* Read the input datasets and do the basic checks.*/
  ui_preparations_read_input(p);

  /* If a convolved image was given, read it in. Otherwise, read the given
     kernel. */
  if(p->convolvedname)
    {
      /* Read the input convolved image. */
      p->conv = gal_array_read_one_ch_to_type(p->convolvedname, p->chdu,
                                              NULL, GAL_TYPE_FLOAT32,
                                              p->cp.minmapsize,
                                              p->cp.quietmmap);

      /* Make sure the convolved image is the same size as the input. */
      if( gal_dimension_is_different(p->input, p->conv) )
        error(EXIT_FAILURE, 0, "%s (hdu %s), given to '--convolved' and "
              "'--convolvehdu', is not the same size as NoiseChisel's "
              "input: %s (hdu: %s)", p->convolvedname, p->chdu,
              p->inputname, p->cp.hdu);
    }
  else
    ui_prepare_kernel(p);

  /* Check for blank values to help later processing.  */
  gal_blank_present(p->input, 1);

  /* Prepare the tessellation. */
  ui_prepare_tiles(p);

  /* Allocate space for the over-all necessary arrays. */
  p->binary=gal_data_alloc(NULL, GAL_TYPE_UINT8, p->input->ndim,
                           p->input->dsize, p->input->wcs, 0,
                           p->cp.minmapsize, p->cp.quietmmap, NULL,
                           "binary", NULL);
  p->olabel=gal_data_alloc(NULL, GAL_TYPE_INT32, p->input->ndim,
                           p->input->dsize, p->input->wcs, 0,
                           p->cp.minmapsize, p->cp.quietmmap, NULL,
                           "labels", NULL);
  p->binary->flag = p->olabel->flag = p->input->flag;
}



















/**************************************************************/
/************     High level reading function     *************/
/**************************************************************/
void
ui_read_check_inputs_setup(int argc, char *argv[],
                           struct noisechiselparams *p)
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


  /* Read the options into the program's structure, and check them and
     their relations prior to printing. */
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


  /* Let the user know that processing has started. */
  if(!p->cp.quiet)
    {
      printf(PROGRAM_NAME" "PACKAGE_VERSION" started on %s",
             ctime(&p->rawtime));
      printf("  - Using %zu CPU thread%s\n", p->cp.numthreads,
             p->cp.numthreads==1 ? "." : "s.");
      printf("  - Input: %s (hdu: %s)\n", p->inputname, p->cp.hdu);
      if(p->convolvedname)
        printf("  - Convolved input: %s (hdu: %s)\n",
               p->convolvedname, p->chdu);
      else
        {
          if(p->kernelname)
            {
              if( strcmp(p->kernelname, UI_NO_CONV_KERNEL_NAME) )
                printf("  - %s: %s (hdu: %s)\n",
                       p->widekernelname ? "Sharp Kernel" : "Kernel",
                       p->kernelname, p->khdu);
              else
                printf("  - No convolution requested.\n");
            }
          else
            printf("  - %s: FWHM=2 pixel Gaussian.\n",
                   p->widekernelname ? "Sharp Kernel" : "Kernel");
        }
      if(p->widekernelname)
        printf("  - Wide Kernel: %s (hdu: %s)\n", p->widekernelname,
               p->whdu);
    }
}




















/**************************************************************/
/************     Pre-finish/abort operations     *************/
/**************************************************************/
void
ui_abort_after_check(struct noisechiselparams *p, char *filename,
                     char *file2name, char *description)
{
  char *name;

  if(file2name)
    {
      if( asprintf(&name, "'%s' and '%s'", filename, file2name)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
    }
  else
    {
      if( asprintf(&name, "'%s'", filename)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
    }

  /* Let the user know that NoiseChisel is aborting. */
  fprintf(stderr,
          "------------------------------------------------\n"
          "%s aborted for a check\n"
          "------------------------------------------------\n"
          "%s (%s) has been created.\n\n"
          "If you want %s to continue its processing AND save any "
          "requested check outputs, please run it again with "
          "'--continueaftercheck'.\n"
          "------------------------------------------------\n",
          PROGRAM_NAME, name, description, PROGRAM_NAME);

  /* Clean up. */
  free(name);
  ui_free_report(p, NULL);

  /* Abort. */
  exit(EXIT_SUCCESS);
}





void
ui_free_report(struct noisechiselparams *p, struct timeval *t1)
{
  /* Free the simply allocated spaces. */
  free(p->cp.hdu);
  free(p->maxtsize);
  free(p->maxltsize);
  free(p->cp.output);
  if(p->skyname)          free(p->skyname);
  if(p->detskyname)       free(p->detskyname);
  if(p->qthreshname)      free(p->qthreshname);
  if(p->detsn_s_name)     free(p->detsn_s_name);
  if(p->detsn_d_name)     free(p->detsn_d_name);
  if(p->detectionname)    free(p->detectionname);

  /* Free the allocated datasets. */
  gal_data_free(p->sky);
  gal_data_free(p->std);
  gal_data_free(p->wconv);
  gal_data_free(p->input);
  gal_data_free(p->kernel);
  gal_data_free(p->binary);
  gal_data_free(p->olabel);
  gal_data_free(p->widekernel);
  if(p->conv!=p->input) gal_data_free(p->conv);

  /* Clean up the tile structure. */
  p->ltl.numchannels=NULL;
  gal_tile_full_free_contents(&p->ltl);
  gal_tile_full_free_contents(&p->cp.tl);

  /* Print the final message. */
  if(!p->cp.quiet && t1)
    gal_timing_report(t1, PROGRAM_NAME" finished in: ", 0);
}
