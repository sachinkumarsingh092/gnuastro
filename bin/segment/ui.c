/*********************************************************************
Segment - Segment initial labels based on signal structure.
Segment is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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
#include <gnuastro/binary.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>
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
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will segment an initially "
  "labeled region based on structure with the signal. It will first find "
  "true clumps (local maxima), estimate which ones have strong connections, "
  "and then grow them to cover the full area of each detection.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;




















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct segmentparams *p,
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

  p->medstd              = NAN;
  p->minstd              = NAN;
  p->maxstd              = NAN;
  p->snquant             = NAN;
  p->clumpsnthresh       = NAN;

  /* Modify common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    {
      /* Select individually. */
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_HDU:
          cp->coptions[i].doc="HDU containing values (science image).";
          break;

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
error_t
parse_opt(int key, char *arg, struct argp_state *state)
{
  struct segmentparams *p = state->input;

  /* Pass 'ygal_options_common_params' into the child parser.  */
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
ui_read_check_only_options(struct segmentparams *p)
{
  /* If the full area is to be used as a single detection, we can't find
     the S/N value from the un-detected regions, so the user must have
     given the 'clumpsnthresh' option. */
  if( p->detectionname
      && !strcmp(p->detectionname, DETECTION_ALL)
      && isnan(p->clumpsnthresh) )
    error(EXIT_FAILURE, 0, "'--clumpsnthresh' ('-%c') not given.\n\n"
          "When '--detection=all' (the whole input dataset is assumed to "
          "be a detection), Segment can't use the undetected pixels to find "
          "the signal-to-noise ratio of true clumps. Therefore it is "
          "mandatory to provide a signal-to-noise ratio manually",
          UI_KEY_CLUMPSNTHRESH);

  /* If the convolved HDU is given. */
  if(p->convolvedname && p->chdu==NULL)
    error(EXIT_FAILURE, 0, "no value given to '--convolvedhdu'. When the "
          "'--convolved' option is called (to specify a convolved dataset "
          "and avoid convolution) it is mandatory to also specify a HDU "
          "for it");

  /* For the options that make tables, the table format option is
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
ui_check_options_and_arguments(struct segmentparams *p)
{
  /* Make sure an input file name was given and if it was a FITS file, that
     a HDU is also given. */
  if(p->inputname)
    {
      /* Check if it exists. */
      gal_checkset_check_file(p->inputname);

      /* If it is FITS, a HDU is also mandatory. */
      if( gal_fits_name_is_fits(p->inputname) && p->cp.hdu==NULL )
        error(EXIT_FAILURE, 0, "no HDU specified. When the input is a FITS "
              "file, a HDU must also be specified, you can use the '--hdu' "
              "('-h') option and give it the HDU number (starting from "
              "zero), extension name, or anything acceptable by CFITSIO");

    }
  else
    error(EXIT_FAILURE, 0, "no input file is specified");
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static void
ui_set_used_names(struct segmentparams *p)
{
  p->useddetectionname = p->detectionname ? p->detectionname : p->inputname;

  p->usedstdname = ( p->stdname
                     ? p->stdname
                     : ( ( p->detectionname
                           && strcmp(p->detectionname, DETECTION_ALL) )
                         ? p->detectionname
                         : p->inputname ) );
}





static void
ui_set_output_names(struct segmentparams *p)
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
         directory. */
      p->cp.keepinputdir=1;
    }
  else
    p->cp.output=gal_checkset_automatic_output(&p->cp, p->inputname,
                                               "_segmented.fits");

  /* Tile check. */
  if(p->cp.tl.checktiles)
    p->cp.tl.tilecheckname=gal_checkset_automatic_output(&p->cp, basename,
                                                         "_tiles.fits");

  /* Clump S/N values. */
  if(p->checksn)
    {
      p->clumpsn_s_name=gal_checkset_automatic_output(&p->cp, basename,
                 ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                   ? "_clumpsn_sky.txt" : "_clumpsn.fits") );
      p->clumpsn_d_name=gal_checkset_automatic_output(&p->cp, basename,
                 ( p->cp.tableformat==GAL_TABLE_FORMAT_TXT
                   ? "_clumpsn_det.txt" : "_clumpsn.fits") );
    }

  /* Segmentation steps. */
  if(p->checksegmentation)
    p->segmentationname=gal_checkset_automatic_output(&p->cp, basename,
                                                      "_segcheck.fits");
}





static void
ui_prepare_inputs(struct segmentparams *p)
{
  int32_t *i, *ii;
  gal_data_t *maxd, *ccin, *blankflag, *ccout=NULL;

  /* Read the input as a single precision floating point dataset. */
  p->input = gal_array_read_one_ch_to_type(p->inputname, p->cp.hdu,
                                           NULL, GAL_TYPE_FLOAT32,
                                           p->cp.minmapsize, p->cp.quietmmap);
  p->input->wcs = gal_wcs_read(p->inputname, p->cp.hdu, 0, 0,
                               &p->input->nwcs);
  p->input->ndim=gal_dimension_remove_extra(p->input->ndim,
                                            p->input->dsize,
                                            p->input->wcs);

  /* Set the name. */
  if(p->input->name) free(p->input->name);
  gal_checkset_allocate_copy("INPUT", &p->input->name);


  /* Check for blank values to help later processing.  */
  gal_blank_present(p->input, 1);


  /* Segment currently only works on 2D datasets (images). */
  if(p->input->ndim!=2)
    error(EXIT_FAILURE, 0, "%s (hdu: %s) has %zu dimensions but Segment "
          "can only operate on 2D datasets (images)", p->inputname, p->cp.hdu,
          p->input->ndim);


  /* If a convolved image is given, read it. */
  if(p->convolvedname)
    {
      /* Read the input convolved image. */
      p->conv = gal_array_read_one_ch_to_type(p->convolvedname, p->chdu,
                                              NULL, GAL_TYPE_FLOAT32,
                                              p->cp.minmapsize,
                                              p->cp.quietmmap);
      p->conv->ndim=gal_dimension_remove_extra(p->conv->ndim,
                                               p->conv->dsize,
                                               p->conv->wcs);
      p->conv->wcs=gal_wcs_copy(p->input->wcs);

      /* Make sure it is the same size as the input. */
      if( gal_dimension_is_different(p->input, p->conv) )
        error(EXIT_FAILURE, 0, "%s (hdu %s), given to '--convolved' and "
              "'--chdu', is not the same size as the input (%s, hdu: %s)",
              p->convolvedname, p->chdu, p->inputname, p->cp.hdu);
    }


  /* Read the detected label image and check its size. When the user gives
     '--detection=all', then the whole input is assumed to be a single
     detection. */
  if( strcmp(p->useddetectionname, DETECTION_ALL) )
    {
      /* Read the dataset into memory. */
      p->olabel = gal_array_read_one_ch(p->useddetectionname, p->dhdu,
                                        NULL, p->cp.minmapsize,
                                        p->cp.quietmmap);
      p->olabel->ndim=gal_dimension_remove_extra(p->olabel->ndim,
                                                 p->olabel->dsize, NULL);
      if( gal_dimension_is_different(p->input, p->olabel) )
        error(EXIT_FAILURE, 0, "'%s' (hdu: %s) and '%s' (hdu: %s) have a"
              "different dimension/size", p->useddetectionname, p->dhdu,
              p->inputname, p->cp.hdu);

      /* Make sure the detected labels are not floating point. */
      if(p->olabel->type==GAL_TYPE_FLOAT32
         || p->olabel->type==GAL_TYPE_FLOAT64)
        error(EXIT_FAILURE, 0, "%s (hdu: %s) has a '%s' type. The detection "
              "(labeled) map must have an integer type (labels/classes can "
              "only be integers). If the pixel values are integers, but only "
              "the numerical type of the image is floating-point, you can "
              "use the command below to convert it to a 32-bit (signed) "
              "integer type:\n\n"
              "    $ astarithmetic %s int32 -h%s\n\n", p->useddetectionname,
              p->dhdu, gal_type_name(p->olabel->type, 1),
              p->useddetectionname, p->dhdu);

      /* If the input has blank values, set them to blank values in the
         labeled image too. It doesn't matter if the labeled image has
         blank pixels that aren't blank on the input image. */
      if(gal_blank_present(p->input, 1))
        {
          blankflag=gal_blank_flag(p->input);
          gal_blank_flag_apply(p->olabel, blankflag);
          gal_data_free(blankflag);
        }

      /* Get the maximum value of the input (total number of labels if they
         are separate). If the maximum is 1 (the image is a binary image),
         then apply the connected components algorithm to separate the
         connected regions. The user is allowed to supply a simple binary
         image.*/
      maxd=gal_statistics_maximum(p->olabel);
      maxd=gal_data_copy_to_new_type_free(maxd, GAL_TYPE_INT64);
      p->numdetections = *((uint64_t *)(maxd->array));
      if( p->numdetections == 1 )
        {
          ccin=gal_data_copy_to_new_type_free(p->olabel, GAL_TYPE_UINT8);
          p->numdetections=gal_binary_connected_components(ccin, &ccout,
                                                           ccin->ndim);
          gal_data_free(ccin);
          p->olabel=ccout;
        }
      else
        p->olabel = gal_data_copy_to_new_type_free(p->olabel, GAL_TYPE_INT32);


      /* Write the WCS into the objects dataset too. */
      p->olabel->wcs=gal_wcs_copy(p->input->wcs);
    }
  else
    {
      /* Set the total number of detections to 1. */
      p->numdetections=1;

      /* Allocate the array. */
      p->olabel=gal_data_alloc(NULL, GAL_TYPE_INT32, p->input->ndim,
                               p->input->dsize, p->input->wcs, 0,
                               p->cp.minmapsize, p->cp.quietmmap,
                               NULL, NULL, NULL);

      /* Initialize it to 1. */
      ii=(i=p->olabel->array)+p->olabel->size; do *i++=1; while(i<ii);
    }
}





/* Prepare the input kernel. */
static void
ui_prepare_kernel(struct segmentparams *p)
{
  float *f, *ff, *k;

/* Import the default kernel. */
#include "kernel-2d.h"

  /* If a kernel file is given, then use it. Otherwise, use the default
     kernel. */
  if(p->kernelname)
    {
      if( strcmp(p->kernelname, UI_NO_CONV_KERNEL_NAME) )
        {
          p->kernel=gal_fits_img_read_kernel(p->kernelname, p->khdu,
                                             p->cp.minmapsize,
                                             p->cp.quietmmap);
          p->kernel->ndim=gal_dimension_remove_extra(p->kernel->ndim,
                                                     p->kernel->dsize,
                                                     NULL);
        }
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
}





/* Set up the tessellation. */
static void
ui_prepare_tiles(struct segmentparams *p)
{
  gal_data_t *check;
  struct gal_tile_two_layer_params *tl=&p->cp.tl, *ltl=&p->ltl;


  /* Check the tile parameters for the small tile sizes and make the tile
     structure.  */
  gal_tile_full_sanity_check(p->inputname, p->cp.hdu, p->input, tl);
  gal_tile_full_two_layers(p->input, tl);
  gal_tile_full_permutation(tl);


  /* Make the large tessellation, except for the size, the rest of the
     parameters are the same as the small tile sizes. */
  ltl->numchannels    = tl->numchannels;
  ltl->remainderfrac  = tl->remainderfrac;
  ltl->workoverch     = tl->workoverch;
  ltl->checktiles     = tl->checktiles;
  ltl->oneelempertile = tl->oneelempertile;
  gal_tile_full_sanity_check(p->inputname, p->cp.hdu, p->input, ltl);
  gal_tile_full_two_layers(p->input, ltl);
  gal_tile_full_permutation(ltl);


  /* If the input has blank elements, then set the appropriate flag for
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
      tl->tilecheckname=NULL;
    }
}





static void
ui_check_size(gal_data_t *base, gal_data_t *comp, size_t numtiles,
              char *bname, char *bhdu, char *cname, char *chdu)
{
  if( gal_dimension_is_different(base, comp) && numtiles!=comp->size )
    error(EXIT_FAILURE, 0, "%s (hdu: %s): doesn't have the right size "
          "(%zu elements or pixels).\n\n"
          "It must either be the same size as '%s' (hdu: '%s'), or "
          "it must have the same number of elements as the total "
          "number of tiles in the tessellation (%zu). In the latter "
          "case, each pixel is assumed to be a fixed value for a "
          "complete tile.\n\n"
          "Run with '-P' to see the (tessellation) options/settings "
          "and their values). For more information on tessellation in "
          "Gnuastro, please run the following command (use the arrow "
          "keys for up and down and press 'q' to return to the "
          "command-line):\n\n"
          "    $ info gnuastro tessellation",
          cname, chdu, comp->size, bname, bhdu, numtiles);
}





/* Subtract 'sky' from the input dataset depending on its size (it may be
   the whole array or a tile-values array).. */
static void
ui_subtract_sky(gal_data_t *in, gal_data_t *sky,
                struct gal_tile_two_layer_params *tl)
{
  size_t tid;
  gal_data_t *tile;
  float *s, *f, *ff, *skyarr=sky->array;

  /* It is the same size as the input or a single value. */
  if( gal_dimension_is_different(in, sky)==0 || sky->size==1)
    {
      s=sky->array;
      ff=(f=in->array)+in->size;
      if(sky->size==1) { if(*s!=0.0) do *f-=*s;   while(++f<ff); }
      else                           do *f-=*s++; while(++f<ff);
    }

  /* It is the same size as the number of tiles. */
  else if( tl->tottiles==sky->size )
    {
      /* Go over all the tiles. */
      for(tid=0; tid<tl->tottiles; ++tid)
        {
          /* For easy reading. */
          tile=&tl->tiles[tid];

          /* Subtract the Sky value from the input image. */
          GAL_TILE_PARSE_OPERATE(tile, NULL, 0, 0, {*i-=skyarr[tid];});
        }
    }

  /* The size must have been checked before, so if control reaches here, we
     have a bug! */
  else
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
          "the problem. For some reason, the size doesn't match", __func__,
          PACKAGE_BUGREPORT);
}





/* The Sky and Sky standard deviation images can be a 'oneelempertile'
   image (only one element/pixel for a tile). So we need to do some extra
   checks on them (after reading the tessellation). */
static float
ui_read_std_and_sky(struct segmentparams *p)
{
  size_t one=1;
  char *tailptr;
  float tmpval, skyval=NAN;
  struct gal_tile_two_layer_params *tl=&p->cp.tl;
  gal_data_t *sky, *keys=gal_data_array_calloc(3);

  /* See if the name used for the standard deviation is a filename or a
     value. When the string is only a number (and nothing else), 'tailptr'
     will point to the end of the string ('\0'). When the string doesn't
     start with a number, it will point to the start of the
     string. However, file names might also be things like '1_std.fits'. In
     such cases, 'strtod' will return '1.0' and 'tailptr' will be
     '_std.fits'. Thus the most robust test is to see if 'tailptr' is the
     NULL string character. */
  tmpval=strtod(p->usedstdname, &tailptr);
  if(*tailptr=='\0')
    {
      /* Allocate the dataset to keep the value and write it in. */
      p->std=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL, 0,
                            -1, 1, NULL, NULL, NULL);
      *(float *)(p->std->array) = tmpval;
    }
  else
    {
      /* Make sure a HDU is also given. */
      if(p->stdhdu==NULL)
        error(EXIT_FAILURE, 0, "no value given to '--stdhdu'.\n\n"
              "When the Sky standard deviation is a dataset, it is mandatory "
              "specify which HDU/extension it is present in. The file can "
              "be specified explicitly with '--std'. If not, segment will "
              "use the file given to '--detection'. If that is also not "
              "called, it will look into the main input file (with no "
              "option)");

      /* Read the STD image. */
      p->std=gal_array_read_one_ch_to_type(p->usedstdname, p->stdhdu,
                                           NULL, GAL_TYPE_FLOAT32,
                                           p->cp.minmapsize, p->cp.quietmmap);
      p->std->ndim=gal_dimension_remove_extra(p->std->ndim,
                                              p->std->dsize, NULL);

      /* Make sure it has the correct size. */
      ui_check_size(p->input, p->std, tl->tottiles, p->inputname, p->cp.hdu,
                    p->usedstdname, p->stdhdu);
    }

  /* When the Standard deviation dataset (not single value) is made by
     NoiseChisel, it puts three basic statistics of the pre-interpolation
     distribution of standard deviations in 'MEDSTD', 'MINSTD' and
     'MAXSTD'. The 'MEDSTD' in particular is most important because it
     can't be inferred after the interpolations and it can be useful in
     MakeCatalog later to give a more accurate estimate of the noise
     level. So if they are present, we will read them here and write them
     to the STD output (which is created when '--rawoutput' is not
     given). */
  if(!p->rawoutput && p->std->size>1)
    {
      keys[0].next=&keys[1];
      keys[1].next=&keys[2];
      keys[2].next=NULL;
      keys[0].array=&p->medstd;     keys[0].name="MEDSTD";
      keys[1].array=&p->minstd;     keys[1].name="MINSTD";
      keys[2].array=&p->maxstd;     keys[2].name="MAXSTD";
      keys[0].type=keys[1].type=keys[2].type=GAL_TYPE_FLOAT32;
      gal_fits_key_read(p->usedstdname, p->stdhdu, keys, 0, 0);
      if(keys[0].status) p->medstd=NAN;
      if(keys[1].status) p->minstd=NAN;
      if(keys[2].status) p->maxstd=NAN;
      keys[0].name=keys[1].name=keys[2].name=NULL;
      keys[0].array=keys[1].array=keys[2].array=NULL;
      gal_data_array_free(keys, 3, 1);
    }

  /* Similar to '--std' above. */
  if(p->skyname)
    {
      tmpval=strtod(p->skyname, &tailptr);
      if(*tailptr=='\0')
        {
          sky=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &one, NULL, 0, -1,
                             1, NULL, NULL, NULL);
          *(float *)(sky->array) = skyval = tmpval;
        }
      else
        {
          /* Make sure a HDU is also given. */
          if(p->skyhdu==NULL)
            error(EXIT_FAILURE, 0, "no value given to '--skyhdu'.\n\n"
                  "When the Sky is a dataset, it is mandatory specify "
                  "which HDU/extension it is present in. The file can be "
                  "specified explicitly with '--sky'. If it is a single "
                  "value, you can just pass the value to '--sky' and no "
                  "HDU will be necessary");

          /* Read the Sky dataset. */
          sky=gal_array_read_one_ch_to_type(p->skyname, p->skyhdu,
                                            NULL, GAL_TYPE_FLOAT32,
                                            p->cp.minmapsize, p->cp.quietmmap);
          sky->ndim=gal_dimension_remove_extra(sky->ndim, sky->dsize,
                                               NULL);

          /* Check its size. */
          ui_check_size(p->input, sky, tl->tottiles, p->inputname, p->cp.hdu,
                        p->skyname, p->skyhdu);
        }

      /* Subtract the sky from the input. */
      ui_subtract_sky(p->input, sky, tl);

      /* If a convolved image is given, subtract the Sky from that too. */
      if(p->conv)
        ui_subtract_sky(p->conv, sky, tl);

      /* Clean up. */
      gal_data_free(sky);
    }

  /* Return the sky value (possibly necessary in verbose mode). */
  return skyval;
}





static float
ui_preparations(struct segmentparams *p)
{
  /* Set the input names. */
  ui_set_used_names(p);

  /* Prepare the names of the outputs. */
  ui_set_output_names(p);

  /* Read the input datasets. */
  ui_prepare_inputs(p);

  /* If a convolved image was given, read it in. Otherwise, read the given
     kernel. */
  if(p->conv==NULL)
    ui_prepare_kernel(p);

  /* Prepare the tessellation. */
  ui_prepare_tiles(p);

  /* Prepare the (optional Sky, and) Sky Standard deviation image. */
  return ui_read_std_and_sky(p);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
ui_read_check_inputs_setup(int argc, char *argv[], struct segmentparams *p)
{
  float sky;
  char *stdunit;
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
  sky=ui_preparations(p);


  /* Let the user know that processing has started. */
  if(!p->cp.quiet)
    {
      /* Basic inputs. */
      printf(PROGRAM_NAME" "PACKAGE_VERSION" started on %s",
             ctime(&p->rawtime));
      printf("  - Using %zu CPU thread%s\n", p->cp.numthreads,
             p->cp.numthreads==1 ? "." : "s.");
      printf("  - Input: %s (hdu: %s)\n", p->inputname, p->cp.hdu);

      /* Sky value information. */
      if(p->skyname)
        {
          if( isnan(sky) )
            printf("  - Sky: %s (hdu: %s)\n", p->skyname, p->skyhdu);
          else
            printf("  - Sky: %g\n", sky);
        }

      /* Sky Standard deviation information. */
      stdunit = p->variance ? "VAR" : "STD";
      if(p->std->size>1)
        printf("  - Sky %s: %s (hdu: %s)\n", stdunit, p->usedstdname,
               p->stdhdu);
      else
        printf("  - Sky %s: %g\n", stdunit, *(float *)(p->std->array) );

      /* Convolution information. */
      if(p->convolvedname)
        printf("  - Convolved input: %s (hdu: %s)\n", p->convolvedname,
               p->chdu);
      else
        {
          if(p->kernelname)
            {
              if( strcmp(p->kernelname, UI_NO_CONV_KERNEL_NAME) )
                printf("  - Kernel: %s (hdu: %s)\n", p->kernelname, p->khdu);
              else
                printf("  - No convolution requested.\n");
            }
          else
            printf("  - Kernel: FWHM=1.5 pixel Gaussian.\n");
        }
      printf("  - Detection: %s (hdu: %s)\n", p->useddetectionname, p->dhdu);
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_abort_after_check(struct segmentparams *p, char *filename,
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
ui_free_report(struct segmentparams *p, struct timeval *t1)
{
  /* Free the allocated arrays: */
  free(p->cp.hdu);
  free(p->cp.output);
  gal_data_free(p->std);
  gal_data_free(p->input);
  gal_data_free(p->kernel);
  gal_data_free(p->binary);
  gal_data_free(p->olabel);
  gal_data_free(p->clabel);
  if(p->khdu) free(p->khdu);
  if(p->chdu) free(p->chdu);
  if(p->dhdu) free(p->dhdu);
  if(p->skyhdu) free(p->skyhdu);
  if(p->stdhdu) free(p->stdhdu);
  if(p->stdname) free(p->stdname);
  if(p->kernelname) free(p->kernelname);
  if(p->detectionname) free(p->detectionname);
  if(p->convolvedname) free(p->convolvedname);
  if(p->conv!=p->input) gal_data_free(p->conv);
  if(p->clumpsn_s_name) free(p->clumpsn_s_name);
  if(p->clumpsn_d_name) free(p->clumpsn_d_name);
  if(p->segmentationname) free(p->segmentationname);

  /* Print the final message. */
  if(!p->cp.quiet && t1)
    gal_timing_report(t1, PROGRAM_NAME" finished in: ", 0);
}
