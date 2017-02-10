/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <string.h>

#include <gnuastro/txt.h>
#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/table.h>
#include <gnuastro/arithmetic.h>

#include <timing.h>
#include <options.h>
#include <checkset.h>
#include <fixedstringmacros.h>

#include "main.h"

#include "ui.h"
#include "eps.h"
#include "jpeg.h"
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
args_doc[] = "InputFile1 [InputFile2] ... [InputFile4]";

const char
doc[] = GAL_STRINGS_TOP_HELP_INFO PROGRAM_NAME" will convert any of the "
  "known input formats to any other of the known formats. The output file "
  "will have the same number of pixels.\n"
  GAL_STRINGS_MORE_HELP_INFO
  /* After the list of options: */
  "\v"
  PACKAGE_NAME" home page: "PACKAGE_URL;





/* Option groups particular to this program */
enum program_args_groups
{
  ARGS_GROUP_FLUX = GAL_OPTIONS_GROUP_AFTER_COMMON,
};




/* Available letters for short options:

   d e f g j k l n p r s t v y z
   A B E F G I J M O Q R T U W X Y Z      */
enum option_keys_enum
{
  /* With short-option version. */
  ARGS_OPTION_KEY_QUALITY             = 'u',
  ARGS_OPTION_KEY_WIDTHINCM           = 'w',
  ARGS_OPTION_KEY_BORDERWIDTH         = 'b',
  ARGS_OPTION_KEY_HEX                 = 'x',
  ARGS_OPTION_KEY_FLUXLOW             = 'L',
  ARGS_OPTION_KEY_FLUXHIGH            = 'H',
  ARGS_OPTION_KEY_HIGH                = 'H',
  ARGS_OPTION_KEY_MAXBYTE             = 'm',
  ARGS_OPTION_KEY_FLMINBYTE           = 'a',
  ARGS_OPTION_KEY_FHMAXBYTE           = 'b',
  ARGS_OPTION_KEY_CHANGE              = 'c',
  ARGS_OPTION_KEY_CHANGEAFTERTRUNC    = 'C',
  ARGS_OPTION_KEY_INVERT            = 'i',

  /* Only with long version (start with a value 1000, the rest will be set
     automatically). */
};



















/**************************************************************/
/*********    Initialize & Parse command-line    **************/
/**************************************************************/
static void
ui_initialize_options(struct converttparams *p,
                      struct argp_option *program_options,
                      struct argp_option *gal_commonopts_options)
{
  size_t i;
  struct gal_options_common_params *cp=&p->cp;


  /* Set the necessary common parameters structure. */
  cp->poptions           = program_options;
  cp->program_name       = PROGRAM_NAME;
  cp->program_exec       = PROGRAM_EXEC;
  cp->program_bibtex     = PROGRAM_BIBTEX;
  cp->program_authors    = PROGRAM_AUTHORS;
  cp->coptions           = gal_commonopts_options;

  /* Program specific non-zero values. */
  p->maxbyte             = UINT8_MAX;

  /* Edit the common options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    switch(cp->coptions[i].key)
      {
      case GAL_OPTIONS_KEY_HDU:
        cp->coptions[i].value=&p->hdus;
        cp->coptions[i].type=GAL_DATA_TYPE_STRLL;
        cp->coptions[i].doc="FITS input HDU, multiple calls possible.";
        break;

      case GAL_OPTIONS_KEY_OUTPUT:
        cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
        cp->coptions[i].doc="Output filename or suffix.";
        break;

      case GAL_OPTIONS_KEY_MINMAPSIZE:
        cp->coptions[i].mandatory=GAL_OPTIONS_MANDATORY;
        break;

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
  struct converttparams *p = state->input;

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
      gal_linkedlist_add_to_stll(&p->inputnames, arg, 0);
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
   check in `ui_check_options_and_arguments'. */
static void
ui_read_check_only_options(struct converttparams *p)
{
  gal_data_t *cond;


  /* Check if quality is smaller and equal to 100 (we already checked if it
     was greater than zero before. */
  if(p->quality>100)
    error(EXIT_FAILURE, 0, "`%u' is larger than 100. The value to the "
          "`--quality' (`-u') option must be between 1 and 100 (inclusive)",
          p->quality);


  /* Read the truncation values into a data structure and see if flux low
     is indeed smaller than fluxhigh. */
  if(p->fluxlowstr)
    {
      p->fluxlow=gal_data_string_to_number(p->fluxlowstr);
      if(p->fluxlow==NULL)
        error(EXIT_FAILURE, 0, "value to the `--fluxlow' (`-L', %s) "
              "couldn't be read as a number", p->fluxlowstr);
    }

  if(p->fluxhighstr)
    {
      p->fluxhigh=gal_data_string_to_number(p->fluxhighstr);
      if(p->fluxhigh==NULL)
        error(EXIT_FAILURE, 0, "value to the `--fluxhigh' (`-H', %s) "
              "couldn't be read as a number", p->fluxhighstr);
    }

  if(p->fluxhighstr && p->fluxlowstr)
    {
      cond=gal_arithmetic(GAL_ARITHMETIC_OP_GT, GAL_ARITHMETIC_NUMOK,
                          p->fluxhigh, p->fluxlow);

      if( *((unsigned char *)cond->array) == 0 )
        error(EXIT_FAILURE, 0, "The value of `--fluxlow' must be less "
              "than `--fluxhigh'");

      gal_data_free(cond);
    }
}





static void
ui_check_options_and_arguments(struct converttparams *p)
{
  /* Check if there was any inputs. */
  if(p->inputnames==NULL)
    error(EXIT_FAILURE, 0, "no input files given");

  /* Reverse the `inputnames' linked list, note that the `hdu' linked list
     was reversed during option parsing.*/
  gal_linkedlist_reverse_stll(&p->inputnames);

}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
static struct change *
ui_make_change_struct(char *arg)
{
  char *p=arg;
  gal_data_t *data;
  size_t len=0, counter=0;
  struct change *out, *last=NULL, *ch;

  /* First set all the delimiters to `\0' and count the number of
     characters in the full string. */
  while(*p!='\0')
    {
      if( isspace(*p) || *p==':' || *p==',' ) *p='\0';
      ++p;
    }
  len=p-arg;

  /* Now, go through the string and read everything that remains. */
  p=arg;
  while(p<arg+len)
    {
      if(*p=='\0')
        ++p;
      else
        {
          /* Read the number and increment the counter. */
          ++counter;
          data=gal_data_string_to_number(p);
          if(data==NULL)
            error(EXIT_FAILURE, 0, "`%s' (input number %zu to the "
                  "`--change' option) couldn't be read as a number", p,
                  counter);

          /* Go to the end of this number (until you reach a `\0'). */
          while(*p!='\0') {++p; continue;}

          /* Put the data structure in the correct place. When the counter is
             an odd number, we have just started a new set of changes.*/
          if(counter%2)             /* Odd. */
            {
              /* Allocate space for the new structure. */
              errno=0;
              ch=malloc(sizeof *ch);
              if(ch==NULL)
                error(EXIT_FAILURE, errno, "%zu bytes for ch in "
                      "`ui_make_change_struct'", sizeof *ch);

              /* If the last structure has already been defined (!=NULL)
                 then we should set its next element to `ch' and change it
                 to point to `ch'. On the other hand, when this is the
                 first structure to be created, then `last==NULL', so to
                 start off the process, we should put `ch' into both the
                 `out' and `last' lists.. */
              if(last)
                {
                  last->next=ch;
                  last=ch;
                }
              else
                out=last=ch;

              /* Put `data' in the `from' element, and since this is the
                 last structure, set its next element to NULL. */
              last->from=data;
              last->next=NULL;
            }
          else                      /* Even. */
            last->to=data;
        }
    }

  /*
    {
    struct change *tmp;
    for(tmp=out;tmp!=NULL;tmp=tmp->next)
    printf("%f --> %f\n", tmp->from, tmp->to);
    }
  */
  return out;
}





/* Go through the input files and make a linked list of all the channels
   that exist in them. When this function finishes the list of channels
   will be filled in the same order as they were read from the inputs. */
static void
ui_make_channels_ll(struct converttparams *p)
{
  char *hdu;
  size_t dsize=0;
  gal_data_t *data;
  struct gal_linkedlist_stll *name;

  /* Go through the input files and add the channel(s). */
  p->numch=0;
  for(name=p->inputnames; name!=NULL; name=name->next)
    {
      /* Check if p->numch has not exceeded 4. */
      if(p->numch>=4)
        error(EXIT_FAILURE, 0, "the number of input color channels (not "
              "necessarily files) has exceeded 4! Note that one file can "
              "contain more than one color channel (for example a JPEG "
              "file in RGB has 3 channels)");

      /* Make sure this input file exists (if it isn't blank). */
      if(strcmp(name->v, "blank")) gal_checkset_check_file(name->v);

      /* FITS: */
      if( gal_fits_name_is_fits(name->v) )
        {
          /* Get the HDU value for this channel. */
          if(p->hdus)
            gal_linkedlist_pop_from_stll(&p->hdus, &hdu);
          else
            error(EXIT_FAILURE, 0, "not enough HDUs. Every input FITS image "
                  "needs a HDU, you can use the `--hdu' (`-h') option once "
                  "for each input FITS image (in the same order)");

          /* Read in the array and its WCS information. */
          data=gal_fits_img_read(name->v, hdu, p->cp.minmapsize);
          gal_wcs_read(name->v, hdu, 0, 0, &data->nwcs, &data->wcs);
          gal_data_add_existing_to_ll(&p->chll, data);

          /* A FITS file only has one channel. */
          ++p->numch;
        }



      /* JPEG: */
      else if ( nameisjpeg(name->v) )
        {
#ifndef HAVE_LIBJPEG
          error(EXIT_FAILURE, 0, "you are giving a JPEG input, however, "
                "when %s was configured libjpeg was not available. To read "
                "from (and write to) JPEG files, libjpeg is required. "
                "Please install it and configure, make and install %s "
                "again", PACKAGE_STRING, PACKAGE_STRING);
#else
          p->numch += jpeg_read_to_ll(name->v, &p->chll, p->cp.minmapsize);
#endif
        }



      /* Blank: */
      else if(strcmp(name->v, BLANK_CHANNEL_NAME)==0)
        {
          gal_data_add_to_ll(&p->chll, NULL, GAL_DATA_TYPE_INVALID, 0,
                             &dsize, NULL, 0, p->cp.minmapsize, "blank",
                             NULL, NULL);
          ++p->numch;
        }



      /* EPS:  */
      else if ( nameiseps(name->v) )
        error(EXIT_FAILURE, 0, "EPS files cannot be used as input. Since "
              "EPS files are not raster graphics, they are only used as "
              "output");



      /* PDF:  */
      else if ( nameispdf(name->v) )
        error(EXIT_FAILURE, 0, "PDF files cannot be used as input. Since "
              "PDF files are not raster graphics, they are only used as "
              "output");


      /* Text: */
      else
        {
          data=gal_txt_image_read(name->v, p->cp.minmapsize);
          gal_data_add_existing_to_ll(&p->chll, data);
          ++p->numch;
        }
    }


  /* Reverse the list of channels into the input order. */
  gal_data_reverse_ll(&p->chll);
}





/* Read the input(s)/channels. */
static void
ui_prepare_input_channels(struct converttparams *p)
{
  struct wcsprm *wcs=NULL;
  size_t i, ndim, *dsize=NULL;
  gal_data_t *tmp, *blank, *prev;

  /* Fill in the channels linked list. */
  ui_make_channels_ll(p);


  /* Make sure there are 1 (for grayscale), 3 (for RGB) or 4 (for CMYK)
     color channels. */
  if(p->numch!=1 && p->numch!=3 && p->numch!=4)
    error(EXIT_FAILURE, 0, "the number of input color channels has to be "
          "1 (for non image data, grayscale or only K channel in CMYK), "
          "3 (for RGB) and 4 (for CMYK). You have given %zu color channels. "
          "Note that some file formats (for example JPEG in RGB mode) can "
          "contain more than one color channel", p->numch);


  /* Go over the channels and make the proper checks/corrections. We won't
     be checking blank channels here, recall that blank channels had a
     dimension of zero. */
  for(tmp=p->chll; tmp!=NULL; tmp=tmp->next)
    if(tmp->ndim>0)
      {
        /* Set the reference size (to check and also use for the blank
           channels). */
        if(dsize==NULL)
          {
            ndim=tmp->ndim;
            dsize=tmp->dsize;
          }
        else
          {
            if(tmp->ndim!=ndim)
              error(EXIT_FAILURE, 0, "All channels must have the same "
                    "number of dimensions, the first input channel had "
                    "%zu dimensions while atleast one other has %zu",
                    ndim, tmp->ndim);
            for(i=0;i<ndim;++i)
              if(dsize[i]!=tmp->dsize[i])
                error(EXIT_FAILURE, 0, "The length along each dimension of "
                      "the channels must be the same");
          }

        /* Incase there is WCS information, also keep a pointer to the
           first WCS information encountered. */
        if(wcs==NULL && tmp->wcs)
          wcs=tmp->wcs;
      }


  /* If ndim is still NULL, then there were no non-blank inputs, so print
     an error. */
  if(dsize==NULL)
    error(EXIT_FAILURE, 0, "all the input(s) are of type blank");


  /* Now, fill in the blank channels with zero valued arrays. */
  prev=NULL;
  for(tmp=p->chll; tmp!=NULL; tmp=tmp->next)
    {
      /* If this is a blank structure, then set it to a zero valued
         array. */
      if(tmp->ndim==0)
        {
          /* Make the blank data structure. */
          blank=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, ndim, dsize,
                               wcs, 1, p->cp.minmapsize, "blank channel",
                               NULL, NULL);

          /* We will use the status value of the data structuer to mark it
             as one that was originally blank. */
          blank->status=1;

          /* If a previous node pointed to this old blank structure, then
             correct it. */
          if(prev) prev->next=blank; else p->chll=blank;

          /* Set the next pointer of this one to same pointer that the old
             blank pointer pointed to. */
          blank->next=tmp->next;

          /* Free the old data structure and put this one in its place. */
          gal_data_free(tmp);
          tmp=blank;
        }

      /* This is the final (to be used) data structure, so keep its pointer
         in case the next one is blank and this structure's `next' element
         must be corrected. */
      prev=tmp;
    }
}





/* We know cp->output is a known suffix, we just don't know if it has a `.`
   before it or not. If it doesn't, one will be added to it and the output
   name will be set using the automatic output function. */
void
ui_add_dot_use_automatic_output(struct converttparams *p)
{
  struct gal_linkedlist_stll *stll;
  char *tmp, *firstname="output.txt", *suffix=p->cp.output;

  /* Find the first non-blank file name in the input(s). */
  for(stll=p->inputnames; stll!=NULL; stll=stll->next)
    if(strcmp(stll->v, BLANK_CHANNEL_NAME))
      {
        firstname=stll->v;
        break;
      }

  /* If the suffix does not start with a `.', put one there. */
  if(suffix[0]!='.')
    {
      asprintf(&tmp, ".%s", suffix);
      free(suffix);
      suffix=tmp;
    }

  /* Set the automatic output and make sure we have write access. */
  p->cp.output=gal_checkset_automatic_output(&p->cp, firstname, suffix);
}





/* Set output name, not that for ConvertType, the output option value is
   mandatory (in `args.h'). So by the time the program reaches here, we
   know it exists. */
static void
ui_set_output(struct converttparams *p)
{
  struct gal_options_common_params *cp=&p->cp;

  /* Determine the type */
  if(gal_fits_name_is_fits(cp->output))
    {
      p->outformat=OUT_FORMAT_FITS;
      if( gal_fits_suffix_is_fits(cp->output) )
        ui_add_dot_use_automatic_output(p);
    }
  else if(nameisjpeg(cp->output))
    {
#ifndef HAVE_LIBJPEG
      error(EXIT_FAILURE, 0, "you have asked for a JPEG output, "
            "however, when %s was configured libjpeg was not "
            "available. To write to JPEG files, libjpeg is required. "
            "Please install it and configure, make and install %s "
            "again", PACKAGE_STRING, PACKAGE_STRING);
#else
      p->outformat=OUT_FORMAT_JPEG;
      if( nameisjpegsuffix(cp->output) )
        ui_add_dot_use_automatic_output(p);
      if(p->quality==0)
        error(EXIT_FAILURE, 0, "no quality specified for %s, please use the "
              "`--quality' (`-q') option with a value between 1 (low "
              "quality) and 100 (high quality) to specify a level",
              p->cp.output);
#endif
    }
  else if(nameiseps(cp->output))
    {
      p->outformat=OUT_FORMAT_EPS;
      if( nameisepssuffix(cp->output) )
        ui_add_dot_use_automatic_output(p);
    }
  else if(nameispdf(cp->output))
    {
      p->outformat=OUT_FORMAT_PDF;
      if( nameispdfsuffix(cp->output) )
        ui_add_dot_use_automatic_output(p);
    }
  else
    {
      p->outformat=OUT_FORMAT_TXT;

      /* If the given value is `stdout', then set p->cp.output to NULL, so
         the result will be printed to the standard output. */
      if( !strcmp(p->cp.output, "stdout") )
        {
          free(p->cp.output);
          p->cp.output=NULL;
        }
      else
        {
          /* Plain text files don't have any unique set of suffixes. So,
             here, we will just adopt two of the most common ones: `txt' or
             `dat'. If the output is just one of these two suffixes, then
             we will use automatic output to generate the full name,
             otherwise, we'll just take the user's given value as the
             filename. */
          if( !strcmp(cp->output, "txt") || !strcmp(cp->output, ".txt")
              || !strcmp(cp->output, "dat") || !strcmp(cp->output, ".dat") )
            ui_add_dot_use_automatic_output(p);

          /* If output type is not an image, there should only be one color
             channel: */
          if(p->numch>1)
            error(EXIT_FAILURE, 0, "text output (`--output=%s`) can only be "
                  "completed with one input color channel. You have given "
                  "%zu. Note that some formats (for example JPEG) can have "
                  "more than one color channel in each file. You can first "
                  "convert the file to FITS, then convert the desired "
                  "channel to text by specifying the HDU",
                  cp->output, p->numch);
        }
    }

  /* Check if the output already exists. */
  gal_checkset_check_remove_file(cp->output, cp->dontdelete);
}





void
ui_preparations(struct converttparams *p)
{
  /* Convert the change string into the proper list. */
  if(p->changestr)
    p->change=ui_make_change_struct(p->changestr);

  /* Read the input channels. */
  ui_prepare_input_channels(p);

  /* Set the output name. */
  ui_set_output(p);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/

void
ui_read_check_inputs_setup(int argc, char *argv[], struct converttparams *p)
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


  /* Read the options into the program's structure, and check them and
     their relations prior to printing. */
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
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
ui_free_report(struct converttparams *p)
{
  /* Free the allocated spaces */
  gal_data_free(p->fluxlow);
  gal_data_free(p->fluxhigh);
  if(p->cp.output) free(p->cp.output);
  gal_linkedlist_free_stll(p->hdus, 1);
  gal_linkedlist_free_stll(p->inputnames, 0);
}
