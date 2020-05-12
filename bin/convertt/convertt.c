/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/eps.h>
#include <gnuastro/pdf.h>
#include <gnuastro/txt.h>
#include <gnuastro/fits.h>
#include <gnuastro/jpeg.h>
#include <gnuastro/arithmetic.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"
#include "color.h"










/**************************************************************/
/**************      Modifying pixel values     ***************/
/**************************************************************/
static void
convertt_change(struct converttparams *p)
{
  gal_data_t *channel, *cond;
  struct change *change, *tmp;
  unsigned char flags = ( GAL_ARITHMETIC_NUMOK
                          | GAL_ARITHMETIC_FREE
                          | GAL_ARITHMETIC_INPLACE );

  /* In case there is no value to convert. */
  if(p->change==NULL) return;

  /* Do the conversion on all channels for each change. */
  for(change=p->change; change!=NULL; change=change->next)
    for(channel=p->chll; channel!=NULL; channel=channel->next)
      {
        /* Make a condition array: all pixels with a value equal to
           'change->from' will be set as 1 in this array. */
        cond=gal_arithmetic(GAL_ARITHMETIC_OP_EQ, 1, GAL_ARITHMETIC_NUMOK,
                            channel, change->from);

        /* Now, use the condition array to set the proper values. */
        channel=gal_arithmetic(GAL_ARITHMETIC_OP_WHERE, 1, flags, channel,
                               cond, change->to);

        /* Clean up, since we set the free flag, all extra arrays have been
           freed.*/
        gal_data_free(change->from);
      }

  /* Free the channels linked list. */
  change=p->change;
  while(change!=NULL) { tmp=change->next; free(change); change=tmp; }
}




static void
convertt_trunc_function(int operator, gal_data_t *data, gal_data_t *value)
{
  gal_data_t *cond, *out;


  /* Note that we need the fluxlow and fluxhigh values later. */
  unsigned char flags = ( GAL_ARITHMETIC_NUMOK | GAL_ARITHMETIC_INPLACE );


  /* Make a condition array: all pixels with a value equal to
     'change->from' will be set as 1 in this array. */
  cond=gal_arithmetic(operator, 1, GAL_ARITHMETIC_NUMOK, data, value);


  /* Now, use the condition array to set the proper values. */
  out=gal_arithmetic(GAL_ARITHMETIC_OP_WHERE, 1, flags, data, cond, value);


  /* A small sanity check. The process must be in-place so the original
     data structure must not have changed. */
  if(out!=data)
    error(EXIT_FAILURE, 0, "%s: a bug! please contact us at %s to solve the "
          "problem. The 'out' and 'data' pointers are the same", __func__,
          PACKAGE_BUGREPORT);


  /* Clean up. */
  gal_data_free(cond);
}





static void
convertt_truncate(struct converttparams *p)
{
  gal_data_t *channel;

  /* Return if no truncation is desired. */
  if(p->fluxhigh==NULL && p->fluxlow==NULL)
    return;

  /* Do the truncation for each channel. */
  for(channel=p->chll; channel!=NULL; channel=channel->next)
    {
      if(p->fluxlow)
        convertt_trunc_function(GAL_ARITHMETIC_OP_LT, channel, p->fluxlow);
      if(p->fluxhigh)
        convertt_trunc_function(GAL_ARITHMETIC_OP_GT, channel, p->fluxhigh);
    }
}




















/**************************************************************/
/**************       convert to 8 bit          ***************/
/**************************************************************/
void
convertt_scale_to_uchar(struct converttparams *p)
{
  size_t size=p->chll->size;
  unsigned char *u, *fu, maxbyte=p->maxbyte;
  gal_data_t *channel, *prev, *copied, *mind, *maxd;
  float *f, *ff, m, tmin, tmax, min=FLT_MAX, max=-FLT_MAX;

  /* Convert everything to single precision floating point type and find
     the minimum and maximum values of all the channels in the process. */
  prev=NULL;
  for(channel=p->chll; channel!=NULL; channel=channel->next)
    {
      /* Only for channels that weren't originally blank. */
      if(channel->status==0)
        {
          /* If the type isn't float, then convert it to float. */
          if(channel->type!=GAL_TYPE_FLOAT32)
            {
              /* Change the type to float. */
              copied=gal_data_copy_to_new_type(channel, GAL_TYPE_FLOAT32);

              /* Correct the pointers. */
              copied->next=channel->next;
              if(prev) prev->next=copied; else p->chll=copied;

              /* Clean the old data structure and put in the new one */
              gal_data_free(channel);
              channel=copied;
            }

          /* Calculate the minimum and maximum. */
          mind = gal_arithmetic(GAL_ARITHMETIC_OP_MINVAL, 1, 0, channel);
          maxd = gal_arithmetic(GAL_ARITHMETIC_OP_MAXVAL, 1, 0, channel);
          tmin = *((float *)(mind->array));
          tmax = *((float *)(maxd->array));
          gal_data_free(mind);
          gal_data_free(maxd);

          /* See the over-all minimum and maximum values. */
          if(tmin<min) min=tmin;
          if(tmax>max) max=tmax;
        }

      /* Set the prev pointer. */
      prev=channel;
    }


  /* Change the minimum and maximum if desired, Note that this is only
     non-redundant when fluxhigh and fluxlow are more or less than the
     maximum and minimum values in the image.*/
  if(p->fluxlow || p->fluxhigh)
    {
      if(p->forcemin)
        {
          /* Convert the fluxlow value to float and put it in min. */
          copied=gal_data_copy_to_new_type(p->fluxlow, GAL_TYPE_FLOAT32);
          min = *((float *)(copied->array));
          gal_data_free(copied);
        }
      if(p->forcemax)
        {
          /* Convert the fluxhigh value to float and put it in min. */
          copied=gal_data_copy_to_new_type(p->fluxhigh, GAL_TYPE_FLOAT32);
          max = *((float *)(copied->array));
          gal_data_free(copied);
        }
    }
  m=(float)maxbyte/(max-min);


  /* Convert all the non-blank channels to unsigned char. */
  prev=NULL;
  for(channel=p->chll; channel!=NULL; channel=channel->next)
    {
      if(channel->status==0)
        {
          /* Convert the values into a range between '0' and 'maxbyte'. */
          ff=(f=channel->array)+size;
          if(p->invert)
            {
              do
                *f = isnan(*f) ? maxbyte : maxbyte-(*f-min)*m;
              while(++f<ff);
            }
          else
            {
              do
                *f = isnan(*f) ? 0 : (*f-min)*m;
              while(++f<ff);
            }

          /* Change the type to unsigned char. */
          copied=gal_data_copy_to_new_type(channel, GAL_TYPE_UINT8);

          /* Correct the pointers. */
          copied->next=channel->next;
          if(prev) prev->next=copied; else p->chll=copied;

          /* Clean the old data structure and put in the new one */
          gal_data_free(channel);
          channel=copied;
        }
      else
        {
          /* In CMYK, a blank channel should have a maximum value. */
          if(p->numch==4)
            {fu=(u=channel->array)+size; do *u=UINT8_MAX; while(++u<fu);}
        }

      /* Set the prev pointer. */
      prev=channel;
    }
}




















/**************************************************************/
/**************           Main function         ***************/
/**************************************************************/
void
convertt(struct converttparams *p)
{
  gal_data_t *channel;

  /* Make any of the desired changes to the data. */
  if(p->changeaftertrunc)
    {
      convertt_truncate(p);
      convertt_change(p);
    }
  else
    {
      convertt_change(p);
      convertt_truncate(p);
    }

  /* Save the outputs: */
  switch(p->outformat)
    {
    /* FITS: a FITS file can have many extensions (channels). */
    case OUT_FORMAT_FITS:
      if(p->numch==3 && p->rgbtohsv)
        color_rgb_to_hsv(p);
      for(channel=p->chll; channel!=NULL; channel=channel->next)
        gal_fits_img_write(channel, p->cp.output, NULL, PROGRAM_NAME);
      break;

    /* Plain text: only one channel is acceptable. */
    case OUT_FORMAT_TXT:
      gal_checkset_writable_remove(p->cp.output, 0, p->cp.dontdelete);
      gal_txt_write(p->chll, NULL, p->cp.output, 0);
      break;

    /* JPEG: */
    case OUT_FORMAT_JPEG:
      if(p->colormap) color_map_prepare(p); else convertt_scale_to_uchar(p);
      gal_jpeg_write(p->chll, p->cp.output, p->quality, p->widthincm);
      break;

    /* EPS. */
    case OUT_FORMAT_EPS:
      if(p->colormap) color_map_prepare(p); else convertt_scale_to_uchar(p);
      gal_eps_write(p->chll, p->cp.output, p->widthincm, p->borderwidth,
                    p->hex, p->forcemin || p->forcemax, 0);
      break;

    /* PDF */
    case OUT_FORMAT_PDF:
      if(p->colormap) color_map_prepare(p); else convertt_scale_to_uchar(p);
      gal_pdf_write(p->chll, p->cp.output, p->widthincm, p->borderwidth,
                    p->forcemin || p->forcemax);
      break;

    /* Not recognized. */
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us so we can find "
            "the problem and fix it The internal type of the output is "
            "not recognized. ", __func__);
    }
}
