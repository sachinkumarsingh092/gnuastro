/*********************************************************************
ConvertType - Convert between various types of files.
ConvertType is part of GNU Astronomy Utilities (gnuastro) package.

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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include "timing.h"
#include "statistics.h"
#include "txtarrayvv.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "jpeg.h"
#include "eps.h"


/**************************************************************/
/**************      Modifying pixel values     ***************/
/**************************************************************/
void
changevalue(struct converttparams *p)
{
  struct change *tmp, *ttmp;
  size_t numchange=0, i, j, size;
  double *from, *to, **ch=p->ch, *d, *df, f, t;

  /* In case there is no value to convert. */
  if(p->change==NULL) return;

  /* Count the number of conversions, allocate the arrays for them
     and put the conversion values in it. */
  for(tmp=p->change;tmp!=NULL;tmp=tmp->next) ++numchange;
  errno=0;
  to=malloc(numchange*sizeof *to);
  from=malloc(numchange*sizeof *from);
  if(to==NULL || from==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for `to` or `from` in convert.",
          numchange*sizeof *from);
  i=1;
  tmp=p->change;
  while(tmp!=NULL)
    {
      to[numchange-i]=tmp->to;
      from[numchange-i]=tmp->from;
      ttmp=tmp->next;
      free(tmp);
      tmp=ttmp;
      ++i;
    }

  /* Convert the pixel values in the given order: */
  size=p->s0[0]*p->s1[0];
  for(i=0;i<p->numch;++i)
    {
      if(p->isblank[i]) continue;
      for(j=0;j<numchange;++j)
        {
          f=from[j]; t=to[j];
          df=(d=ch[i])+size;
          do *d = (*d==f) ? t : *d; while(++d<df);
        }
    }

  free(to);
  free(from);
}





void
truncateflux(struct converttparams *p)
{
  size_t i, size;
  double *d, *df, **ch=p->ch;
  double fluxlow=p->fluxlow, fluxhigh=p->fluxhigh;

  /* If the flux truncation values are equal, there is no truncation.*/
  if(fluxlow==fluxhigh) return;

  size=p->s0[0]*p->s1[0];
  for(i=0;i<p->numch;++i)
    {
      if(p->isblank[i]) continue;
      df=(d=ch[i])+size;
      do
        {
          if(*d<fluxlow) *d=fluxlow;
          if(*d>fluxhigh) *d=fluxhigh;
        }
      while(++d<df);
    }
}





/* Change the array to the logarithm of the array.

   Note that the logarithm of a zero or negative value is not
   defined. We assume here that the full range of flux is desired. So
   if the minimum value in the array is zero or negative, all the
   elements of the array are added to a constant such that the minimum
   value is only slightly larger than zero. The "slight"ness is
   defined in terms of the distance between the minimum and maximum
   for each color channel.*/
void
takelog(struct converttparams *p)
{
  size_t i, size;
  double **ch=p->ch;
  double min, max, toadd, *d, *df;

  size=p->s0[0]*p->s1[0];
  for(i=0;i<p->numch;++i)
    {
      dminmax(ch[i], size, &min, &max);
      if(p->isblank[i]) continue;
      if(min<=0.0f)
        {
          toadd=-1.0*min+(max-min)/10000.0f;
          df=(d=ch[i])+size;
          do *d+=toadd; while(++d<df);
        }
      df=(d=ch[i])+size;
      do *d=log10(*d); while(++d<df);
    }
}




















/**************************************************************/
/**************       Save text and FITS        ***************/
/**************************************************************/
void
savetxt(struct converttparams *p)
{
  char comments[10000];
  int iprec[]={0, 8}, fprec[]={6, 8};
  int int_cols[]={-1}, accu_cols[]={-1};
  int ispace[]={1, 10, 15}, fspace[]={1, 15, 15};


  /* Find the input file name and  */
  sprintf(comments, "# Pixel values of %s (%lu x %lu pixels).\n"
          "# Created on %s# %s", p->names[0], p->s0[0], p->s1[0],
          ctime(&p->rawtime), SPACK_STRING);

  if(p->bitpixs[0]==BYTE_IMG || p->bitpixs[0]==SHORT_IMG
     || p->bitpixs[0]==LONG_IMG || p->bitpixs[0]==LONGLONG_IMG )
    arraytotxt(p->ch[0], p->s0[0], p->s1[0], comments, int_cols,
               accu_cols, ispace, iprec, 'f', p->cp.output);
  else
    arraytotxt(p->ch[0], p->s0[0], p->s1[0], comments, int_cols,
               accu_cols, fspace, fprec, 'g', p->cp.output);

}





void
savefits(struct converttparams *p)
{
  void *array;
  size_t i, size;
  char hdu[1000];

  size=p->s0[0]*p->s1[0];
  for(i=0;i<p->numch;++i)
    {
      /* Make sure array is in the correct format. */
      if(p->bitpixs[i]!=DOUBLE_IMG)
        changetype(p->ch[i], DOUBLE_IMG, size, p->numnul[i],
                   &array, p->bitpixs[i]);
      else
        array=p->ch[i];

      /* Write array to a FITS file.*/
      sprintf(hdu, "Channel%lu", i+1);
      arraytofitsimg(p->cp.output, hdu, p->bitpixs[i], array,
                     p->s0[i], p->s1[i], 0, NULL, SPACK_STRING);

      /* If array was allocated separately, free it. */
      if(p->bitpixs[i]!=DOUBLE_IMG)
        free(array);
    }
}




















/**************************************************************/
/**************       convert to 8 bit          ***************/
/**************************************************************/
void
doubleto8bit(struct converttparams *p)
{
  size_t i, size, numch=p->numch;
  uint8_t *u, *fu, **ech=p->ech, maxbyte=p->maxbyte;
  double *d, m, tmin, tmax, min=FLT_MAX, max=-FLT_MAX;

  /* Allocate space for all the channels. */
  size=p->s0[0]*p->s1[0];
  for(i=0;i<numch;++i)
    {
      errno=0;
      ech[i]=malloc(size*sizeof *ech[i]);
      if(ech[i]==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for ech[%lu]",
              size*sizeof *ech[i], i);
    }

  /* Find the minimum and maximum of all the color channels. */
  for(i=0;i<numch;++i)
    {
      dminmax(p->ch[i], size, &tmin, &tmax);
      if(i)
        {
          if(tmin<min) min=tmin;
          if(tmax>max) max=tmax;
          continue;
        }
      min=tmin;
      max=tmax;
    }

  /* Change the minimum and maximum if desired, Note that this
     is only non-redundant when fluxhigh and fluxlow are more
     or less than the maximum and minimum values in the
     image.*/
  if(p->fluxlow!=p->fluxhigh)
    {
      if(p->flminbyte) max=p->fluxlow;
      if(p->fhmaxbyte) max=p->fluxhigh;
    }
  m=(double)maxbyte/(max-min);

  /* Write the values in the array: */
  for(i=0;i<numch;++i)
    {
      if(p->isblank[i])
        {
          if(numch==3)          /* RGB  */
            {fu=(u=ech[i])+size; do *u=0; while(++u<fu);}
          else if(numch==4)     /* CMYK */
            {fu=(u=ech[i])+size; do *u=UINT8_MAX; while(++u<fu);}
          else error(EXIT_FAILURE, 0, "A bug! The number of channels in "
                     "doubleto8bit is not 3 or 4 when there is a blank "
                     "channel, this should not happen. Please contact us "
                     "So we can fix it.");
        }
      else
        {
          d=p->ch[i];
          fu=(u=ech[i])+size;
          if(p->invert)
            {
              do
                {
                  *u = isnan(*d) ? maxbyte : maxbyte-(*d-min)*m;
                  ++d;
                }
              while(++u<fu);
            }
          else
            {
              do
                {
                  *u = isnan(*d) ? 0 : (*d-min)*m;
                  ++d;
                }
              while(++u<fu);
            }
        }
    }
  /* Check all channels (practical for a small input file).
  {
    size_t j;
    for(j=0;j<size;++j)
      {
        for(i=0;i<numch;++i)
          printf("%d:", (int)p->ech[i][j]);
        printf("\b.\n");
      }
  }
  */
}



















/**************************************************************/
/**************           Main function         ***************/
/**************************************************************/
void
convertt(struct converttparams *p)
{
  size_t i;

  /* Make any of the desired changes to the data. */
  if(p->changeaftertrunc)
    {
      truncateflux(p);
      changevalue(p);
    }
  else
    {
      changevalue(p);
      truncateflux(p);
    }
  if(p->log) takelog(p);



  /* Save the outputs: */
  switch(p->outputtype)
    {
    case TXTFORMAT:
      savetxt(p);
      break;
    case FITSFORMAT:
      savefits(p);
      break;
    case JPEGFORMAT:
      doubleto8bit(p);
#ifdef HAS_LIBJPEG
      savejpeg(p);
#else
      error(EXIT_FAILURE, 0, "You have asked for a JPEG output, however, "
            "when %s was configured libjpeg was not available. To write "
            "to JPEG files, libjpeg is required. Please install it and "
            "configure, make and install %s again.", PACKAGE_STRING,
            PACKAGE_STRING);
#endif
      break;
    case EPSFORMAT: case PDFFORMAT:
      doubleto8bit(p);
      saveepsorpdf(p);
      break;
    default:
      error(EXIT_FAILURE, 0, "A bug! The internal type of the output is "
            "not recognized. Please contact us so we can find the problem "
            "and fix it.");
    }

  /* Free the ech arrays. Note that if they have not been
     allocated, they are NULL and free won't complain. */
  for(i=0;i<p->numch;++i)
    free(p->ech[i]);
}
