/*********************************************************************
SpatialConvolve - Convolve an image in the spatial domain.
This is part of GNU Astronomy Utilities (gnuastro) package.

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
#include <stdlib.h>

#include "box.h"
#include "astrthreads.h"



/*******************************************************************/
/**************           Input structure           ****************/
/*******************************************************************/
struct sconvparams
{
  /* General input parameters: */
  float           *input;     /* Input image array.               */
  float          *kernel;     /* Kernel image array.              */
  float             *out;     /* Output image.                    */
  size_t             is0;     /* Image size along first C axis.   */
  size_t             is1;     /* Image size along second C axis.  */
  size_t             ks0;     /* Kernel size along first C axis.  */
  size_t             ks1;     /* Kernel size along second C axis. */
  float      blankweight;     /* Maximum blank pixel weight.      */
  int        inputhasnan;     /* Input image has NAN pixels.      */

  /* Thread parameters. */
  size_t      numthreads;     /* Number of threads.               */
  size_t         *indexs;     /* Indexs to be used in this thread.*/
  pthread_barrier_t   *b;     /* Barrier to keep threads waiting. */
};




















/*******************************************************************/
/**************             Each thread             ****************/
/*******************************************************************/
void
nonusedpixels(struct sconvparams *scp)
{

}





void *
sconvonthread(void *inparam)
{
  struct sconvparams *scp=(struct sconvparams *)inparam;

  double sum;
  size_t numrows;
  int inputhasnan=scp->inputhasnan;
  long naxes[2]={scp->is1, scp->is0};
  float blankweight, *istart, *kstart;
  float *f, *fp, *k, maxbw=scp->blankweight;
  size_t i, j, os1, ind, *indexs=scp->indexs;
  size_t is1=scp->is1, ks0=scp->ks0, ks1=scp->ks1;
  long fpixel_i[2], lpixel_i[2], fpixel_o[2], lpixel_o[2];
  float *input=scp->input, *kernel=scp->kernel, *out=scp->out;

  /* Go over all the pixels associated with this thread. */
  for(i=0;indexs[i]!=NONTHRDINDEX;++i)
    {
      /* Set the starting and ending pixels on the kernel, note that
         the overlap function in box.c uses FITS coordinates. */
      ind=indexs[i];
      fpixel_o[0]=1;                fpixel_o[1]=1;
      lpixel_o[0]=ks1;              lpixel_o[1]=ks0;
      fpixel_i[0]=ind%is1+1-ks1/2;  fpixel_i[1]=ind/is1+1-ks0/2;
      lpixel_i[0]=ind%is1+1+ks1/2;  lpixel_i[1]=ind/is1+1+ks0/2;
      overlap(naxes, fpixel_i, lpixel_i, fpixel_o, lpixel_o);

      /* fpixels and lpixels now point to the overlap's starting and
         ending both on the image and on the kernel. */
      sum=0.0f;
      blankweight=0.0f;
      os1=lpixel_i[0]-fpixel_i[0];
      numrows=lpixel_i[1]-fpixel_i[1]+1; /* fpixel is inside the box */
      istart= &input[ (fpixel_i[1]-1) * is1 + fpixel_i[0]-1 ];
      kstart=&kernel[ (fpixel_o[1]-1) * ks1 + fpixel_o[0]-1 ];
      for(j=0;j<numrows;++j)
        {
          k=kstart+j*ks1;
          fp = ( f=istart+j*is1 ) + os1;
          do
            {
              if( inputhasnan && isnan(*f) )
                {
                  blankweight+=*k;
                  if(blankweight>maxbw)
                    {
                      out[ind]=NAN;
                      break;
                    }
                }
              else
                sum += *k * *f;
              ++k;
            }
          while(++f<fp);
          if isnan(out[ind]) break;
          out[ind]=sum;
        }
    }

  /* Wait until all other threads finish. */
  if(scp->numthreads>1)
    pthread_barrier_wait(scp->b);

  return NULL;
}



















/*******************************************************************/
/**************           Outside function          ****************/
/*******************************************************************/
void
spatialconvolve(float *input, size_t is0, size_t is1,
                float *kernel, size_t ks0, size_t ks1,
                size_t nt, float blankweight, int inputhasnan,
                float **out)
{
  int err;
  pthread_t t;          /* All thread ids saved in this, not used. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct sconvparams *scp;
  size_t i, nb, *indexs, thrdcols;


  /* Allocate the arrays to keep the thread and parameters for each
     thread. */
  errno=0;
  scp=malloc(nt*sizeof *scp);
  if(scp==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes in spatialconvolve "
          "(spatialconvolve.c) for scp", nt*sizeof *scp);


  /* Allocate the output array: */
  errno=0;
  *out=malloc(is0*is1*sizeof **out);
  if(*out==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for convolution output.",
          is0*is1*sizeof **out);


  /* Distribute the image pixels into the threads: */
  distinthreads(is0*is1, nt, &indexs, &thrdcols);

  if(nt==1)
    {
      scp[0].input=input; scp[0].kernel=kernel; scp[0].out=*out;
      scp[0].is0=is0; scp[0].is1=is1; scp[0].ks0=ks0;
      scp[0].ks1=ks1; scp[0].indexs=indexs; scp[0].numthreads=nt;
      scp[0].blankweight=blankweight; scp[0].inputhasnan=inputhasnan;
      sconvonthread(&scp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
	 (that spinns off the nt threads) is also a thread, so the
	 number the barrier should be one more than the number of
	 threads spinned off. */
      if(is0*is1<nt) nb=is0*is1+1;
      else nb=nt+1;
      attrbarrierinit(&attr, &b, nb);

      /* Spin off the threads: */
      for(i=0;i<nt;++i)
        if(indexs[i*thrdcols]!=NONTHRDINDEX)
          {
            scp[i].b=&b;
            scp[i].inputhasnan=inputhasnan;
            scp[i].blankweight=blankweight;
            scp[i].indexs=&indexs[i*thrdcols];
            scp[i].input=input; scp[i].kernel=kernel;
            scp[i].out=*out; scp[i].is0=is0; scp[i].is1=is1;
            scp[i].ks0=ks0; scp[i].ks1=ks1; scp[i].numthreads=nt;
	    err=pthread_create(&t, &attr, sconvonthread, &scp[i]);
	    if(err)
	      error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }

  /* Free the allocated spaces: */
  free(scp);
  free(indexs);
}
