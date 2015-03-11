/*********************************************************************
SpatialConvolve - Convolve an image in the spatial domain.
This is part of GNU Astronomy Utilities (gnuastro) package.

Copyright (C) 2013-2015 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

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

  /* Thread parameters. */
  size_t      numthreads;     /* Number of threads.               */
  size_t         *indexs;     /* Indexs to be used in this thread.*/
  pthread_barrier_t   *b;     /* Barrier to keep threads waiting. */
};




















/*******************************************************************/
/**************           Outside function          ****************/
/*******************************************************************/
void *
sconvonthread(void *inparam)
{
  struct sconvparams *scp=(struct sconvparams *)inparam;


  return NULL;
}



















/*******************************************************************/
/**************           Outside function          ****************/
/*******************************************************************/
void
spatialconvolve(float *input, size_t is0, size_t is1,
                float *kernel, size_t ks0, size_t ks1,
                size_t nt, float **out)
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
