/*********************************************************************
SpatialConvolve - Convolve an image in the spatial domain.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/box.h>
#include <gnuastro/spatialconvolve.h>










/*******************************************************************/
/**************             Each thread             ****************/
/*******************************************************************/
void
gal_spatialconvolve_pparams(float *input, size_t is0, size_t is1,
                            float *kernel, size_t ks0, size_t ks1,
                            size_t nt, int edgecorrection, float *out,
                            size_t *indexs,
                            struct gal_spatialconvolve_params *scp)
{
  /* Put the simple values in: */
  scp->is0=is0;
  scp->is1=is1;
  scp->ks0=ks0;
  scp->ks1=ks1;
  scp->out=out;
  scp->input=input;
  scp->kernel=kernel;
  scp->indexs=indexs;
  scp->numthreads=nt;
  scp->edgecorrection=edgecorrection;
}






void *
gal_spatialconvolve_thread(void *inparam)
{
  struct gal_spatialconvolve_params *scp=
    (struct gal_spatialconvolve_params *)inparam;

  double sum, ksum;
  long naxes[2]={scp->is1, scp->is0};
  float *f, *fp, *k, *istart, *kstart;
  int edgecorrection=scp->edgecorrection;
  size_t is1=scp->is1, ks0=scp->ks0, ks1=scp->ks1;
  long *fpixel_i=scp->fpixel_i, *lpixel_i=scp->lpixel_i;
  long *fpixel_o=scp->fpixel_o, *lpixel_o=scp->lpixel_o;
  size_t i, j, ind, *indexs=scp->indexs, numrows, numcols;
  float *input=scp->input, *kernel=scp->kernel, *out=scp->out;

  /* Go over all the pixels associated with this thread. */
  for(i=0;indexs[i]!=GAL_THREADS_NON_THRD_INDEX;++i)
    {
      /* Set the index, if it is a NaN pixel, then set the output to
         be NaN too. */
      ind=indexs[i];
      if(isnan(input[ind]))
        {
          out[ind]=NAN;
          continue;
        }

      /* Set the starting and ending pixels on the kernel, note that
         the overlap function in box.c uses FITS coordinates. */
      fpixel_o[0]=1;                fpixel_o[1]=1;
      lpixel_o[0]=ks1;              lpixel_o[1]=ks0;
      fpixel_i[0]=ind%is1+1-ks1/2;  fpixel_i[1]=ind/is1+1-ks0/2;
      lpixel_i[0]=ind%is1+1+ks1/2;  lpixel_i[1]=ind/is1+1+ks0/2;
      gal_box_overlap(naxes, fpixel_i, lpixel_i, fpixel_o, lpixel_o);

      /* fpixels and lpixels now point to the overlap's starting and
         ending both on the image and on the kernel. */
      sum=0.0f;
      numcols=lpixel_i[0]-fpixel_i[0]+1; /* lpixel is inside the box. */
      numrows=lpixel_i[1]-fpixel_i[1]+1; /* lpixel is inside the box. */
      ksum = edgecorrection ? 0.0f : 1.0f;
      istart =  &input[ (fpixel_i[1]-1) * is1 + fpixel_i[0]-1 ];
      kstart = &kernel[ (fpixel_o[1]-1) * ks1 + fpixel_o[0]-1 ];
      for(j=0;j<numrows;++j)
        {
          k = kstart + j*ks1;
          fp = (f = istart + j*is1) + numcols;
          do
            {
              if( isnan(*f)==0 )
                {
                  sum += *k * *f;
                  if(edgecorrection) ksum+=*k;
                }
              ++k;
            }
          while(++f<fp);
          out[ind]=sum/ksum;
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
gal_spatialconvolve_convolve(float *input, size_t is0, size_t is1,
                             float *kernel, size_t ks0, size_t ks1,
                             size_t numthreads, int edgecorrection,
                             float **out)
{
  int err;
  pthread_t t;          /* All thread ids saved in this, not used. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct gal_spatialconvolve_params *scp;
  size_t i, nb, *indexs, thrdcols;


  /* Array keeping thread parameters for each thread. */
  errno=0;
  scp=malloc(numthreads*sizeof *scp);
  if(scp==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes in gal_spatialconvolve_convolve "
          "(spatialconvolve.c) for scp", numthreads*sizeof *scp);


  /* Allocate the output array: */
  errno=0;
  *out=malloc(is0*is1*sizeof **out);
  if(*out==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for convolution output",
          is0*is1*sizeof **out);


  /* Distribute the image pixels into the threads: */
  gal_threads_dist_in_threads(is0*is1, numthreads, &indexs, &thrdcols);

  /* Start the convolution. */
  if(numthreads==1)
    {
      gal_spatialconvolve_pparams(input, is0, is1, kernel, ks0, ks1,
                                  numthreads, edgecorrection, *out, indexs,
                                  &scp[0]);
      gal_spatialconvolve_thread(&scp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread (that
         spinns off the numthreads threads) is also a thread, so the number
         the barrier should be one more than the number of threads spinned
         off. */
      if(is0*is1<numthreads) nb=is0*is1+1;
      else nb=numthreads+1;
      gal_threads_attr_barrier_init(&attr, &b, nb);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
        if(indexs[i*thrdcols]!=GAL_THREADS_NON_THRD_INDEX)
          {
            scp[i].b=&b;
            gal_spatialconvolve_pparams(input, is0, is1, kernel, ks0,
                                        ks1, numthreads, edgecorrection,
                                        *out, &indexs[i*thrdcols], &scp[i]);
            err=pthread_create(&t, &attr, gal_spatialconvolve_thread,
                               &scp[i]);
            if(err)
              error(EXIT_FAILURE, 0, "can't create thread %zu", i);
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
