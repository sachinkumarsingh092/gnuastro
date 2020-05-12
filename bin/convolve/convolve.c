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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

#include <gnuastro/wcs.h>
#include <gnuastro/tile.h>
#include <gnuastro/fits.h>
#include <gnuastro/pointer.h>
#include <gnuastro/threads.h>
#include <gnuastro/convolve.h>

#include <gnuastro-internal/timing.h>

#include "main.h"
#include "convolve.h"



/******************************************************************/
/*************           Complex numbers          *****************/
/******************************************************************/

/* We have a complex (R+iI) array and we want to display it. But we
   can only do that either with the spectrum, or the phase:

   Spectrum: sqrt(R^2+I^2)
   Phase:    arctan(I/R)
*/
void
complextoreal(double *c, size_t size, int action, double **output)
{
  double *out, *o, *of;

  /* Allocate the space for the real array. */
  *output=out=gal_pointer_allocate(GAL_TYPE_FLOAT64, size, 0, __func__,
                                   "output");

  /* Fill the real array with the derived value from the complex array. */
  of=(o=out)+size;
  switch(action)
    {
    case COMPLEX_TO_REAL_SPEC:
      do { *o++ = sqrt( *c**c + *(c+1)**(c+1) ); c+=2; } while(o<of);
      break;
    case COMPLEX_TO_REAL_PHASE:
      do { *o++ = atan2( *(c+1), *c );           c+=2; } while(o<of);
      break;
    case COMPLEX_TO_REAL_REAL:
      do { *o++ = *c;                            c+=2; } while(o<of);
      break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so we can "
            "correct it. The 'action' code %d is not recognized", __func__,
            PACKAGE_BUGREPORT, action);
    }
}





/* Multily two complex arrays and save the result in the first:

   (a+ib)*(c+id)=ac+iad+ibc-bd=(ac-bd)+i(ad-bc)

   The loop is easy to understand: we want to replace two
   variables. But changing one, will affect the other. So what we do,
   is to store the final value of one, then replace the second, then
   finally replace the first one.

   Here, we first get the real component but don't put it in the
   output. Then we find and replace the imaginary component, finally,
   we put the new real component in the image.
 */
void
complexarraymultiply(double *a, double *b, size_t size)
{
  double r, *af;

  af=a+2*size;
  do
    {
      r      = (*a * *b) - (*(a+1) * *(b+1));
      *(a+1) = (*(a+1) * *b) + (*a * *(b+1));
      *a++=r;            /* Go onto (set) the imaginary part of a. */
      b+=2;
    }
  while(++a<af);  /* Go onto the next complex number. */
}





/* Divide the elements of the first array by the elements of the second
   array and put the result into the elements of the first array.

   (a+ib)/(c+id)=[(a+ib)*(c-id)]/[(c+id)*(c-id)]
                =(ac-iad+ibc+bd)/(c^2+d^2)
                =[(ac+bd)+i(bc-ad)]/(c^2+d^2)

   See the explanations above complexarraymultiply for an explanation
   on the loop.
 */
void
complexarraydivide(double *a, double *b, size_t size, double minsharpspec)
{
  double r, *af;

  af=a+2*size;
  do
    {
      if (sqrt(*b**b + *(b+1)**(b+1))>minsharpspec)
        {
          r      = ( ( (*a * *b) + (*(a+1) * *(b+1)) )
                     / ( *b * *b + *(b+1) * *(b+1) ) );
          *(a+1) = ( ( (*(a+1) * *b) - (*a * *(b+1)) )
                     / ( *b * *b + *(b+1) * *(b+1) ) );
          *a=r;

          /* Just as a sanity check (the result should never be larger than
             one. */
          if(sqrt(*a**a + *(a+1)**(a+1))>1.00001f)
            *a=*(a+1)=0.0f;
        }
      else
        {
          *a=0;
          *(a+1)=0;
        }

      a+=2;
      b+=2;
    }
  while(a<af);  /* Go onto the next complex number. */
}




















/******************************************************************/
/*************      Padding and initializing      *****************/
/******************************************************************/
void
frequency_make_padded_complex(struct convolveparams *p)
{
  size_t i, ps0, ps1;
  double *o, *op, *pimg, *pker;
  size_t is0=p->input->dsize[0],  is1=p->input->dsize[1];
  size_t ks0=p->kernel->dsize[0], ks1=p->kernel->dsize[1];
  float *f, *ff, *input=p->input->array, *kernel=p->kernel->array;


  /* Find the sizes of the padded image, note that since the kernel
     sizes are always odd, the extra padding on the input image is
     always going to be an even number (clearly divisable). */
  ps0=p->ps0 = p->makekernel ? is0 : is0 + ks0 - 1;
  ps1=p->ps1 = p->makekernel ? is1 : is1 + ks1 - 1;


  /* The Discrete Fourier transforms operate faster on even-sized
     arrays. So if the padded sides are not even, make them so: */
  if(ps0%2) ps0=p->ps0=ps0+1;
  if(ps1%2) ps1=p->ps1=ps1+1;


  /* Allocate the space for the padded input image and fill it. */
  pimg=p->pimg=gal_pointer_allocate(GAL_TYPE_FLOAT64, 2*ps0*ps1, 0,
                                    __func__, "pimg");
  for(i=0;i<ps0;++i)
    {
      op=(o=pimg+i*2*ps1)+2*ps1; /* pimg is complex.            */
      if(i<is0)
        {
          ff=(f=input+i*is1)+is1;
          do {*o++=*f; *o++=0.0f;} while(++f<ff);
        }
      do *o++=0.0f; while(o<op);
    }


  /* Allocate the space for the padded Kernel and fill it. */
  pker=p->pker=gal_pointer_allocate(GAL_TYPE_FLOAT64, 2*ps0*ps1, 0,
                                    __func__, "pker");
  for(i=0;i<ps0;++i)
    {
      op=(o=pker+i*2*ps1)+2*ps1; /* pker is complex.            */
      if(i<ks0)
        {
          ff=(f=kernel+i*ks1)+ks1;
          do {*o++=*f; *o++=0.0f;} while(++f<ff);
        }
      do *o++=0.0f; while(o<op);
    }
}





/*  Remove the padding from the final convolved image and also correct for
    roundoff errors.

    NOTE: The padding to the input image (on the first axis for example)
          was 'p->kernel->dsize[0]-1'. Since 'p->kernel->dsize[0]' is
          always odd, the padding will always be even.  */
void
removepaddingcorrectroundoff(struct convolveparams *p)
{
  size_t ps1=p->ps1;
  size_t *isize=p->input->dsize;
  float *o, *input=p->input->array;
  double *d, *df, *start, *rpad=p->rpad;
  size_t i, hi0, hi1, mkwidth=2*p->makekernel-1;

  /* Set all the necessary parameters to crop the desired region. hi0 and
     hi1 are the coordinates of the first pixel in the output image. In the
     case of deconvolution, if the maximum radius is larger than the input
     image, we will also only be using region that contains non-zero rows
     and columns.*/
  if(p->makekernel)
    {
      hi0      = mkwidth < isize[0] ? p->ps0/2-p->makekernel : 0;
      hi1      = mkwidth < isize[1] ? p->ps1/2-p->makekernel : 0;
      isize[0] = mkwidth < isize[0] ? 2*p->makekernel-1 : isize[0];
      isize[1] = mkwidth < isize[1] ? 2*p->makekernel-1 : isize[1];
    }
  else
    {
      hi0 = ( p->kernel->dsize[0] - 1 )/2;
      hi1 = ( p->kernel->dsize[1] - 1 )/2;
    }

  /* To start with, 'start' points to the first pixel in the final
     image: */
  start=&rpad[hi0*ps1+hi1];
  for(i=0;i<isize[0];++i)
    {
      o = &input[ i * isize[1] ];

      df = ( d = start + i * ps1 ) + isize[1];
      do
        *o++ = ( *d<-CONVFLOATINGPOINTERR || *d>CONVFLOATINGPOINTERR )
          ? *d
          : 0.0f;
      while (++d<df);
    }
}





/* Allocate the necessary arrays, note that we put everything in the
   first element of the fftonthreadparams structure array. All the
   other elements will point to this one later. This structure will be
   given to threads to run two times with a fixed set of parameters,
   that is why we are doing this here to facilitate the job. */
void
fftinitializer(struct convolveparams *p, struct fftonthreadparams **outfp)
{
  size_t i;
  struct fftonthreadparams *fp;

  /* Allocate the fftonthreadparams array.  */
  errno=0;
  *outfp=fp=malloc(p->cp.numthreads*sizeof *fp);
  if(fp==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for fp",
          __func__, p->cp.numthreads*sizeof *fp);

  /* Initialize the gsl_fft_wavetable structures (these are thread
     safe): */
  fp[0].ps0wave=gsl_fft_complex_wavetable_alloc(p->ps0);
  fp[0].ps1wave=gsl_fft_complex_wavetable_alloc(p->ps1);

  /* Set the values for all the other threads: */
  for(i=0;i<p->cp.numthreads;++i)
    {
      fp[i].p=p;
      fp[i].ps0wave=fp[0].ps0wave;
      fp[i].ps1wave=fp[0].ps1wave;
      fp[i].ps0work=gsl_fft_complex_workspace_alloc(p->ps0);
      fp[i].ps1work=gsl_fft_complex_workspace_alloc(p->ps1);
    }
}





void
freefp(struct fftonthreadparams *fp)
{
  size_t i;
  gsl_fft_complex_wavetable_free(fp[0].ps0wave);
  gsl_fft_complex_wavetable_free(fp[0].ps1wave);
  for(i=0;i<fp->p->cp.numthreads;++i)
    {
      gsl_fft_complex_workspace_free(fp[i].ps0work);
      gsl_fft_complex_workspace_free(fp[i].ps1work);
    }
  free(fp);
}





/* Unfortunately I don't understand why the division operation in
   deconvolution (makekernel) does not produce a centered image, the
   image is translated by half the input size in both dimensions. So I
   am correcting this in the spatial domain here. */
void
correctdeconvolve(struct convolveparams *p, double **spatial)
{
  double r, *s, *n, *d, *df, sum=0.0f;
  size_t i, j, ps0=p->ps0, ps1=p->ps1;
  int ii, jj, ci=p->ps0/2-1, cj=p->ps1/2-1;

  /* Check if the image has even sides. */
  if(ps0%2 || ps1%2)
    error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s. The padded "
          "image sides are not an even number", __func__, PACKAGE_BUGREPORT);

  /* First convert the complex image to a real image: */
  complextoreal(p->pimg, ps0*ps1, COMPLEX_TO_REAL_SPEC, &s);

  /* Allocate the array to keep the new values */
  errno=0;
  n=malloc(ps0*ps1*sizeof *n);
  if(n==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'n'",
          __func__, ps0*ps1*sizeof *n);


  /* Put the elements in their proper place: For example in one
     dimension where the values are actually the true distances:

        s[0]=0, s[1]=1, s[2]=2, s[3]=3, s[4]=4, s[5]=5

     We want the value 0 to be in the 'center'. Note that 's' is
     periodic, for example the next 6 elements have distances:

        s[6]=0, s[7]=1, s[8]=2, s[9]=3, s[10]=4, s[11]=5

     So a 'center'ed array would be like:

        s[0]=4, s[1]=5, s[2]=0, s[3]=1, s[4]=2, s[5]=3

     The relations between the old (i and j) and new (ii and jj) come
     from something like the above line.
   */
  for(i=0;i<ps0;++i)
    {
      ii= i>ps0/2 ? i-(ps0/2+1) : i+ps0/2-1;
      for(j=0;j<ps1;++j)
        {
          jj = j>ps1/2 ? j-(ps1/2+1) : j+ps1/2-1;

          r=sqrt( (ii-ci)*(ii-ci) + (jj-cj)*(jj-cj) );
          sum += n[ii*ps1+jj] = r < p->makekernel ? s[i*ps1+j] : 0;

          /*printf("(%zu, %zu) --> (%zu, %zu)\n", i, j, ii, jj);*/
        }
    }


  /* Divide all elements by the sum so the kernel is normalized: */
  df=(d=n)+ps0*ps1; do *d++/=sum; while(d<df);


  /* Clean up: */
  free(s);
  *spatial=n;
}


















/******************************************************************/
/*************    Frequency domain convolution    *****************/
/******************************************************************/
/* The indexs array specifies the row or column numbers for this
  thread to work on. If forward1backwardn1 is one, then this is the
  forward transform, meaning that in convolution there are two
  images. If it is -1, then this is the final backward transform and
  there is only one image to run FFTW on and the values in indexs will
  always be smaller than p->s0 and p->s1. When there are two images,
  then the index numbers are going to be at most double p->s0 and
  p->s1. In this case, those index values which are smaller than p->s0
  or p->s1 belong to the input image and those which are equal or
  larger than larger belong to the kernel image (after subtraction for
  p->s0 or p->s1).*/
void *
onedimensionfft(void *inparam)
{
  struct fftonthreadparams *fp = (struct fftonthreadparams *)inparam;
  struct convolveparams *p=fp->p;

  double *d, *df;
  size_t indmultip, maxindex;
  gsl_fft_complex_workspace *work;
  gsl_fft_complex_wavetable *wavetable;
  double *data, *pimg=p->pimg, *pker=p->pker;
  int forward1backwardn1=fp->forward1backwardn1;
  size_t i, size, stride=fp->stride, *indexs=fp->indexs;

  /* Set the number of points to transform,

     indmultip: The value to be multiplied by the value in indexs to
     specify the first pixel of the row or column.
   */
  if(stride==1)
    { size=p->ps1; wavetable=fp->ps1wave; work=fp->ps1work;
      maxindex=p->ps0; indmultip=p->ps1; }
  else
    { size=p->ps0; wavetable=fp->ps0wave; work=fp->ps0work;
      maxindex=p->ps1; indmultip=1;      }


  /* Go over all the rows or columns given for this thread.

     NOTE: The final array (after the two FFT'd arrays are multiplied
     by each other) is stored in p->pimg. So the check below works
     both in the forward and the backward transformation.
  */
  for(i=0; indexs[i]!=GAL_BLANK_SIZE_T; ++i)
    {
      data = ( indexs[i]<maxindex
               ? &pimg[ 2*indexs[i]*indmultip ]   /* *2 because complex. */
               : &pker[ 2*(indexs[i]-maxindex)*indmultip ] );

      gsl_fft_complex_transform(data, stride, size, wavetable, work,
                                forward1backwardn1);

      /* Normalize in the backward transform: */
      if(forward1backwardn1==-1)
        {
          df=(d=data)+2*size*stride;
          do {*d/=size; *(d+1)/=size; d+=2*stride;} while(d<df);
        }
    }

  /* Wait until all other threads finish. */
  if(p->cp.numthreads>1)
    pthread_barrier_wait(fp->b);
  return NULL;
}





/* Do the forward Fast Fourier Transform either on two input images
   (the padded image and kernel) or on one image (the multiplication
   of the FFT of the two). In the second case, it is assumed that we
   are looking at the complex conjugate of the array so in practice
   this will be a backward transform. */
void
twodimensionfft(struct convolveparams *p, struct fftonthreadparams *fp,
                int forward1backwardn1)
{
  int err;
  pthread_t t;          /* All thread ids saved in this, not used. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  size_t i, nb, *indexs, thrdcols;
  size_t nt=p->cp.numthreads, multiple=0;

  /* First we are going to get the 1D fourier transform on the rows of
     both images. */
  if(forward1backwardn1==1)       multiple=2;
  else if(forward1backwardn1==-1) multiple=1;
  else
    error(EXIT_FAILURE, 0, "%s: a bug! The value of the variable "
          "'forward1backwardn1' is %d not 1 or 2. Please contact us at %s "
          "so we can find the cause of the problem and fix it", __func__,
          forward1backwardn1, PACKAGE_BUGREPORT);


  /* ==================== */
  /* 1D FFT on each row. */
  /* ==================== */
  gal_threads_dist_in_threads(multiple*p->ps0, nt, &indexs, &thrdcols);
  if(nt==1)
    {
      fp[0].stride=1;
      fp[0].indexs=&indexs[0];
      fp[0].forward1backwardn1=forward1backwardn1;
      onedimensionfft(&fp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
         (that spinns off the nt threads) is also a thread, so the
         number the barrier should be one more than the number of
         threads spinned off. */
      if( multiple*p->ps0 < nt ) nb=multiple*p->ps0+1;
      else nb=nt+1;
      gal_threads_attr_barrier_init(&attr, &b, nb);

      /* Spin off the threads: */
      for(i=0;i<nt;++i)
        if(indexs[i*thrdcols]!=GAL_BLANK_SIZE_T)
          {
            fp[i].id=i;
            fp[i].b=&b;
            fp[i].stride=1; /* On each row, stride=1 */
            fp[i].indexs=&indexs[i*thrdcols];
            fp[i].forward1backwardn1=forward1backwardn1;
            err=pthread_create(&t, &attr, onedimensionfft, &fp[i]);
            if(err)
              error(EXIT_FAILURE, 0, "%s: can't create thread %zu for rows",
                    __func__, i);
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }
  free(indexs);



  /* ====================== */
  /* 1D FFT on each column. */
  /* ====================== */
  /* No comments, exact duplicate, except the p->ps1s! */
  gal_threads_dist_in_threads(multiple*p->ps1, nt, &indexs, &thrdcols);
  if(nt==1)
    {
      fp[0].stride=p->ps1;
      fp[0].indexs=indexs;
      fp[0].forward1backwardn1=forward1backwardn1;
      onedimensionfft(&fp[0]);
    }
  else
    {
      if( multiple*p->ps1 < nt ) nb=multiple*p->ps1+1;
      else nb=nt+1;
      gal_threads_attr_barrier_init(&attr, &b, nb);
      for(i=0;i<nt;++i)
        if(indexs[i*thrdcols]!=GAL_BLANK_SIZE_T)
          {
            fp[i].b=&b;
            fp[i].stride=p->ps1; /* On each column, stride is p->ps1 */
            fp[i].indexs=&indexs[i*thrdcols];
            fp[i].forward1backwardn1=forward1backwardn1;
            err=pthread_create(&t, &attr, onedimensionfft, &fp[i]);
            if(err)
              error(EXIT_FAILURE, 0, "%s: can't create thread %zu for columns",
                    __func__, i);
          }
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }
  free(indexs);
}





void
convolve_frequency(struct convolveparams *p)
{
  double *tmp;
  size_t dsize[2];
  struct timeval t1;
  gal_data_t *data=NULL;
  struct fftonthreadparams *fp;


  /* Make the padded arrays. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  frequency_make_padded_complex(p);
  if(!p->cp.quiet)
    gal_timing_report(&t1, "Input and Kernel images padded.", 1);
  if(p->checkfreqsteps)
    {
      /* Prepare the data structure for viewing the steps, note that we
         don't need the array that is initially made. */
      dsize[0]=p->ps0; dsize[1]=p->ps1;
      data=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 2, dsize, NULL, 0,
                          p->cp.minmapsize, p->cp.quietmmap,
                          NULL, NULL, NULL);
      free(data->array);

      /* Save the padded input image. */
      complextoreal(p->pimg, p->ps0*p->ps1, COMPLEX_TO_REAL_REAL, &tmp);
      data->array=tmp; data->name="input padded";
      gal_fits_img_write(data, p->freqstepsname, NULL, PROGRAM_NAME);
      free(tmp); data->name=NULL;

      /* Save the padded kernel image. */
      complextoreal(p->pker, p->ps0*p->ps1, COMPLEX_TO_REAL_REAL, &tmp);
      data->array=tmp; data->name="kernel padded";
      gal_fits_img_write(data, p->freqstepsname, NULL, PROGRAM_NAME);
      free(tmp); data->name=NULL;
    }


  /* Initialize the structures: */
  fftinitializer(p, &fp);


  /* Forward 2D FFT on each image. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  twodimensionfft(p, fp, 1);
  if(!p->cp.quiet)
    gal_timing_report(&t1, "Images converted to frequency domain.", 1);
  if(p->checkfreqsteps)
    {
      complextoreal(p->pimg, p->ps0*p->ps1, COMPLEX_TO_REAL_SPEC, &tmp);
      data->array=tmp; data->name="input transformed";
      gal_fits_img_write(data, p->freqstepsname, NULL, PROGRAM_NAME);
      free(tmp); data->name=NULL;

      complextoreal(p->pker, p->ps0*p->ps1, COMPLEX_TO_REAL_SPEC, &tmp);
      data->array=tmp; data->name="kernel transformed";
      gal_fits_img_write(data, p->freqstepsname, NULL, PROGRAM_NAME);
      free(tmp); data->name=NULL;
    }

  /* Multiply or divide the two arrays and save them in the output.*/
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  if(p->makekernel)
    {
      complexarraydivide(p->pimg, p->pker, p->ps0*p->ps1, p->minsharpspec);
      if(!p->cp.quiet)
        gal_timing_report(&t1, "Divided in the frequency domain.", 1);
    }
  else
    {
      complexarraymultiply(p->pimg, p->pker, p->ps0*p->ps1);
      if(!p->cp.quiet)
        gal_timing_report(&t1, "Multiplied in the frequency domain.", 1);
    }
  if(p->checkfreqsteps)
    {
      complextoreal(p->pimg, p->ps0*p->ps1, COMPLEX_TO_REAL_SPEC, &tmp);
      data->array=tmp; data->name=p->makekernel ? "Divided" : "Multiplied";
      gal_fits_img_write(data, p->freqstepsname, NULL, PROGRAM_NAME);
      free(tmp); data->name=NULL;
    }

  /* Forward (in practice inverse) 2D FFT on each image. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  twodimensionfft(p, fp, -1);
  if(p->makekernel)
    correctdeconvolve(p, &p->rpad);
  else
    complextoreal(p->pimg, p->ps0*p->ps1, COMPLEX_TO_REAL_REAL, &p->rpad);
  if(!p->cp.quiet)
    gal_timing_report(&t1, "Converted back to the spatial domain.", 1);
  if(p->checkfreqsteps)
    {
      data->array=p->rpad; data->name="padded output";
      gal_fits_img_write(data, p->freqstepsname, NULL, PROGRAM_NAME);
      data->name=NULL; data->array=NULL;
    }

  /* Free the padded arrays (they are no longer needed) and put the
     converted array (that is real, not complex) in p->pimg. */
  gal_data_free(data);
  free(p->pimg);
  free(p->pker);

  /* Crop out the center, numbers smaller than 10^{-17} are errors,
     remove them. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);
  removepaddingcorrectroundoff(p);
  if(!p->cp.quiet) gal_timing_report(&t1, "Padded parts removed.", 1);


  /* Free all the allocated space. */
  freefp(fp);
}



















/******************************************************************/
/*************          Outside function          *****************/
/******************************************************************/
void
convolve(struct convolveparams *p)
{
  gal_data_t *out, *check;
  int multidim=p->input->ndim>1;
  struct gal_options_common_params *cp=&p->cp;


  /* Do the convolution. */
  if(p->domain==CONVOLVE_DOMAIN_SPATIAL)
    {
      /* Prepare the mesh structure. */
      if(multidim) gal_tile_full_two_layers(p->input, &cp->tl);

      /* Save the tile IDs if they are requested. */
      if(multidim && cp->tl.tilecheckname)
        {
          check=gal_tile_block_check_tiles(cp->tl.tiles);
          gal_fits_img_write(check, cp->tl.tilecheckname, NULL, PROGRAM_NAME);
          gal_data_free(check);
        }

      /* Do the spatial convolution. One of the main reason someone would
         want to do spatial domain convolution with this Convolve program
         is edge correction. So by default we assume it and will only
         ignore it if the user asks.*/
      out=gal_convolve_spatial(multidim ? cp->tl.tiles : p->input, p->kernel,
                               cp->numthreads,
                               multidim ? !p->noedgecorrection : 1,
                               multidim ? cp->tl.workoverch : 1 );

      /* Clean up: free the actual input and replace it's pointer with the
         convolved dataset to save as output. */
      gal_tile_full_free_contents(&cp->tl);
      gal_data_free(p->input);
      p->input=out;
    }
  else
    convolve_frequency(p);

  /* Save the output (which is in p->input) array. */
  if(p->input->ndim==1)
    gal_table_write(p->input, NULL, p->cp.tableformat, p->cp.output,
                    "CONVOLVED", 0);
  else
    gal_fits_img_write_to_type(p->input, cp->output, NULL, PROGRAM_NAME,
                               cp->type);

  /* Write Convolve's parameters as keywords into the first extension of
     the output. */
  if( gal_fits_name_is_fits(p->cp.output) )
    {
      gal_fits_key_write_filename("input", p->filename, &cp->okeys, 1);
      gal_fits_key_write_config(&cp->okeys, "Convolve configuration",
                                "CONVOLVE-CONFIG", cp->output, "0");
    }
  printf("  - Output: %s\n", p->cp.output);
}
