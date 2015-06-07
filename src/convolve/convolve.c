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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>

#include "timing.h"
#include "astrthreads.h"
#include "fitsarrayvv.h"
#include "spatialconvolve.h"

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

  errno=0;
  *output=out=malloc(size*sizeof *out);
  if(out==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for out in complextoreal",
          size*sizeof *out);

  of=(o=out)+size;
  switch(action)
    {
    case COMPLEXTOREALSPEC:
      do { *o++ = sqrt( *c**c + *(c+1)**(c+1) ); c+=2; } while(o<of);
      break;
    case COMPLEXTOREALPHASE:
      do { *o++ = atan2( *(c+1), *c );           c+=2; } while(o<of);
      break;
    case COMPLEXTOREALREAL:
      do { *o++ = *c;                            c+=2; } while(o<of);
      break;
    default:
      error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can "
            "correct it. For some reason, the action code for complextoreal "
            "in convolve.c (%d) is not recognized.", PACKAGE_BUGREPORT,
            action);
    }
}





/* Multily two complex arrays and save the result in the first:

   (a+ib)*(c+id)=ac+iad+ibc-bd=(ac-bd)+i(ad-bc)
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
      *a++=r;     /* Go onto the imaginary part of a. */
      b+=2;
    }
  while(++a<af);  /* Go onto the next complex number. */
}




















/******************************************************************/
/*************      Padding and initializing      *****************/
/******************************************************************/
void
makepaddedcomplex(struct convolveparams *p)
{
  size_t i, ps0, ps1;
  double *o, *op, *pimg, *pker;
  float *f, *fp, *input=p->input, *kernel=p->kernel;
  size_t is0=p->is0, is1=p->is1, ks0=p->ks0, ks1=p->ks1;

  /* Find the sizes of the padded image, note that since the kernel
     sizes are always odd, the extra padding on the input image is
     always going to be an even number (clearly divisable). */
  ps0=p->ps0=p->is0+p->ks0-1;   /* Find the padded image size. */
  ps1=p->ps1=p->is1+p->ks1-1;

  errno=0;                      /* Pad the input image.        */
  pimg=p->pimg=malloc(2*ps0*ps1*sizeof *pimg);
  if(pimg==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for pimg.", ps0*ps1*sizeof *pimg);
  for(i=0;i<ps0;++i)
    {
      op=(o=pimg+i*2*ps1)+2*ps1; /* pimg is complex.            */
      if(i<is0)
        {
          fp=(f=input+i*is1)+is1;
          do {*o++=*f; *o++=0.0f;} while(++f<fp);
        }
      do *o++=0.0f; while(o<op);
    }


  errno=0;                      /* Pad the Kernel.             */
  pker=p->pker=malloc(2*ps0*ps1*sizeof *pker);
  if(pker==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for pimg.", ps0*ps1*sizeof *pker);
  for(i=0;i<ps0;++i)
    {
      op=(o=pker+i*2*ps1)+2*ps1; /* pker is complex.            */
      if(i<ks0)
        {
          fp=(f=kernel+i*ks1)+ks1;
          do {*o++=*f; *o++=0.0f;} while(++f<fp);
        }
      do *o++=0.0f; while(o<op);
    }
}





/*  Remove the padding from the final convolved image and also correct
    for roundoff errors. Notice that we are putting the pixels in the
    input image of float type.

    NOTE: The padding to the input image (on the first axis for
          example) was p->ks0-1. Since p->ks0 is always odd, the
          padding will alwys be even.  */
void
removepaddingcorrectroundoff(struct convolveparams *p)
{
  size_t ps1=p->ps1;
  float *o, *input=p->input;
  double *d, *df, *start, *pimg=p->pimg;
  size_t i, hi0, hi1, is0=p->is0, is1=p->is1;

  hi0=(p->ks0-1)/2;
  hi1=(p->ks1-1)/2;

  /* To start with, `start' points to the first pixel in the final
     image: */
  start=&pimg[hi0*ps1+hi1];
  for(i=0;i<is0;++i)
    {
      o=&input[i*is1];
      df=(d=start+i*ps1)+is1;
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
    error(EXIT_FAILURE, errno, "%lu bytes in forward2dfft for fp",
          p->cp.numthreads*sizeof *fp);

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
  for(i=0;indexs[i]!=NONTHRDINDEX;++i)
    {
      data = ( indexs[i]<maxindex
               ? &pimg[ 2*indexs[i]*indmultip ]   /* *2 because complex. */
               : &pker[ 2*(indexs[i]-maxindex)*indmultip ]);

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
    error(EXIT_FAILURE, 0, "A bug! In forward2dfft, the value of the "
          "variable forward1backwardn1 is somehow not 1 or 2, but %d. "
          "Please contact us at "PACKAGE_BUGREPORT" so we can find the "
          "cause of the problem and fix it.", forward1backwardn1);


  /* ==================== */
  /* 1D FFT on each row. */
  /* ==================== */
  distinthreads(multiple*p->ps0, nt, &indexs, &thrdcols);
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
      attrbarrierinit(&attr, &b, nb);

      /* Spin off the threads: */
      for(i=0;i<nt;++i)
        if(indexs[i*thrdcols]!=NONTHRDINDEX)
          {
            fp[i].id=i;
            fp[i].b=&b;
            fp[i].stride=1; /* On each row, stride=1 */
            fp[i].indexs=&indexs[i*thrdcols];
            fp[i].forward1backwardn1=forward1backwardn1;
	    err=pthread_create(&t, &attr, onedimensionfft, &fp[i]);
	    if(err)
	      error(EXIT_FAILURE, 0, "Can't create thread %lu for rows.",
                    i);
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
  distinthreads(multiple*p->ps1, nt, &indexs, &thrdcols);
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
      attrbarrierinit(&attr, &b, nb);
      for(i=0;i<nt;++i)
        if(indexs[i*thrdcols]!=NONTHRDINDEX)
          {
            fp[i].b=&b;
            fp[i].stride=p->ps1; /* On each column, stride is p->ps1 */
            fp[i].indexs=&indexs[i*thrdcols];
            fp[i].forward1backwardn1=forward1backwardn1;
	    err=pthread_create(&t, &attr, onedimensionfft, &fp[i]);
	    if(err)
	      error(EXIT_FAILURE, 0, "Can't create thread %lu for columns.",
                    i);
          }
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }
  free(indexs);
}





void
frequencyconvolve(struct convolveparams *p)
{
  double *tmp;
  struct timeval t1;
  int verb=p->cp.verb;
  struct fftonthreadparams *fp;

  /* Make the padded arrays. */
  if(verb) gettimeofday(&t1, NULL);
  makepaddedcomplex(p);
  if(verb) reporttiming(&t1, "Input and Kernel images padded.", 1);
  if(p->viewfreqsteps)
    {
      complextoreal(p->pimg, p->ps0*p->ps1, COMPLEXTOREALREAL, &tmp);
      arraytofitsimg(p->up.freqstepsname, "Input transform", DOUBLE_IMG,
                     tmp, p->ps0, p->ps1, 0, NULL, NULL, SPACK_STRING);
      free(tmp);
      complextoreal(p->pker, p->ps0*p->ps1, COMPLEXTOREALREAL, &tmp);
      arraytofitsimg(p->up.freqstepsname, "Kernel transform", DOUBLE_IMG,
                     tmp, p->ps0, p->ps1, 0, NULL, NULL, SPACK_STRING);
      free(tmp);
    }


  /* Initialize the structures: */
  fftinitializer(p, &fp);


  /* Forward 2D FFT on each image. */
  if(verb) gettimeofday(&t1, NULL);
  twodimensionfft(p, fp, 1);
  if(verb) reporttiming(&t1, "Images converted to frequency domain.", 1);
  if(p->viewfreqsteps)
    {
      complextoreal(p->pimg, p->ps0*p->ps1, COMPLEXTOREALSPEC, &tmp);
      arraytofitsimg(p->up.freqstepsname, "Input transform", DOUBLE_IMG,
                     tmp, p->ps0, p->ps1, 0, NULL, NULL, SPACK_STRING);
      free(tmp);
      complextoreal(p->pker, p->ps0*p->ps1, COMPLEXTOREALSPEC, &tmp);
      arraytofitsimg(p->up.freqstepsname, "Kernel transform", DOUBLE_IMG,
                     tmp, p->ps0, p->ps1, 0, NULL, NULL, SPACK_STRING);
      free(tmp);
    }


  /* Multiply the two arrays and save them in the output.*/
  if(verb) gettimeofday(&t1, NULL);
  complexarraymultiply(p->pimg, p->pker, p->ps0*p->ps1);
  if(verb) reporttiming(&t1, "Multiplied in the frequency domain.", 1);
  if(p->viewfreqsteps)
    {
      complextoreal(p->pimg, p->ps0*p->ps1, COMPLEXTOREALSPEC, &tmp);
      arraytofitsimg(p->up.freqstepsname, "Multiplied", DOUBLE_IMG,
                     tmp, p->ps0, p->ps1, 0, NULL, NULL, SPACK_STRING);
      free(tmp);
    }



  /* Forward (in practice inverse) 2D FFT on each image. */
  if(verb) gettimeofday(&t1, NULL);
  twodimensionfft(p, fp, -1);
  complextoreal(p->pimg, p->ps0*p->ps1, COMPLEXTOREALREAL, &tmp);
  if(verb) reporttiming(&t1, "Converted back to the spatial domain.", 1);
  if(p->viewfreqsteps)
    {
      arraytofitsimg(p->up.freqstepsname, "Spatial", DOUBLE_IMG,
                     tmp, p->ps0, p->ps1, 0, NULL, NULL, SPACK_STRING);
      free(tmp);
    }


  /* Free the padded arrays (they are no longer needed) and put the
     converted array (that is real, not complex) in p->pimg. */
  free(p->pimg);
  free(p->pker);
  p->pimg=tmp;


  /* Crop out the center, numbers smaller than 10^{-17} are errors,
     remove them. */
  if(verb) gettimeofday(&t1, NULL);
  removepaddingcorrectroundoff(p);
  if(verb) reporttiming(&t1, "Padded parts removed.", 1);




  /* Free all the allocated space. */
  freefp(fp);
}




















/******************************************************************/
/*************          Outside function          *****************/
/******************************************************************/
void
convolve(struct convolveparams *p)
{
  float *convolved;
  long *meshindexs;
  struct meshparams *mp=&p->mp;

  /* Do the convolution. */
  if(p->spatial)
    {
      /* Prepare the mesh structure: */
      mp->img=p->input;      mp->s0=p->is0;      mp->s1=p->is1;
      mp->kernel=p->kernel;  mp->ks0=p->ks0;     mp->ks1=p->ks1;
      mp->numthreads=p->cp.numthreads;
      makemesh(mp);
      if(p->meshname)
        {
          checkmeshid(mp, &meshindexs);
          arraytofitsimg(p->meshname, "Input", FLOAT_IMG, p->mp.img,
                         mp->s0, mp->s1, p->numblank, p->wcs, NULL,
                         SPACK_STRING);
          arraytofitsimg(p->meshname, "MeshIndexs", LONG_IMG, meshindexs,
                         mp->s0, mp->s1, 0, p->wcs, NULL, SPACK_STRING);
          free(meshindexs);
        }

      /* Do the spatial convolution on the mesh: */
      spatialconvolveonmesh(mp, &convolved);

      /* Replace the input image array with the convolved array: */
      free(p->input);
      p->input=convolved;
    }
  else
    frequencyconvolve(p);

  /* Save the output (which is in p->input) array. Note that p->input
     will be freed in ui.c. */
  arraytofitsimg(p->cp.output, "Convolved", FLOAT_IMG, p->input,
                 p->is0, p->is1, p->numblank, p->wcs, NULL, SPACK_STRING);
}
