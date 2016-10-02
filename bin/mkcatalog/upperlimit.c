/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <stdlib.h>
#include <sys/time.h>

#include <gsl/gsl_rng.h>
#include <gnuastro/fits.h>
#include <gnuastro/threads.h>
#include <gnuastro/txtarray.h>
#include <gnuastro/statistics.h>
#include <gsl/gsl_statistics_double.h>

#include <timing.h>

#include "upperlimit.h"




/*********************************************************************/
/************           Structures and macros         ****************/
/*********************************************************************/
/* These structures and macros are defined in the C file since they are
   only necessary in the processing here, they do not need to be given to
   the user through the header, it will just complicate their
   name-space. */
struct upperlimitparams
{
  size_t                  s0;   /* Number of rows in image.          */
  size_t                  s1;   /* Number of columns in image.       */
  int                envseed;   /* ==1: evironment for rng setup.    */
  long                minlab;   /* Minimum label in image.           */
  long               numlabs;   /* Number of labels in image.        */
  float                 *img;   /* The input image.                  */
  size_t                *box;	/* The box parameters of each label. */
  unsigned char        **pix;   /* The pixels of each label.         */
  float                 *std;	/* Standard deviation of the sums.   */
  size_t               upnum;   /* Number of random samples.         */
  size_t          numthreads;   /* Number of threads to use.         */
  float          sclipmultip;   /* Multiple of STD for sigma clip.   */
  float            sclipaccu;   /* Accuracy to stop sigma clipping.  */
};





/* Information for each thread. */
struct tupperlimitparams
{
  size_t                  id;   /* Id of this thread.                */
  size_t             *indexs;   /* Indexes for this thread.          */
  pthread_barrier_t       *b;   /* Barrier for all threads.          */
  struct upperlimitparams *p;   /* Pointer to main program struct.   */
};





/* Columns of labeled region information. Note that the width columns
   (WIDCOL) are actually the same as the max columns (MAXCOL). The
   width values are written over the maximum columns because we don't
   need the maximum values any more. */
#define XMINCOL 0
#define YMINCOL 1
#define XMAXCOL 2
#define YMAXCOL 3
#define XWIDCOL XMAXCOL
#define YWIDCOL YMAXCOL




















/*********************************************************************/
/************               Preparations              ****************/
/*********************************************************************/
static void
fillseginfo(struct upperlimitparams *p, long *seg)
{
  float *f, *ff;
  long *l, *ll, maxlab;
  unsigned char *u, *uu;
  size_t xw, yw, xmin, ymin;
  size_t i, j, asize, label, size=p->s0*p->s1;


  /* Find the minimum and maximum labels: */
  maxlab    = INT32_MIN;
  p->minlab = INT32_MAX;
  ll=(l=seg)+size;
  do
    /* We don't want to look at pixels with value zero or blank pixels, so
       ignore them.  */
    if( *l!=0 && *l!=GAL_FITS_LONG_BLANK )
      {
        if(*l>maxlab)    maxlab=*l;
        if(*l<p->minlab) p->minlab=*l;
      }
  while(++l<ll);


  /* Do a small sanity check: */
  if(maxlab<0 || p->minlab<0)
    error(EXIT_FAILURE, 0, "the labeled image must not contain negative "
          "pixels");


  /* Allocate the space to keep the object information. This array will
     keep the minimum and maximum coordinate for each labeled region. Note
     that max is in the image, and we want to use the labels as indexs (row
     numbers), so we must allocate max+1 rows.

     Note: we are not assuming that the first index is 1. Because many
     cases can arise where the lowest index can be a large number (for
     example when the user can be working on a sub-set of objects selected
     from a large catalog), so we don't want to allocate a huge amount of
     memory and waste thread resources on them. We have stored the minimum
     label and will use that in reporting the final output to make a full
     array. */
  p->numlabs = maxlab - p->minlab + 1;
  errno=0;
  p->box=malloc(p->numlabs * 4 * sizeof *p->box);
  if(p->box==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->box",
          p->numlabs * 4 * sizeof *p->box);

  /* Initialize the information array*/
  for(i=0;i<p->numlabs;++i)
    {
      p->box[ i*4 + XMINCOL ] = p->s0;
      p->box[ i*4 + YMINCOL ] = p->s1;
      p->box[ i*4 + XMAXCOL ] = p->box[ i*4 + YMAXCOL ] = 0;
    }

  /* Fill it in with the minimum and maximum positions. */
  for(i=0;i<p->s0;++i)
    for(j=0;j<p->s1;++j)
      if(seg[i*p->s1+j]>0)
	{
	  label = seg[i*p->s1+j] - p->minlab;
	  if(i<p->box[ label*4 + XMINCOL ]) p->box[ label*4 + XMINCOL ] = i;
	  if(j<p->box[ label*4 + YMINCOL ]) p->box[ label*4 + YMINCOL ] = j;
	  if(i>p->box[ label*4 + XMAXCOL ]) p->box[ label*4 + XMAXCOL ] = i;
	  if(j>p->box[ label*4 + YMAXCOL ]) p->box[ label*4 + YMAXCOL ] = j;
	}

  /* For a check. Note that in DS9 the coordinates are inverted and
     start from 1, not zero.
  {
    i=0;
    printf("%lu: (%lu, %lu) --> (%lu, %lu)\n", i+p->minlab,
	   p->box[i*4+YMINCOL]+1, p->box[i*4+XMINCOL]+1,
	   p->box[i*4+YMAXCOL]+1, p->box[i*4+XMAXCOL]+1);
  }
  */


  /* Allocate the pointers to the object's pixels and allocate space
     for the area of each detection. In the meantime, Change the xmax
     and ymax to width in X and width in Y*/
  errno=0;
  p->pix=malloc(p->numlabs * sizeof *p->pix);
  if(p->pix==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->pix",
          p->numlabs * sizeof *p->pix);
  for(i=0; i<p->numlabs; ++i)
    {
      /* Replace the max values with the width of each box. */
      p->box[i*4+XWIDCOL] = p->box[i*4+XMAXCOL] - p->box[i*4+XMINCOL] + 1;
      p->box[i*4+YWIDCOL] = p->box[i*4+YMAXCOL] - p->box[i*4+YMINCOL] + 1;

      /* Allocate space for the pixels of this label. */
      errno=0;
      asize = p->box[i*4+XWIDCOL] * p->box[i*4+YWIDCOL] * sizeof **p->pix;
      p->pix[i] = malloc(asize);
      if(p->pix[i]==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for p->pix[i] (i=%lu)",
              asize, i);

      /* Set the array values. Any pixel that is equal to the object's
	 label over the box is given a value of 1.0f, other pixels are
	 given a value of 0.0f. We set the starting pointers for the
	 larger and smaller arrays in this row, then go along the row
	 until it finshes. `j' will take us to the next row after this
	 row finishes. */
      xmin = p->box[ i*4+XMINCOL ];
      ymin = p->box[ i*4+YMINCOL ];
      xw   = p->box[ i*4+XWIDCOL ];
      yw   = p->box[ i*4+YWIDCOL ];
      for(j=0; j<xw; ++j)
	{
	  uu = ( u = p->pix[i] + yw*j ) + yw;
	  l  = seg + (xmin+j)*p->s1 + ymin;
	  do *u = (*l++ == i+p->minlab); while(++u<uu);
	}
    }


  /* For a check:
  i=0;
  gal_fits_array_to_file("tmpf.fits", "check", BYTE_IMG, p->pix[i],
			 p->box[i*4+XWIDCOL], p->box[i*4+YWIDCOL], 0,
			 NULL, NULL, "myprog");
  exit(0);
  */

  /* Allocate an array for the output values, note that this is for the
     full list of IDs, not just those in the image. */
  p->std = malloc( (maxlab+1) * sizeof *p->std );
  if(p->std==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->std (upperlimit.c)",
          (maxlab+1) * sizeof *p->std);
  ff=(f=p->std)+(maxlab+1); do *f++=NAN; while(f<ff);
}





/* Subtract the Sky and set all the segmentation and mask pixels to
   NaN. Note that the mask and sky images are optional. If the caller
   didn't want them, they have to be set to NULL.

   We could also be doing this during the actual flux calculation, but
   since many iterations of many objects will be done, it is more efficient
   to do all the preparations (setting blank pixels and subtracting the
   sky. */
void
prepareimg(struct upperlimitparams *p, float *img, float *sky, long *seg,
           long *mask)
{
  float *f, *ff;
  size_t size=p->s0*p->s1;

  /* Allocate the space for the image. */
  errno=0;
  p->img=malloc(size * sizeof *p->img);
  if(p->img==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->img (upperlimit.c)",
          size * sizeof *p->img);

  /* Set the values. */
  ff=(f=p->img)+size;
  do
    {
      /* Note that the Sky and mask are optional (can be NULL), so we are
         checking their existance (not NULL) first.  */
      if(mask)
        {
          *f = *img++ - (( *seg || *mask )  ? NAN :  ( sky ? *sky : 0 ));
          ++mask;
        }
      else
        *f = *img++ - (*seg  ? NAN :  ( sky ? *sky : 0 ));

      /* Increment the pointers (when they exist) to the next pixel. */
      ++seg;
      if(sky) ++sky;
    }
  while(++f<ff);

  /* For a check:
  gal_fits_array_to_file("in.fits", "check", FLOAT_IMG, p->img,
			 p->s0, p->s1, 0, NULL, NULL, "myprog");
  exit(0);
  */
}




















/*********************************************************************/
/************       Actual operation on threads       ****************/
/*********************************************************************/
static void *
upperlimit_on_thread(void *inparam)
{
  /* The first thing to do is to say what the input pointer actually is. */
  struct tupperlimitparams *tp = (struct tupperlimitparams *)inparam;
  struct upperlimitparams *p = tp->p;

  /* Now you can go onto do defining the function like any other
     function: first you define the variables and so on... */
  double *sum;
  gsl_rng *rng;
  unsigned char *u, *uu;
  float *l, ave, med, *fsum;
  size_t i, j, c, xw, yw, lab, xmin, ymin;

  /* Allcate space to keep the total sums */
  errno=0;
  sum=malloc(p->upnum*sizeof *sum);
  if(sum==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for sum (upperlimit.c)",
          p->upnum*sizeof *sum);

  /* Allocate the random number generator for this thread. Note that when
     envseed is non-zero, then we need to use the given environment
     value. However, we cannot use the same value for all threads because
     all the positions of all threads will be the same, so we add the
     thread-id (which is also a known paramter) to the seed. A given number
     of segments will be distributed between a given number of threads in a
     reproducible manner, so with this seed value for this thread, we will
     get reproducible results.*/
  rng=gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, gsl_rng_default_seed+tp->id);


  /* Go over the jobs indexed for this thread: */
  for(i=0; tp->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    {
      /* Setup this labels parameters */
      c=0;
      lab=tp->indexs[i];
      xw = p->box[ lab*4 + XWIDCOL ];
      yw = p->box[ lab*4 + YWIDCOL ];

      /* The starting point should be placed between 0 and Image_width
	 - box_width. */
      while(c<p->upnum)
	{
	  /* Get a randomly placed start for this label */
	  sum[c] = 0.0f;
	  xmin = gsl_rng_uniform_int(rng, p->s0-xw);
	  ymin = gsl_rng_uniform_int(rng, p->s1-yw);

	  /* Go over the pixels and get the sum of the pixel values,
	     this is very similar to the loop in fillseginfo. */
	  for(j=0; j<xw; ++j)
	    {
	      uu = ( u = p->pix[lab] + yw*j ) + yw;
	      l  = p->img + (xmin+j)*p->s1 + ymin;

	      /* We can't simply multiply, since NaN values might be
		 involved over regions that we don't need.*/
	      do {sum[c] += *u==1.0f ? *l : 0.0f; ++l;} while(++u<uu);
	    }

	  /* Only use this sum if it is not a NaN. */
          /*printf("%lu: %f\n", c, sum[c]);*/
	  if(!isnan(sum[c])) ++c;
	}

      /* Get the standard deviation by sigma-clipping */
      gal_fits_change_type (sum, DOUBLE_IMG, c, 0, (void **)&fsum, FLOAT_IMG);
      gal_statistics_sigma_clip_converge(fsum, 0, c, p->sclipmultip,
                                         p->sclipaccu, &ave, &med,
					 &p->std[p->minlab+lab], 0);

      /*
      printf("%f\n", p->std[ lab*p->numoutcols + p->nimg ]);
      exit(0);

      p->std[ lab*p->numoutcols + p->nimg ] = gsl_stats_sd(sum, 1, c);
      */
    }

  /* Free the random number generator space. */
  free(sum);
  gsl_rng_free(rng);

  /* Wait until all other threads finish. When there was only one thread,
     we explicitly set the pointer to the barrier structure to NULL, so
     only wait when a barrier is actually defined.*/
  if(tp->b)
    pthread_barrier_wait(tp->b);

  /* Return the NULL pointer. */
  return NULL;
}





/* This is the thread spinner function. For each image we will divide
   all the objects between the threads and the  */
void
upperlimit_manager(struct upperlimitparams *p)
{
  int err;
  pthread_t t;          /* All thread ids saved in this, not used. */
  size_t numbarriers;
  pthread_attr_t attr;
  pthread_barrier_t b;
  size_t i, *indexs, thrdcols;
  struct tupperlimitparams *tp;

  /* Allocate the array of `param' structures. Note that in most cases, the
     number of threads will not be a constant like this simple case, it
     will be a variable passed to the thread-spinner. So we are using
     dynamic allocation for more general use as a tutorial. */
  errno=0;
  tp = malloc(p->numthreads*sizeof *tp);
  if( tp==NULL )
    error(EXIT_FAILURE, errno, "%lu bytes for tp (upperlimit.c)",
          p->numthreads*sizeof *tp);

  /* Distribute the actions into the threads: */
  gal_threads_dist_in_threads(p->numlabs, p->numthreads, &indexs, &thrdcols);

  /* Do the job: when only one thread is necessary, there is no need to
     spin off one thread, just call the function directly (spinning off
     threads is expensive). This is for the generic thread spinner
     function, not this simple function where `p->numthreads' is a
     constant. */
  if(p->numthreads==1)
    {
      tp[0].p=p;
      tp[0].id=0;
      tp[0].b=NULL;
      tp[0].indexs=indexs;
      upperlimit_on_thread(&tp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
         (that spinns off the nt threads) is also a thread, so the
         number the barriers should be one more than the number of
         threads spinned off. */
      numbarriers = ( (p->numlabs<p->numthreads ? p->numlabs : p->numthreads)
                      + 1 );
      gal_threads_attr_barrier_init(&attr, &b, numbarriers);

      /* Spin off the threads: */
      for(i=0;i<p->numthreads;++i)
        if(indexs[i*thrdcols]!=GAL_THREADS_NON_THRD_INDEX)
          {
            tp[i].p=p;
            tp[i].id=i;
            tp[i].b=&b;
            tp[i].indexs=&indexs[i*thrdcols];
            err=pthread_create(&t, &attr, upperlimit_on_thread, &tp[i]);
            if(err)
              {
                fprintf(stderr, "can't create thread %lu", i);
                exit(EXIT_FAILURE);
              }
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }

  /* Clean up. */
  free(tp);
  free(indexs);
}




















/*********************************************************************/
/************             External function           ****************/
/*********************************************************************/
float *
upperlimit(float *img, float *sky, long *seg, long *mask, size_t s0,
           size_t s1, size_t upnum, size_t numthreads, int envseed,
           float sclipmultip, float sclipaccu)
{
  size_t i;
  struct upperlimitparams p;

  /* Put the input parameters in the structure. */
  p.s0=s0;
  p.s1=s1;
  p.upnum=upnum;
  p.envseed=envseed;
  p.sclipaccu=sclipaccu;
  p.numthreads=numthreads;
  p.sclipmultip=sclipmultip;

  /* Fill the information for each segment. */
  fillseginfo(&p, seg);

  /* Prepare the image. */
  prepareimg(&p, img, sky, seg, mask);

  /* Find the upper limit for all the objects on a thread. */
  upperlimit_manager(&p);

  /* For a check:
  for(i=1;i<p.minlab+p.numlabs;++i)
    printf("i=%lu: %-10.5f\n", i, p.std[i]);
  exit(0);
  */

  /* Clean up and return. */
  for(i=0; i<p.numlabs; ++i)
    free(p.pix[i]);
  free(p.img);
  free(p.pix);
  free(p.box);
  return p.std;
}
