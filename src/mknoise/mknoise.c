/*********************************************************************
MakeNoise - Add noise to a dataset.
MakeNoise is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <sys/time.h>		 /* Generate random seed. */
#include <gsl/gsl_rng.h>	 /* Used in setrandoms.   */
#include <gsl/gsl_randist.h>	 /* To make noise.        */

#include "timing.h"
#include "fitsarrayvv.h"

#include "main.h"









unsigned long int
random_seed()
{
  struct timeval tv;
  gettimeofday(&tv,0);
  return tv.tv_sec+tv.tv_usec;
}





void
convertsaveoutput(struct mknoiseparams *p)
{
  void *array;

  /* Convert the output to the input image format: */
  if(p->inputbitpix==DOUBLE_IMG || p->doubletype)
    {
      array=p->input;
      p->inputbitpix=DOUBLE_IMG; /* In case it wasn't and p->doubletype==1 */
    }
  else
    changetype((void **)p->input, DOUBLE_IMG, p->is0*p->is1,
               p->numblank, &array, p->inputbitpix);

  /* Save the output: */
  arraytofitsimg(p->cp.output, "NoiseAdded", p->inputbitpix, array,
                 p->is0, p->is1, p->numblank, p->wcs,
                 SPACK_STRING);

  if(array!=p->input)
    free(array);
}





void
mknoise(struct mknoiseparams *p)
{
  gsl_rng *r;
  unsigned long seed;
  const gsl_rng_type *T;
  char message[VERBMSGLENGTH_V];
  double *d, *df, background=p->background, stdadd=p->stdadd;

  /* Set the random number generator parameters. */
  gsl_rng_env_setup();
  T=gsl_rng_default;
  r=gsl_rng_alloc(T);
  if(p->envseed)
    seed=gsl_rng_default_seed;
  else
    {
      seed=random_seed();
      gsl_rng_set(r,seed);
    }

  if(p->cp.verb)
    {
      sprintf(message, "Generator type: %s", gsl_rng_name(r));
      reporttiming(NULL, message, 1);
      sprintf(message, "Generator seed: %lu", seed);
      reporttiming(NULL, message, 1);
    }

  /* Add the noise: */
  df=(d=p->input)+p->is0*p->is1;
  if(p->backgroundinmean)
    {
      do *d+=background+gsl_ran_gaussian(r,sqrt(stdadd+background+*d));
      while(++d<df);
    }
  else
    { do *d+=gsl_ran_gaussian(r,sqrt(stdadd+background+*d)); while(++d<df); }

  /* Convert and save the output in the proper format: */
  convertsaveoutput(p);

  /* Free what ever was allocated here. */
  gsl_rng_free(r);
}
