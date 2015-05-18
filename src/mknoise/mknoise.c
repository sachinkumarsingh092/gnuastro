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
#include <string.h>
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
  struct fitsheaderll *headers=NULL;
  char keyname1[FLEN_KEYWORD], keyname2[FLEN_KEYWORD];
  char keyname3[FLEN_KEYWORD], keyname4[FLEN_KEYWORD], keyname5[FLEN_KEYWORD];

  /* Convert the output to the input image format: */
  if(p->inputbitpix==DOUBLE_IMG || p->doubletype)
    {
      array=p->input;
      p->inputbitpix=DOUBLE_IMG; /* In case it wasn't and p->doubletype==1 */
    }
  else
    changetype((void **)p->input, DOUBLE_IMG, p->is0*p->is1,
               p->numblank, &array, p->inputbitpix);

  /* Add the proper information to the header of the output: */
  filenameinkeywords("INF", p->up.inputname, &headers);
  strcpy(keyname1, "BCKGRND");
  add_to_fitsheaderllend(&headers, TDOUBLE, keyname1, 0, &p->mbackground, 0,
                         "Background value (in magnitude) for noise.",
                         0, NULL);
  strcpy(keyname2, "BZRPNT");
  add_to_fitsheaderllend(&headers, TDOUBLE, keyname2, 0, &p->zeropoint, 0,
                         "Zeropoint magnitude of image.", 0, NULL);
  strcpy(keyname3, "STDADD");
  add_to_fitsheaderllend(&headers, TDOUBLE, keyname3, 0, &p->stdadd, 0,
                         "Instrumental noise in units of flux.",
                         0, NULL);
  strcpy(keyname4, "RNGTYPE");
  add_to_fitsheaderllend(&headers, TSTRING, keyname4, 0, &p->rng_type, 0,
                         "Random number generator (by GSL) type.", 0, NULL);
  strcpy(keyname5, "RNGSEED");
  add_to_fitsheaderllend(&headers, TLONG, keyname5, 0, &p->rng_seed, 0,
                         "Random number generator (by GSL) seed.", 0, NULL);

  /* Save the output: */
  arraytofitsimg(p->cp.output, "NoiseAdded", p->inputbitpix, array,
                 p->is0, p->is1, p->numblank, p->wcs, headers,
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
      do
        *d+=background+gsl_ran_gaussian(r,sqrt(stdadd+background+*d));
      while(++d<df);
    }
  else
    {
      do
        *d+=gsl_ran_gaussian(r,sqrt(stdadd+background+*d));
      while(++d<df);
    }

  /* Convert and save the output in the proper format: */
  p->rng_seed=seed;
  strcpy(p->rng_type, gsl_rng_name(r));
  convertsaveoutput(p);

  /* Free what ever was allocated here. */
  gsl_rng_free(r);
}
