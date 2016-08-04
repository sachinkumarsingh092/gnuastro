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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>            /* Generate random seed. */
#include <gsl/gsl_rng.h>         /* Used in setrandoms.   */
#include <gsl/gsl_randist.h>     /* To make noise.        */

#include <gnuastro/timing.h>
#include <gnuastro/fits.h>

#include "main.h"














void
convertsaveoutput(struct mknoiseparams *p)
{
  void *array;
  char keyname1[FLEN_KEYWORD];
  struct gal_fits_header_ll *headers=NULL;
  char keyname2[FLEN_KEYWORD], keyname3[FLEN_KEYWORD];
  char keyname4[FLEN_KEYWORD], keyname5[FLEN_KEYWORD];

  /* Convert the output to the input image format: */
  if(p->inputbitpix==DOUBLE_IMG || p->doubletype)
    {
      array=p->input;
      p->inputbitpix=DOUBLE_IMG; /* Not converted and p->doubletype==1 */
    }
  else
    gal_fits_change_type((void **)p->input, DOUBLE_IMG, p->is0*p->is1,
                              p->anyblank, &array, p->inputbitpix);

  /* Add the proper information to the header of the output: */
  gal_fits_file_name_in_keywords("INF", p->up.inputname, &headers);
  strcpy(keyname1, "BCKGRND");
  gal_fits_add_to_fits_header_ll_end(&headers, TDOUBLE, keyname1, 0,
                                     &p->mbackground, 0, "Background "
                                     "value (in magnitude) for noise.",
                                     0, NULL);
  strcpy(keyname2, "BZRPNT");
  gal_fits_add_to_fits_header_ll_end(&headers, TDOUBLE, keyname2, 0,
                                     &p->zeropoint, 0, "Zeropoint "
                                     "magnitude of image.", 0, NULL);
  strcpy(keyname3, "STDADD");
  gal_fits_add_to_fits_header_ll_end(&headers, TDOUBLE, keyname3, 0,
                                     &p->stdadd, 0, "Instrumental noise "
                                     "in units of flux.", 0, NULL);
  strcpy(keyname4, "RNGTYPE");
  gal_fits_add_to_fits_header_ll_end(&headers, TSTRING, keyname4, 0,
                                     &p->rng_type, 0, "Random number "
                                     "generator (by GSL) type.", 0, NULL);
  strcpy(keyname5, "RNGSEED");
  gal_fits_add_to_fits_header_ll_end(&headers, TLONG, keyname5, 0,
                                     &p->rng_seed, 0, "Random number "
                                     "generator (by GSL) seed.", 0, NULL);

  /* Save the output: */
  gal_fits_array_to_file(p->cp.output, "NoiseAdded", p->inputbitpix,
                         array, p->is0, p->is1, p->anyblank, p->wcs,
                         headers, SPACK_STRING);

  if(array!=p->input)
    free(array);
}





void
mknoise(struct mknoiseparams *p)
{
  double *d, *df, background=p->background, stdadd=p->stdadd;

  /* Add the noise: */
  df=(d=p->input)+p->is0*p->is1;
  do
    *d+=background+gsl_ran_gaussian(p->rng, sqrt(stdadd+background+*d));
  while(++d<df);

  /* Convert and save the output in the proper format: */
  convertsaveoutput(p);
}
