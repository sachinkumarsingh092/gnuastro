/*********************************************************************
MakeNoise - Add noise to a dataset.
MakeNoise is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>            /* Generate random seed. */

#include <gsl/gsl_rng.h>         /* Used in setrandoms.   */
#include <gnuastro/fits.h>
#include <gsl/gsl_randist.h>     /* To make noise.        */

#include <gnuastro-internal/timing.h>

#include "main.h"














void
convertsaveoutput(struct mknoiseparams *p)
{
  char keyname1[FLEN_KEYWORD];
  gal_fits_list_key_t *headers=NULL;
  char keyname2[FLEN_KEYWORD], keyname3[FLEN_KEYWORD];
  char keyname4[FLEN_KEYWORD], keyname5[FLEN_KEYWORD];


  /* Add the proper information to the header of the output: */
  gal_fits_key_write_filename("INF", p->inputname, &headers, 0);
  if( !isnan(p->background_mag) )
    {
      strcpy(keyname1, "BCKGRND");
      gal_fits_key_list_add_end(&headers, GAL_TYPE_FLOAT64, keyname1, 0,
                                &p->background_mag, 0, "Background "
                                "value (in magnitude) for noise.",
                                0, NULL);
      strcpy(keyname2, "BZRPNT");
      gal_fits_key_list_add_end(&headers, GAL_TYPE_FLOAT64, keyname2, 0,
                                &p->zeropoint, 0,
                                "Zeropoint magnitude of image.", 0, NULL);
      strcpy(keyname3, "INSTRU");
      gal_fits_key_list_add_end(&headers, GAL_TYPE_FLOAT64, keyname3, 0,
                                &p->instrumental, 0,
                                "Instrumental noise in units of flux.",
                                0, NULL);
    }
  else
    {
      strcpy(keyname1, "SIGMA");
      gal_fits_key_list_add_end(&headers, GAL_TYPE_FLOAT64, keyname1, 0,
                                &p->sigma, 0, "Total noise sigma", 0, NULL);
    }
  strcpy(keyname4, "RNGTYPE");
  gal_fits_key_list_add_end(&headers, GAL_TYPE_STRING, keyname4, 0,
                            (void *)(p->rng_name), 0,
                            "Random number generator (by GSL) type.",
                            0, NULL);
  strcpy(keyname5, "RNGSEED");
  gal_fits_key_list_add_end(&headers, GAL_TYPE_ULONG, keyname5, 0,
                            &p->rng_seed, 0,
                            "Random number generator (by GSL) seed.",
                            0, NULL);

  /* Save the output: */
  p->input=gal_data_copy_to_new_type_free(p->input, p->cp.type);
  p->input->name="NOISED";
  gal_fits_img_write(p->input, p->cp.output, headers, PROGRAM_NAME);
  p->input->name=NULL;

  /* Write the configuration keywords. */
  gal_fits_key_write_filename("input", p->inputname, &p->cp.okeys, 1);
  gal_fits_key_write_config(&p->cp.okeys, "MakeNoise configuration",
                            "MKNOISE-CONFIG", p->cp.output, "0");
}





void
mknoise(struct mknoiseparams *p)
{
  double *d, *df, background=p->background;
  double instpowtwo = p->instrumental*p->instrumental;

  /* Add the noise: */
  df=(d=p->input->array)+p->input->size;
  if( !isnan(p->sigma) )
    {
      do
        *d += gsl_ran_gaussian(p->rng, p->sigma);
      while(++d<df);
    }
  else
    {
      do
        *d += ( background
                + gsl_ran_gaussian(p->rng,
                                   sqrt( instpowtwo + background + *d )) );
      while(++d<df);
    }

  /* Convert and save the output in the proper format: */
  convertsaveoutput(p);
}
