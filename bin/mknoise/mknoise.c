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
#include <gnuastro-internal/checkset.h>

#include "main.h"














void
convertsaveoutput(struct mknoiseparams *p)
{
  double tmp;
  char *keyname;
  gal_fits_list_key_t *headers=NULL;

  /* Add the proper information to the header of the output: */
  gal_fits_key_write_filename("INF", p->inputname, &headers, 0);
  if( !isnan(p->background) )
    {
      gal_checkset_allocate_copy("BCKGRND", &keyname);
      gal_fits_key_list_add_end(&headers, GAL_TYPE_FLOAT64, keyname, 1,
                                &p->background, 0,
                                "Background value for Poisson noise.",
                                0, NULL, 0);
      if( !isnan(p->zeropoint) )
        {
          tmp=-2.5 * log10(p->background) + p->zeropoint;
          gal_checkset_allocate_copy("BCKGMAG", &keyname);
          gal_fits_key_list_add_end(&headers, GAL_TYPE_FLOAT64, keyname, 1,
                                    &tmp, 0,
                                    "Background value in magnitudes",
                                    0, NULL, 0);
          gal_checkset_allocate_copy("BCKGZP", &keyname);
          gal_fits_key_list_add_end(&headers, GAL_TYPE_FLOAT64, keyname, 1,
                                    &p->zeropoint, 0,
                                    "Zeropoint for interpreting magnitudes.",
                                    0, NULL, 0);
        }
      if( !isnan(p->instrumental) )
        {
          gal_checkset_allocate_copy("INSTRU", &keyname);
          gal_fits_key_list_add_end(&headers, GAL_TYPE_FLOAT64, keyname, 1,
                                    &p->instrumental, 0,
                                    "Instrumental noise in units of flux.",
                                    0, NULL, 0);
        }
    }
  else
    {
      gal_checkset_allocate_copy("SIGMA", &keyname);
      gal_fits_key_list_add_end(&headers, GAL_TYPE_FLOAT64, keyname, 1,
                                &p->sigma, 0, "Total noise sigma", 0,
                                NULL, 0);
    }
  gal_checkset_allocate_copy("RNGTYPE", &keyname);
  gal_fits_key_list_add_end(&headers, GAL_TYPE_STRING, keyname, 1,
                            (void *)(p->rng_name), 0,
                            "Random number generator (by GSL) type.",
                            0, NULL, 0);
  gal_checkset_allocate_copy("RNGSEED", &keyname);
  gal_fits_key_list_add_end(&headers, GAL_TYPE_ULONG, keyname, 1,
                            &p->rng_seed, 0,
                            "Random number generator (by GSL) seed.",
                            0, NULL, 0);

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
  double *d, *df, back=p->background;
  double inst = ( isnan(p->instrumental)
                  ? 0.0f
                  : p->instrumental*p->instrumental );

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
        *d += back + gsl_ran_gaussian(p->rng, sqrt( inst + back + *d ));
      while(++d<df);
    }

  /* Convert and save the output in the proper format: */
  convertsaveoutput(p);
}
