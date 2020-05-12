/*********************************************************************
CosmicCalculator - Calculate cosmological parameters
CosmicCalculator is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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

#include <gnuastro/cosmology.h>

#include "main.h"

#include "ui.h"
#include "cosmiccal.h"











/**************************************************************/
/************            Main function            *************/
/**************************************************************/
static void
cosmiccal_print_input(struct cosmiccalparams *p)
{
  printf("%s\n", PROGRAM_STRING);
  printf("\n Input parameters\n");
  printf(  " ----------------\n");
  if( !isnan(p->redshift) )
    printf(FLTFORMAT, "Desired redshift for calculations (z):", p->redshift);
  printf(FLTFORMAT, "Expansion rate (Hubble constant, H0), now:", p->H0);
  printf(FLTFORMAT, "Cosmological constant fractional density, now:",
         p->olambda);
  printf(FLTFORMAT, "Matter fractional density, now:", p->omatter);
  printf(EXPFORMAT, "Radiation fractional density, now:", p->oradiation);
  printf(EXPFORMAT, "Curvature fractional density (from the above):",
         1 - ( p->olambda + p->omatter + p->oradiation ));
}





static void
cosmiccal_printall(struct cosmiccalparams *p)
{
  double ad, ld, vz, pd, absmagconv;
  double curage, ccritd, distmod, outage, zcritd;

  /* The user wants everything, do all the calculations and print
     everything with full descriptions. */
  curage=gal_cosmology_age(0.0f, p->H0, p->olambda, p->omatter,
                           p->oradiation);

  ccritd=gal_cosmology_critical_density(0.0f, p->H0, p->olambda, p->omatter,
                                        p->oradiation);

  pd=gal_cosmology_proper_distance(p->redshift, p->H0, p->olambda, p->omatter,
                                   p->oradiation);

  ad=gal_cosmology_angular_distance(p->redshift, p->H0, p->olambda,
                                    p->omatter, p->oradiation);

  ld=gal_cosmology_luminosity_distance(p->redshift, p->H0, p->olambda,
                                       p->omatter, p->oradiation);

  distmod=gal_cosmology_distance_modulus(p->redshift, p->H0, p->olambda,
                                         p->omatter, p->oradiation);

  absmagconv=gal_cosmology_to_absolute_mag(p->redshift, p->H0, p->olambda,
                                           p->omatter, p->oradiation);

  outage=gal_cosmology_age(p->redshift, p->H0, p->olambda, p->omatter,
                           p->oradiation);

  zcritd=gal_cosmology_critical_density(p->redshift, p->H0, p->olambda,
                                        p->omatter, p->oradiation);

  vz=gal_cosmology_comoving_volume(p->redshift, p->H0, p->olambda, p->omatter,
                                   p->oradiation);

  /* Print out results: */
  cosmiccal_print_input(p);


  printf("\n\n Universe now\n");
  printf(    " ------------\n");
  printf(FLTFORMAT, "Age of Universe now (Ga*):", curage);
  printf(EXPFORMAT, "Critical density now (g/cm^3):",  ccritd);
  printf(FLTFORMAT, "Proper distance to z (Mpc):", pd);
  printf(FLTFORMAT, "Angular diameter distance to z (Mpc):", ad);
  printf(FLTFORMAT, "Tangential distance covered by 1 arcsec at z (Kpc):",
         ad*1000*M_PI/3600/180);
  printf(FLTFORMAT, "Luminosity distance to z (Mpc):", ld);
  printf(FLTFORMAT, "Distance modulus at z (no unit):", distmod);
  printf(FLTFORMAT, "Conversion to absolute magnitude (no unit):",
         absmagconv);


  printf("\n\n Universe at desired redshift z\n");
  printf(    " ------------------------------\n");
  printf(FLTFORMAT, "Age of Universe at z (Ga*):", outage);
  printf(FLTFORMAT, "Look-back time to z (Ga*):", curage-outage);
  printf(EXPFORMAT, "Critical density at z (g/cm^3):",  zcritd);

  printf("\n\n Comoving universe (time independent)\n");
  printf(    " ------------------------------------\n");
  printf(FLTFORMAT, "Comoving volume over 4pi stradian to z (Mpc^3):", vz);

  printf("\n-------\n");
  printf("*: Ga is short for Giga Annum, or billion years (IAU standard).\n");
}





void
cosmiccal(struct cosmiccalparams *p)
{
  gal_list_i32_t *tmp;
  double curage, zage;

  /* If no redshift is given, just print the input parameters along with a
     notice that further calculations are only possible with a redshift and
     abort. */
  if(isnan(p->redshift))
    {
      cosmiccal_print_input(p);
      printf("\n\nPlease specify a redshift with the '--redshift' (or '-z') "
             "option.\n");
      return;
    }

  /* In case the user just wants one number, only print that and
     return. */
  if(p->specific)
    {
      for(tmp=p->specific;tmp!=NULL;tmp=tmp->next)
        {
          switch(tmp->v)
            {
            case UI_KEY_USEDREDSHIFT:
              printf("%g",
                     p->redshift==MAIN_REDSHIFT_ZERO ? 0.0f: p->redshift);
              break;

            case UI_KEY_AGENOW:
              printf("%f", gal_cosmology_age(0.0f, p->H0, p->olambda,
                                              p->omatter, p->oradiation));
              break;

            case UI_KEY_CRITICALDENSITYNOW:
              printf("%e", gal_cosmology_critical_density(0.0f, p->H0,
                                                           p->olambda,
                                                           p->omatter,
                                                           p->oradiation));
              break;

            case UI_KEY_PROPERDISTANCE:
              printf("%f", gal_cosmology_proper_distance(p->redshift, p->H0,
                                                          p->olambda,
                                                          p->omatter,
                                                          p->oradiation));
              break;

            case UI_KEY_ANGULARDIMDIST:
              printf("%f", gal_cosmology_angular_distance(p->redshift, p->H0,
                                                           p->olambda,
                                                           p->omatter,
                                                           p->oradiation));
              break;

            case UI_KEY_ARCSECTANDIST:
              printf("%f", ( gal_cosmology_angular_distance(p->redshift, p->H0,
                                                             p->olambda,
                                                             p->omatter,
                                                             p->oradiation)
                              * 1000 * M_PI / 3600 / 180 ) );
              break;

            case UI_KEY_LUMINOSITYDIST:
              printf("%f", gal_cosmology_luminosity_distance(p->redshift,
                                                              p->H0,
                                                              p->olambda,
                                                              p->omatter,
                                                              p->oradiation));
              break;

            case UI_KEY_DISTANCEMODULUS:
              printf("%f", gal_cosmology_distance_modulus(p->redshift, p->H0,
                                                           p->olambda,
                                                           p->omatter,
                                                           p->oradiation));
              break;

            case UI_KEY_ABSMAGCONV:
              printf("%f", gal_cosmology_to_absolute_mag(p->redshift, p->H0,
                                                          p->olambda,
                                                          p->omatter,
                                                          p->oradiation));
              break;

            case UI_KEY_AGE:
              printf("%f", gal_cosmology_age(p->redshift, p->H0, p->olambda,
                                              p->omatter, p->oradiation));
              break;

            case UI_KEY_LOOKBACKTIME:
              curage=gal_cosmology_age(0.0f, p->H0, p->olambda, p->omatter,
                                       p->oradiation);
              zage=gal_cosmology_age(p->redshift, p->H0, p->olambda, p->omatter,
                                     p->oradiation);
              printf("%f", curage-zage);
              break;

            case UI_KEY_CRITICALDENSITY:
              printf("%e", gal_cosmology_critical_density(p->redshift, p->H0,
                                                           p->olambda,
                                                           p->omatter,
                                                           p->oradiation));
              break;

            case UI_KEY_VOLUME:
              printf("%f", gal_cosmology_comoving_volume(p->redshift, p->H0,
                                                          p->olambda,
                                                          p->omatter,
                                                          p->oradiation));
              break;

            case UI_KEY_LINEATZ:
              printf("%g", gal_list_f64_pop(&p->specific_arg)*(1+p->redshift));
              break;

            default:
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
                    "fix the problem. The code %d is not recognized as a "
                    "single value calculation code", __func__,
                    PACKAGE_BUGREPORT, tmp->v);
            }

          /* Only add a space-character if there are more results to print. */
          if(tmp->next) printf(" ");
        }

      /* Print a new-line character to finish the output. */
      printf("\n");
    }
  else
    cosmiccal_printall(p);
}
