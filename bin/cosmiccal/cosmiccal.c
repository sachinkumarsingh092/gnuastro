/*********************************************************************
CosmicCalculator - Calculate cosmological parameters
CosmicCalculator is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2016, Free Software Foundation, Inc.

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
#include <stdlib.h>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_integration.h>

#include "main.h"

#include "cosmiccal.h"


/**************************************************************/
/************         Integrand functions         *************/
/**************************************************************/
/* In these functions, z as a separate argument, this is not necessarily
   the same z as the redshift in cosmiccalparams. */

double
Ez(double z, void *params)
{
  struct cosmiccalparams *p=(struct cosmiccalparams *)params;
  return sqrt( p->olambda
               + p->ocurv           * (1+z) * (1+z)
               + p->omatter         * (1+z) * (1+z) * (1+z)
               + p->oradiation      * (1+z) * (1+z) * (1+z) * (1+z));
}





double
age(double z, void *params)
{
  return 1 / ( (1+z)*Ez(z,params) );
}





double
propdist(double z, void *params)
{
  return 1 / ( Ez(z,params) );
}




















/**************************************************************/
/************             Integrators             *************/
/**************************************************************/
/* Estimate the age of the universe, note that z might be different
   from the desired redshift. */
double
ageofuniverse(struct cosmiccalparams *p, double z)
{
  gsl_function F;
  double result, error;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(GSLILIMIT);

  /* Set the GSL function parameters */
  F.params=p;
  F.function=&age;

  gsl_integration_qagiu(&F, z, GSLIEPSABS, GSLIEPSREL, GSLILIMIT, w,
                        &result, &error);

  return result / p->H0s / (365*GSL_CONST_MKSA_DAY) / 1e9;
}





/* Return the proper distance to a source at z in units of mega
   parsecs */
double
properdistance(struct cosmiccalparams *p, double z)
{
  size_t neval;
  gsl_function F;
  double result, error;

  /* Set the GSL function parameters */
  F.params=p;
  F.function=&propdist;

  gsl_integration_qng(&F, 0.0f, z, GSLIEPSABS, GSLIEPSREL,
                      &result, &error, &neval);

  return result * p->c / p->H0s / (1e6 * GSL_CONST_MKSA_PARSEC);
}





double
covolume(double z, void *params)
{
  size_t neval;
  gsl_function F;
  double result, error;

  /* Set the GSL function parameters */
  F.params=params;
  F.function=&propdist;

  gsl_integration_qng(&F, 0.0f, z, GSLIEPSABS, GSLIEPSREL,
                      &result, &error, &neval);

  return result * result / ( Ez(z,params) );
}





double
comovingvolume(struct cosmiccalparams *p, double z)
{
  size_t neval;
  gsl_function F;
  double result, error;
  double cH = p->c / p->H0s / (1e6 * GSL_CONST_MKSA_PARSEC);

  /* Set the GSL function parameters */
  F.params=p;
  F.function=&covolume;

  gsl_integration_qng(&F, 0.0f, z, GSLIEPSABS, GSLIEPSREL,
                      &result, &error, &neval);

  return result * 4 * M_PI * cH*cH*cH;
}




















/**************************************************************/
/************        Intermediary functions       *************/
/**************************************************************/
/* Critical density at redshift z in units of gram/cm^3*/
double
criticaldensity(struct cosmiccalparams *p, double z)
{
  double H = p->H0s*Ez(z,p);
  return 3*H*H/(8*M_PI*GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT)/1000;
}

























/**************************************************************/
/************            Main function            *************/
/**************************************************************/
void
cosmiccal(struct cosmiccalparams *p)
{
  double ad, ld, vz, pd, absmagconv;
  double curage, ccritd, distmod, outage, zcritd;

  /* In case the user just wants one number, only print that and
     return. */
  if(p->onlyvolume){
    printf("%f\n", comovingvolume(p,p->redshift));
    return;
  }
  if(p->onlyabsmagconv){
    pd=properdistance(p, p->redshift);
    ld=pd*(1+p->redshift);
    distmod=5*(log10(ld*1000000)-1);
    absmagconv=distmod-2.5*log10(1+p->redshift);
    printf("%f\n", absmagconv);
    return;
  }

  /* The user wants everything, do all the calculations and print
     everything with full descriptions. */
  curage=ageofuniverse(p, 0.0f);
  ccritd=criticaldensity(p, 0.0f);
  vz=comovingvolume(p,p->redshift);
  pd=properdistance(p, p->redshift);
  outage=ageofuniverse(p, p->redshift);
  zcritd=criticaldensity(p, p->redshift);

  ad=pd/(1+p->redshift);
  ld=pd*(1+p->redshift);
  distmod=5*(log10(ld*1000000)-1);
  absmagconv=distmod-2.5*log10(1+p->redshift);

  /* Print out results: */
  printf("%s\n", SPACK_STRING);
  printf("\n Input parameters\n");
  printf(  " ----------------\n");
  printf(FLTFORMAT, "Desired redshift for calculations (z):", p->redshift);
  printf(FLTFORMAT, "Expansion rate (Hubble constant, H0), now:", p->H0);
  printf(FLTFORMAT, "Cosmological constant fractional density, now:",
         p->olambda);
  printf(FLTFORMAT, "Matter fractional density, now:", p->omatter);
  printf(EXPFORMAT, "Radiation fractional density, now:", p->oradiation);
  printf(EXPFORMAT, "Curvatue fractional density (from the above):",
         p->ocurv);


  printf("\n\n Universe now\n");
  printf(    " ------------\n");
  printf(FLTFORMAT, "Age of Universe now (Gyr):", curage);
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
  printf(FLTFORMAT, "Age of Universe at z (Gyr):", outage);
  printf(FLTFORMAT, "Look-back time to z (Gyr):", curage-outage);
  printf(EXPFORMAT, "Critical density at z (g/cm^3):",  zcritd);

  printf("\n\n Comoving universe (time independent)\n");
  printf(    " ------------------------------------\n");
  printf(FLTFORMAT, "Comoving volume over 4pi stradian to z (Mpc^3):", vz);
}
