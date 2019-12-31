/*********************************************************************
MakeProfiles - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_gamma.h>   /* For total Sersic brightness. */

#include "main.h"
#include "mkprof.h"             /* Needs main.h, astrthreads.h */
#include "profiles.h"





/****************************************************************
 *****************         Profiles:         ********************
 ****************************************************************/
/* The Gaussian function at a point. */
double
profiles_radial_distance(struct mkonthread *mkp)
{
  return mkp->r;
}





/* The integral of the Gaussian from -inf to +inf equals the square root of
   PI. So from zero to +inf it equals half of that.*/
double
profiles_gaussian_total(double q)
{
  return q*sqrt(M_PI)/2;
}





/* The Gaussian function at a point. */
double
profiles_gaussian(struct mkonthread *mkp)
{
  return exp( mkp->gaussian_c * mkp->r * mkp->r );
}





/* This function will find the moffat function alpha value based on
   the explantions here:

   http://labs.adsabs.harvard.edu/adsabs/abs/2001MNRAS.328..977T/

   alpha=(FWHM/2)/(2^(1.0/beta)-1)^(0.5). Then the moffat
   function at r is: (1.0 + (r/alpha)^2.0)^(-1.0*beta)*/
double
profiles_moffat_alpha(double fwhm, double beta)
{
  return (fwhm/2)/pow((pow(2, 1/beta)-1), 0.5f);
}





/* Find the total value of the Moffat profile. I am using equation 10
 from Pengetal 2010 (Galfit). In finding the profiles, I am assuming
 \Sigma_0=1. So that is what I put here too.*/
double
profiles_moffat_total(double alpha, double beta, double q)
{
  return M_PI*alpha*alpha*q/(beta-1);
}





/* Find the Moffat profile for a certain radius.

   rda=r/alpha     and nb=-1*b.

   This is done before hand to speed up the process. */
double
profiles_moffat(struct mkonthread *mkp)
{
  return pow(1+mkp->r*mkp->r/mkp->moffat_alphasq, mkp->moffat_nb);
}





/* This approximation of b(n) for n>0.35 is taken from McArthur,
   Courteau and Holtzman 2003:
   http://adsabs.harvard.edu/abs/2003ApJ...582..689 */
double
profiles_sersic_b(double n)
{
  if(n<=0.35f)
    error(EXIT_FAILURE, 0, "the Sersic index cannot be smaller "
          "than 0.35. It is %.3f", n);
  return 2*n-(1/3)+(4/(405*n))+(46/(25515*n*n))+
    (131/(1148175*n*n*n)-(2194697/(30690717750*n*n*n*n)));
}





/* Find the total brightness in a Sersic profile. From equation 4 in
   Peng 2010. This assumes the surface brightness at the effective
   radius is 1.*/
double
profiles_sersic_total(double n, double re, double b, double q)
{
  return (2*M_PI*re*re*exp(b)*n*pow(b, -2*n)*q*
          gsl_sf_gamma(2*n));
}





/* Find the Sersic profile for a certain radius. rdre=r/re, inv_n=1/n,
   nb= -1*b.  */
double
profiles_sersic(struct mkonthread *mkp)
{
  return exp( mkp->sersic_nb
              * ( pow(mkp->r/mkp->sersic_re, mkp->sersic_inv_n) -1 ) );
}





/* Make a circumference (inner to the radius). */
double
profiles_circumference(struct mkonthread *mkp)
{
  return ( (mkp->r > mkp->intruncr && mkp->r <= mkp->truncr)
           ? mkp->fixedvalue : 0.0f );
}





/* Always returns a fixed value: */
double
profiles_flat(struct mkonthread *mkp)
{
  return mkp->fixedvalue;
}
