/*********************************************************************
mode -- Find the mode of a distribution.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <float.h>
#include <stdlib.h>

#include <gnuastro/array.h>
#include <gnuastro/statistics.h>

#include "mode.h"


/****************************************************************
 *****************        Mode plots         ********************
 ****************************************************************/
/*
This is the python code you can put in `plot.py` so the plot
command mentioned here works:

-----------------------------------------------------------------
#! /usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

a=np.loadtxt(sys.argv[1])
b=np.loadtxt(sys.argv[2])
w=a[1,0]-a[0,0]

plt.bar(a[:,0], a[:,1], width=w, color="blue", linewidth=0, alpha=0.6)
plt.bar(a[:,0], a[:,2], width=w, color="green", linewidth=0, alpha=0.4)
plt.plot(b[:,0], b[:,1], linewidth=2, color="blue")
plt.plot(b[:,0], b[:,2], linewidth=2, color="green")

plt.ylim([0,np.amax(a[:,1])])
plt.xlim([np.amin(a[:,0]),np.amax(a[:,0])])

plt.savefig(sys.argv[3])
---------------------------------------------------------

Run it with `./plot.py histsname.txt cfpsname.txt outputpdfname.pdf`
It will plot the corresponding histograms and cumulative frequency
plots. If you like to, this call is available as a system() call in
the functions below.
*/

/* This is used for the plots, it will allocate an array and put the
   mirrored array values in it. `mi` is the index the mirror is to
   be placed on.  */
void
makemirrored(float *in, size_t mi, float **outmirror, size_t *outsize)
{
  size_t i, size;
  float *mirror, zf;

  zf=in[mi];
  size=2*mi+1;

  errno=0;
  mirror=malloc(size*sizeof *mirror);
  if(mirror==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for mirror array in "
          "makemirrored (mode.c)", size*sizeof *mirror);

  for(i=0;i<=mi;++i)
    mirror[i]    = in[i];
  for(i=1;i<=mi;++i)
    mirror[mi+i] = 2*zf - mirror[mi-i];

  *outmirror=mirror;
  *outsize=size;
}




















/****************************************************************
 *****************           Mode            ********************
 ****************************************************************/

/*Find the index where the two CDFs dirverge beyond CDFSDIVERGEDIFF.
  T he input CDF and a CDF of a mirrored distribution about the
  quantile `mirrorquant` are compared.

  The basic idea behind finding the mode is comparing the mirrored CDF
  (about a test for the mode) with the actual CDF for a given
  point. This is the job of this function. It takes the ordered array
  and the quantile that is to be checked, it then finds the maximum
  difference between the mirrored CDF about the given point and the
  input CDF.

  `zf` keeps the flux at the mirror (zero) point.  `i` is used to
  count the pixels before m. So `m+i` is the index of the mirrored
  distribution and mf=zf+(zf-a[m-i])=2*zf-a[m-i] is the mirrored flux
  at this point. `j` is found such that a[m+j] has the nearest flux to
  `mf`.

  The desired difference between the input CDF and the mirrored one
  for each `i` is then simply: `j-i`.

  Once `i` is incremented, `mf` will increase, so to find the new `j`
  we don't need to begin looking from `j=0`. Remember that the array
  is sorted, so the desired `j` is definitely larger than the previous
  `j`. So, if we keep the previous `j` in `prevj` then, all we have to
  do is to start incrementing `j` from `prevj`. This will really help
  in speeding up the job :-D. Only for the first element, `prevj=0`.*/
size_t
mirrormaxdiff(float *a, size_t size, size_t m,
              size_t numcheck, size_t interval, size_t stdm)
{
  /* The variables:
   i:        Index on mirror distribution.
   j:        Index on input distribution.
   prevj:    Index of previously checked point in the actual array.
   mf:       Flux that is approximately equal in both distributions.*/
  float mf, zf;
  size_t  maxdiff=0, errordiff;
  size_t i, j, absdiff, prevj=0;

  /* Find the error and mirror flux value  */
  zf=a[m];
  errordiff=stdm*sqrt(m);

  /*
  printf("###############\n###############\n");
  printf("### Mirror pixel: %zu\n", m);
  printf("###############\n###############\n");
  */
  /* Go over the mirrored points. */
  for(i=1; i<numcheck && i<=m && m+i<size ;i+=interval)
    {
      mf=2*zf-a[m-i];

      /* The moment a[m+j]>mf, we have reached the last pixel to
         check. Now we just have to see if a[m+j-1] is closer to mf or
         if a[m+j]. We change `j` accordingly and break out of the `j`
         loop. */
      for(j=prevj;j<size-m;++j)
        if(a[m+j]>mf)
          {
            if( a[m+j]-mf < mf-a[m+j-1] )
              break;
            else
              {
                j--;
                break;
              }
          }
      /*
      printf("i:%-5zu j:%-5zu diff:%-5d maxdiff: %zu\n",
             i, j, (int)j-(int)i, maxdiff);
      */
      /* The index of the actual CDF corresponding the the mirrored
         flux has been found. We want the mirrored distribution to be
         within the actual distribution, not beyond it, so the only
         acceptable results are when i<j. If i>j+errordiff then the
         result is not acceptable! */
      if(i>j+errordiff)
        {
          maxdiff = MODE_MIRROR_IS_ABOVE_RESULT;
          break;
        }
      absdiff  = i>j ? i-j : j-i;
      if(absdiff>maxdiff)
        maxdiff=absdiff;

      prevj=j;
    }
  return maxdiff;
}





/* Find the mode of a float array of size `size`. I assume that
   mirrormaxdiff() has one minimum (within the statistical errors) in the
   function. To find that minimum, the golden section search algorithm is
   going to used. Read the Wikipedia article for a very nice
   introduction. In summary we will constantly be finding middle points in
   the given interval and thus decreasing the interval until a certain
   tolerance is reached.

   If the input interval is on points `a` and `b`, then the middle point
   (lets call it `c`, where c>a and c<b) to test should be positioned such
   that (b-c)/(c-a)=GAL_MODE_GOLDEN_RATIO. Once we open up this relation,
   we can find c using:

      c=(b+GAL_MODE_GOLDEN_RATIO*a)/(1+GAL_MODE_GOLDEN_RATIO).

   We need a fourth point to be placed between. With this configuration,
   the probing point is located at: */
size_t
modegoldenselection(struct gal_statistics_mode_params *mp)
{
  size_t di, dd;
  /*static int counter=1;*/
  /*------------------------------------------------------------------
  char outname[500], command[1000];
  char histsname[500], cfpsname[500];
  ------------------------------------------------------------------*/

  /* Find the probing point in the larger interval. */
  if(mp->highi-mp->midi > mp->midi-mp->lowi)
    di = mp->midi + MODE_TWO_TAKE_GOLDEN_RATIO *(float)(mp->highi-mp->midi);
  else
    di = mp->midi - MODE_TWO_TAKE_GOLDEN_RATIO * (float)(mp->midi-mp->lowi);

  /* Since these are all indexs (and positive) we don't need an
     absolute value, highi is also always larger than lowi! In some
     cases, the first (standard) condition might be satisfied, while
     highi-lowi<=2. In such cases, also jump out! */
  if( (mp->highi - mp->lowi) < mp->tolerance*(mp->midi+di)
      || (mp->highi - mp->lowi) <= 3)
    return (mp->highi+mp->lowi)/2;

  /* Find the maximum difference for this quantile. */
  dd = mirrormaxdiff(mp->sorted, mp->size, di, mp->numcheck,
                     mp->interval, mp->errorstdm);

  /*------------------------------------------------------------------
  sprintf(outname, "%dcmp.pdf", counter);
  sprintf(cfpsname, "%dcfps.txt", counter);
  sprintf(histsname, "%dhists.txt", counter);
  gal_mode_make_mirror_plots(mp->sorted, mp->size, di, histsname, cfpsname);
  sprintf(command, "./plot.py %s %s %s", histsname, cfpsname, outname);
  system(command);
  -------------------------------------------------------------------*/
  /*
  printf("%-5zu\t%-5zu(%d)\t%-5zu ----> dq: %-5zu di: %d\n",
         mp->lowi, mp->midi, (int)mp->midd, mp->highi,
         di, (int)dd);
  */
  /* +++++++++++++ The mirrored distribution's cumulative frequency plot
     has be lower than the actual's cfp. If it isn't, `di` will be
     GAL_MODE_MIRROR_IS_ABOVE_RESULT. In this case, the normal golden
     section minimization is not going to give us what we want. So I have
     added this modification to it. In such cases, we want the search to go
     to the lower intervals.*/
  if(dd==MODE_MIRROR_IS_ABOVE_RESULT)
    {
      if(mp->midi < di)
        {
          mp->highi=di;
          return modegoldenselection(mp);
        }
      else
        {
          mp->highi=mp->midi;
          mp->midi=di;
          mp->midd=dd;
          return modegoldenselection(mp);
        }
    }
  /* +++++++++++++ End of my addition to the golden section search. */

  /* This is the standard golden section search: */
  if(dd<mp->midd)
    {
      if(mp->highi-mp->midi > mp->midi-mp->lowi)
        {
          mp->lowi  = mp->midi;
          mp->midi  = di;
          mp->midd  = dd;
          return modegoldenselection(mp);
        }
      else
        {
          mp->highi = mp->midi;
          mp->midi  = di;
          mp->midd  = dd;
          return modegoldenselection(mp);
        }
    }
  else
    {
      if(mp->highi-mp->midi > mp->midi-mp->lowi)
        {
          mp->highi = di;
          return modegoldenselection(mp);
        }
      else
        {
          mp->lowi  = di;
          return modegoldenselection(mp);
        }
    }
}





/* Once the mode is found, we need to do a quality control. This
   quality control is the measure of symmetricity. Lets assume the
   mode index is at `m`, the error in `m` can be assumed to be
   sqrt(m). So lets call the first point that the difference between
   the cumulative distribution of the mirror and actual data deviate
   above sqrt(m), is at index `b`. For a scale parameter, lets assume
   that the index of 5% of `m` is `a`. We could have taken the
   distribution minimum, but the scatter in that can be too high! Now
   the symmetricity of the mode can be quantified as: (b-m)/(m-a). For
   a completly symmetric mode, this should be 1. Note that the search
   for `b` only goes to the 95% of the distribution.  */
void
modesymmetricity(float *a, size_t size, size_t mi, float errorstdm,
                 float *sym)
{
  float af, bf, mf, fi;
  size_t i, j, bi=0, topi, errdiff, prevj=0;

  mf=a[mi];
  errdiff=errorstdm*sqrt(mi);
  topi = 2*mi>size-1 ? size-1 : 2*mi;
  af=a[gal_statistics_quantile_index(2*mi+1,
                          GAL_STATISTICS_MODE_SYMMETRICITY_LOW_QUANT)];

  /* This loop is very similar to that of mirrormaxdiff(). It will
     find the index where the difference between the two cumulative
     frequency plots exceeds that of the error in the mirror index. */
  for(i=1; i<topi-mi ;i+=1)
    {
      fi=2*mf-a[mi-i];

      for(j=prevj;j<size-mi;++j)
        if(a[mi+j]>fi)
          {
            if( a[mi+j]-fi < fi-a[mi+j-1] )
              break;
            else
              {
                j--;
                break;
              }
          }

      if(i>j+errdiff || j>i+errdiff)
        {
          bi=mi+i;
          break;
        }
      prevj=j;
    }

  /* bi==0 shows that no point with a larger difference could be
     found. So bi should be set to the end of the search region. */
  if(bi==0) bi=topi;

  bf=a[bi];
  /*
  printf("%f, %f, %f\n", af, mf, bf);
  */

  *sym=(bf-mf)/(mf-af);

  /* This is mainly used for plotting, which subtracts `mf`.
  printf("SymmetricFlux: %f\n", a[bi]-mf);
  printf("symmetricity: %f\n", sym);
  */
}
