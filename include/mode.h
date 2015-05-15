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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef MODE_H
#define MODE_H

#define MODELOWQUANTILE  0.01f
#define MODEHIGHQUANTILE 0.51f

#define MODESYMGOOD       0.2f
#define MODELOWQUANTGOOD  0.02f

#define SYMMETRICITYLOWQUANT 0.01f

#define GOLDENRATIO        1.618034f
#define TWOTAKEGOLDENRATIO 0.38197f

#define MIRRORISABOVERESULT    (size_t)(-1)

struct modeparams
{
  float     *sorted;		/* Sorted array to be used.                */
  size_t       size;		/* Number of elements in the sorted array. */
  size_t       lowi;		/* Lower quantile of interval.             */
  size_t       midi;		/* Middle quantile of interval.            */
  size_t       midd;		/* Middle index of inteval.                */
  size_t      highi;		/* Higher quantile of interval.            */
  float   tolerance;		/* Tolerance level to terminate search.    */
  size_t   numcheck;		/* Number of pixels after mode to check.   */
  size_t   interval;		/* Interval to check pixels.               */
  float   errorstdm;		/* Multiple of standard deviation.         */
};

void
makemirrorplots(float *sorted, size_t size, size_t mirrorindex,
                float min, float max, size_t numbins, char *histsname,
                char *cfpsname, float mirrorplotdist);

float
valuefromsym(float *sorted, size_t size, size_t modeindex, float sym);

void
modeindexinsorted(float *sorted, size_t size, float errorstdm,
                  size_t *modeindex, float *modesym);

#endif
