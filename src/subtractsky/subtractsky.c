/*********************************************************************
SubtractSky - Find and subtract the sky value from an image.
SubtractSky is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"
#include "timing.h"
#include "fitsarrayvv.h"

#include "main.h"

void
subtractsky(struct subtractskyparams *p)
{
  long *meshindexs;
  float *sky, *std;
  struct timeval t1;
  struct meshparams *mp=&p->mp;


  /* Prepare the mesh array. */
  gettimeofday(&t1, NULL);
  makemesh(mp);
  if(p->meshname)
    {
      checkmeshid(mp, &meshindexs);
      arraytofitsimg(p->meshname, "Input", FLOAT_IMG, p->mp.img,
                     mp->s0, mp->s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      arraytofitsimg(p->meshname, "MeshIndexs", LONG_IMG, meshindexs,
                     mp->s0, mp->s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      free(meshindexs);
    }
  if(p->cp.verb) reporttiming(&t1, "Mesh grid ready.", 1);



  /* Find the sky value and its standard deviation on each mesh. */
  fillmesh(mp, MODEEQMED_AVESTD, 0);
  if(p->interpname)
    {
      checkgarray(mp, MODEEQMED_AVESTD, &sky, &std);
      arraytofitsimg(p->interpname, "Sky", FLOAT_IMG, sky,
                     mp->s0, mp->s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      arraytofitsimg(p->interpname, "STD", FLOAT_IMG, std,
                     mp->s0, mp->s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      free(sky); free(std);
    }
  if(p->cp.verb)
    reporttiming(&t1, "Sky and its STD found on some meshes.", 1);


  /* Interpolate over the meshs to fill all the blank ones in both the
     sky and the standard deviation arrays: */
  meshinterpolate(mp, 2);
  if(p->interpname)
    {
      checkgarray(mp, MODEEQMED_AVESTD, &sky, &std);
      arraytofitsimg(p->interpname, "Sky", FLOAT_IMG, sky,
                     mp->s0, mp->s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      arraytofitsimg(p->interpname, "STD", FLOAT_IMG, std,
                     mp->s0, mp->s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      free(sky); free(std);
    }
  if(p->cp.verb) reporttiming(&t1, "All blank meshs filled.", 1);


  /* Smooth the interpolated array:  */



  /* Subtract the sky value */

}
