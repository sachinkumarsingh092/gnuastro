/*********************************************************************
NoiseChisel - Detect and segment signal in noise.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include "timing.h"

#include "main.h"

#include "thresh.h"
#include "noisechisel.h"










/*********************************************************************/
/******************         NoiseChisel        ***********************/
/*********************************************************************/
void
noisechisel(struct noisechiselparams *p)
{
  struct meshparams *smp=&p->smp, *lmp=&p->lmp;

  long *meshindexs;
  struct timeval t1;
  size_t s0=smp->s0, s1=smp->s1;


  /* Prepare the mesh array. */
  gettimeofday(&t1, NULL);
  makemesh(smp);
  makemesh(lmp);
  if(p->meshname)
    {
      arraytofitsimg(p->meshname, "Input", FLOAT_IMG, smp->img, s0, s1,
                     p->numblank, p->wcs, NULL, SPACK_STRING);
      checkmeshid(smp, &meshindexs);
      arraytofitsimg(p->meshname, "SmallMeshIndexs", LONG_IMG, meshindexs,
                     s0, s1, 0, p->wcs, NULL, SPACK_STRING);
      free(meshindexs);
      checkmeshid(lmp, &meshindexs);
      arraytofitsimg(p->meshname, "LargeMeshIndexs", LONG_IMG, meshindexs,
                     s0, s1, 0, p->wcs, NULL, SPACK_STRING);
      free(meshindexs);
    }
  if(p->cp.verb) reporttiming(&t1, "Mesh grids ready.", 1);


  /* Convolve the image: */
  gettimeofday(&t1, NULL);
  spatialconvolveonmesh(smp, &p->conv);
  if(p->detectionname)
    {
      arraytofitsimg(p->detectionname, "Input", FLOAT_IMG, smp->img,
                     s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(p->detectionname, "Convolved", FLOAT_IMG, p->conv,
                     s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);
    }
  if(p->cp.verb) reporttiming(&t1, "Convolved with kernel.", 1);


  /* Find the threshold and apply it: */
  gettimeofday(&t1, NULL);
  findapplythreshold(p);
  if(p->detectionname)
    arraytofitsimg(p->detectionname, "Thresholded", BYTE_IMG, p->byt,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(p->cp.verb)
    reporttiming(&t1, "Quantile threshold found and applied.", 1);

  /* Clean up: */
  free(p->byt);
  free(p->conv);
  freemesh(smp);
  freemesh(lmp);
}
