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
#include "binary.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "label.h"
#include "thresh.h"
#include "detection.h"
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
  int verb=p->cp.verb;
  char report[VERBMSGLENGTH_V];
  size_t i, s0=smp->s0, s1=smp->s1, numlabs;
  char *initdetectionname=p->initdetectionname;


  /* Prepare the mesh array. */
  if(verb) gettimeofday(&t1, NULL);
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
  if(verb) reporttiming(&t1, "Mesh grids ready.", 1);



  /* Convolve the image: */
  if(verb) gettimeofday(&t1, NULL);
  spatialconvolveonmesh(smp, &p->conv);
  if(initdetectionname)
    {
      arraytofitsimg(initdetectionname, "Input", FLOAT_IMG, smp->img,
                     s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(initdetectionname, "Convolved", FLOAT_IMG, p->conv,
                     s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);
    }
  if(verb) reporttiming(&t1, "Convolved with kernel.", 1);



  /* Find the threshold and apply it: */
  if(verb) gettimeofday(&t1, NULL);
  findapplythreshold(p);
  if(initdetectionname)
    arraytofitsimg(initdetectionname, "Thresholded", BYTE_IMG, p->byt,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(verb)
    {
      sprintf(report, "%.2f quantile threshold found and applied.",
              p->qthresh);
      reporttiming(&t1, report, 1);
    }



  /* Erode the thresholded image: */
  if(verb) gettimeofday(&t1, NULL);
  if(p->erodengb==4)
    for(i=0;i<p->numerosion;++i)
      dilate0_erode1_4con(p->byt, s0, s1, 1);
  else
    for(i=0;i<p->numerosion;++i)
      dilate0_erode1_8con(p->byt, s0, s1, 1);
  if(initdetectionname)
    arraytofitsimg(initdetectionname, "Eroded", BYTE_IMG, p->byt,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(verb)
    {
      sprintf(report, "Eroded %lu times (%s connectivity).", p->numerosion,
              p->erodengb==4 ? "4" : "8");
      reporttiming(&t1, report, 1);
    }



  /* Do the opening: */
  if(verb) gettimeofday(&t1, NULL);
  opening(p->byt, s0, s1, p->opening, p->openingngb);
  if(initdetectionname)
    arraytofitsimg(initdetectionname, "Opened", BYTE_IMG, p->byt,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(verb)
    {
      sprintf(report, "Opened (depth: %lu, %s connectivity).",
              p->opening, p->openingngb==4 ? "4" : "8");
      reporttiming(&t1, report, 1);
    }



  /* Label the connected regions: */
  if(verb) gettimeofday(&t1, NULL);
  numlabs=BF_concmp(p->byt, p->olab, s0, s1, 4);
  if(initdetectionname)
    arraytofitsimg(initdetectionname, "Labeled", LONG_IMG, p->olab,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);
  if(verb)
    {
      sprintf(report, "%lu initial detections found.", numlabs-1);
      reporttiming(&t1, report, 1);
    }



  /* Remove the false detections */
  if(verb) gettimeofday(&t1, NULL);
  detectonmesh(p, &numlabs);
  if(verb)
    {
      sprintf(report, "Final number of detections: %lu.", numlabs-1);
      reporttiming(&t1, report, 1);
    }


  /* Clean up: */
  free(p->byt);
  free(p->conv);
  freemesh(smp);
  freemesh(lmp);
}
