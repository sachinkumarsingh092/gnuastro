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
  size_t i, s0=smp->s0, s1=smp->s1;


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
  if(p->detectionname)
    {
      arraytofitsimg(p->detectionname, "Input", FLOAT_IMG, smp->img,
                     s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(p->detectionname, "Convolved", FLOAT_IMG, p->conv,
                     s0, s1, p->numblank, p->wcs, NULL, SPACK_STRING);
    }
  if(verb) reporttiming(&t1, "Convolved with kernel.", 1);



  /* Do the initial detection: */
  if(verb)
    {
      reporttiming(NULL, "Starting to find initial detections.", 1);
      gettimeofday(&t1, NULL);
    }
  initialdetection(p);
  if(verb)
    {
      sprintf(report, "%lu initial detections found.", p->numobjects-1);
      reporttiming(&t1, report, 1);
    }


  /* Remove the false detections */
  if(verb)
    {
      reporttiming(NULL, "Starting to find and remove false detections.", 1);
      gettimeofday(&t1, NULL);
    }
  onlytruedetections(p);
  if(verb)
    {
      sprintf(report, "%lu true detections identified.", p->numobjects-1);
      reporttiming(&t1, report, 1);
    }


  /* Dilate the byt array and find the new number of detections: */
  if(verb) gettimeofday(&t1, NULL);
  if(p->dilate)
    {
      for(i=0;i<p->dilate;++i)
        dilate0_erode1_8con(p->byt, s0, s1, 0);
        p->numobjects=BF_concmp(p->byt, p->olab, s0, s1, 4);
      if(verb)
        {
          sprintf(report, "%lu detections after %lu dilation%s",
                  p->numobjects-1, p->dilate, p->dilate>1 ? "s." : ".");
          reporttiming(&t1, report, 1);
        }
    }
  if(p->detectionname)
    arraytofitsimg(p->detectionname, "Dilated", LONG_IMG, p->olab,
                   s0, s1, 0, p->wcs, NULL, SPACK_STRING);

  /* Clean up: */
  free(p->conv);
  freemesh(smp);
  freemesh(lmp);
}
