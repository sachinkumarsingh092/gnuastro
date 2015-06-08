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
#include <string.h>
#include <stdlib.h>

#include "fitsarrayvv.h"

#include "main.h"

#include "label.h"
#include "clumps.h"
#include "segmentation.h"






/******************************************************************/
/*****************         Segmentation         *******************/
/******************************************************************/
void
segmentdetections(struct noisechiselparams *p, size_t numobjsinit,
                  size_t *labareas, size_t **labinds)
{

}




















/******************************************************************/
/*****************         Main function        *******************/
/******************************************************************/
void
segmentation(struct noisechiselparams *p)
{
  float *f, *fp;
  size_t i, s0=p->smp.s0, s1=p->smp.s1;
  char *segmentationname=p->segmentationname;
  size_t *labareas, **labinds, numobjsinit=p->numobjects;

  /* Start off the counter for the number of objects and clumps. The
     value to these variables will be the label that is given to the
     next clump or object found. Note that we stored a copy of the
     initial number of objects in the numobjsinit variable above.*/
  p->numclumps=1;
  p->numobjects=1;


  /* Start the steps image: */
  if(segmentationname)
    {
      arraytofitsimg(segmentationname, "Input-SkySubtracted",
                     FLOAT_IMG, p->img, s0, s1, p->numblank,
                     p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(segmentationname, "Convolved-SkySubtracted",
                     FLOAT_IMG, p->conv, s0, s1, p->numblank,
                     p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(segmentationname, "InitialLabels",
                     LONG_IMG, p->olab, s0, s1, 0, p->wcs,
                     NULL, SPACK_STRING);
    }

  /* All possible NaN pixels should be given the largest possible
     float flux in the convolved image (which is used for
     over-segmentation). NOTE that the convolved image is used for
     relative pixel values, not absolute ones. This is because NaN
     pixels might be in the centers of stars or bright objects (they
     might slice through a connected region). So we can't allow them
     to cut our objects. The safest approach is to start segmentation
     of each object or noise mesh with the NaNs it might contain. A
     NaN region will never be calculated in any flux measurement any
     way and if it is connecting two bright objects, they will be
     segmented because on the two sides of the NaN region, they have
     different fluxes.*/
  if(p->numblank)
    {
      fp=(f=p->conv)+s0*s1;
      do if(isnan(*f)) *f=FLT_MAX; while(++f<fp);
    }


  /* Find the true clump S/N threshold andput it in p->lmp.garray1. */
  p->b0f1=0;
  clumpsngrid(p);


  /* Save the indexs of all the detected labels. We will be working on
     each label independently: */
  labindexs(p->olab, s0*s1, numobjsinit, &labareas, &labinds);


  /* Set all the elements of p->clab to zero: */
  memset(p->clab, 0, s0*s1*sizeof *p->clab);


  /* Segment the detections: */
  segmentdetections(p, numobjsinit, labareas, labinds);

  /* Clean up */
  for(i=0;i<numobjsinit;++i) free(labinds[i]);
  free(labareas);
  free(labinds);
}
