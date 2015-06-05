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

#include "fitsarrayvv.h"

#include "main.h"

#include "segmentation.h"




/******************************************************************/
/*****************         Main function        *******************/
/******************************************************************/
void
segmentation(struct noisechiselparams *p)
{
  size_t s0=p->smp.s0, s1=p->smp.s1;
  char *segmentationname=p->segmentationname;

  /* Start the steps image: */
  if(segmentationname)
    {
      arraytofitsimg(segmentationname, "Input-SkySubtracted",
                     FLOAT_IMG, p->img, s0, s1, p->numblank,
                     p->wcs, NULL, SPACK_STRING);
      arraytofitsimg(segmentationname, "Convolved-SkySubtracted",
                     FLOAT_IMG, p->conv, s0, s1, p->numblank,
                     p->wcs, NULL, SPACK_STRING);
    }
}
