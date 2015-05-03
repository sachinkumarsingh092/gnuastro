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
  struct timeval t1;
  struct meshparams *mp=&p->mp;

  /* Prepare the mesh array. */
  gettimeofday(&t1, NULL);
  p->mp.numthreads=p->cp.numthreads;
  makemesh(mp);
  if(p->meshname)
    {
      checkmesh(p->img, mp, &meshindexs);
      arraytofitsimg(p->meshname, "Input", FLOAT_IMG, p->img,
                     mp->s0, mp->s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      arraytofitsimg(p->meshname, "MeshIndexs", LONG_IMG, meshindexs,
                     mp->s0, mp->s1, p->numblank, p->wcs, NULL,
                     SPACK_STRING);
      free(meshindexs);
    }
  if(p->cp.verb) reporttiming(&t1, "Mesh grid ready", 1);

}
