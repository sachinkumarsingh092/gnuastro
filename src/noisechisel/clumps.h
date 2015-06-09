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
#ifndef CLUMPS_H
#define CLUMPS_H


struct clumpsthreadparams
{
  /* Main NoiseChisel parameters: */
  struct noisechiselparams *p;  /* Main NoiseChisel structure. */

  /* Threads (only for detections, the noise uses the mesh threads.) */
  size_t       thislabel; /* The initial label of this detection.        */
  size_t              id; /* ID of this thread.                          */
  size_t       *allareas; /* Array keeping the areas of all detections.  */
  size_t    **alllabinds; /* Array of pointers to the indexs of all dets.*/
  size_t         *indexs; /* 2D array of initial indexs for each thread. */
  pthread_barrier_t   *b; /* pthreads barrier for running threads.       */

  /* Box coordinates for this thread: */
  size_t              x0; /* Bottom left corner on x axis.               */
  size_t              y0; /* Bottom left corner on y axis.               */
  size_t              x1; /* Top right corner on x axis.                 */
  size_t              y1; /* Top right corner on y axis.                 */
  double            *xys; /* The light weighted center of each clump.    */

  /* Other basic parameters: */
  float              std; /* Standard deviation on this detection.       */
  size_t        *topinds; /* Indexs of the top flux in each clump.       */
  size_t       numclumps; /* Number of clumps in this set of pixels.     */
  size_t      numobjects; /* Number of objects in this detected region.  */
  size_t            area; /* Array keeping the areas of all detections.  */
  size_t           *inds; /* Array of pointers to the indexs of all dets.*/
  size_t      *blankinds; /* Array of pixels which should be grown.      */
  size_t       numblanks; /* Number of blank pixels.                     */
  long     *segtoobjlabs; /* Convert from grown segments to object label.*/
};


/* Important sizes and values (do not change). */
#define SEGMENTNOOBJ     0
#define SEGMENTMASKED   -4
#define SEGMENTTMPCHECK -3
#define SEGMENTINIT     -2
#define SEGMENTRIVER    -1
#define INFOTABCOLS      5
#define WNGBSIZE        20
#define NOTOPIND        (size_t)(-1)

void
oversegment(struct clumpsthreadparams *ctp);

void
growclumps(struct clumpsthreadparams *ctp, int withrivers);

void
clumpsntable(struct clumpsthreadparams *ctp, float **sntable);

void
clumpsngrid(struct noisechiselparams *p);

void
removefalseclumps(struct clumpsthreadparams *ctp, float *sntable);

#endif
