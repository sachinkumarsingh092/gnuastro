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
#ifndef SEGMENTATION_H
#define SEGMENTATION_H


struct segmentationparams
{
  /* Main NoiseChisel parameters: */
  struct noisechiselparams *p;  /* Main NoiseChisel structure. */

  /* Indexs and areas of all the labels. */
  size_t     thisinitlab; /* Initial label of this whole detection.      */
  size_t       *labareas; /* Array keeping the areas of all detections.  */
  size_t       **labinds; /* Array of pointers to the indexs of all dets.*/

  /* Threads: */
  size_t              id; /* ID of this thread.                          */
  size_t         *indexs; /* 2D array of initial indexs for each thread. */
  pthread_barrier_t   *b; /* pthreads barrier for running threads.       */
};


void
segmentation(struct noisechiselparams *p);

#endif
