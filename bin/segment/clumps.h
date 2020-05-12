/*********************************************************************
Segment - Segment initial labels based on signal structure.
Segment is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

Gnuastro is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

Gnuastro is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef CLUMPS_H
#define CLUMPS_H



/* Parameters for all threads. */
struct clumps_params
{
  /* General */
  int                     step; /* Counter if we want to check steps.      */
  int                sky0_det1; /* If working on the Sky or Detections.    */
  struct segmentparams      *p; /* Pointer to main Segment parameters.     */
  pthread_mutex_t     labmutex; /* Mutex to change the total numbers.      */

  /* For Sky region. */
  gal_data_t               *sn; /* Array of clump S/N tables.              */
  gal_data_t            *snind; /* Array of clump S/N index (for check).   */

  /* For detections. */
  gal_data_t        *labindexs; /* Array of 'gal_data_t' with obj indexs.  */
  size_t            totobjects; /* Total number of objects at any point.   */
  size_t             totclumps; /* Total number of clumps at any point.    */
};


/* Parameters for one thread (a tile or a detected region). */
struct clumps_thread_params
{
  size_t                    id; /* ID of this detection/tile over tile.    */
  size_t              *topinds; /* Indexs of all local maxima.             */
  size_t         numinitclumps; /* Number of clumps in tile/detection.     */
  size_t         numtrueclumps; /* Number of true clumps in tile/detection.*/
  size_t            numobjects; /* Number of objects over this clump.      */
  float                    std; /* Standard deviation of noise on center.  */
  gal_data_t           *indexs; /* Array containing indexs of this det.    */
  gal_data_t    *diffuseindexs; /* Diffuse region (after finding clumps).  */
  gal_data_t             *info; /* Information for all clumps.             */
  gal_data_t               *sn; /* Signal-to-noise ratio for these clumps. */
  gal_data_t            *snind; /* Index of S/N for these clumps.          */
  gal_data_t       *clumptoobj; /* Index of object that a clump belongs to.*/
  struct clumps_params  *clprm; /* Pointer to main structure.              */
};

void
clumps_grow_prepare_initial(struct clumps_thread_params *cltprm);

void
clumps_grow_prepare_final(struct clumps_thread_params *cltprm);

void
clumps_grow(gal_data_t *labels, gal_data_t *diffuseindexs, int withrivers,
            int connectivity);

void
clumps_true_find_sn_thresh(struct segmentparams *p);

void
clumps_make_sn_table(struct clumps_thread_params *cltprm);

gal_data_t *
clumps_det_label_indexs(struct segmentparams *p);

void
clumps_det_keep_true_relabel(struct clumps_thread_params *cltprm);

#endif
