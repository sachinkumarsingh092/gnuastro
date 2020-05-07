/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef ONECROP_H
#define ONECROP_H

#include <fitsio.h>

#include <gnuastro/threads.h>

struct onecropparams
{
  void *array;

  /* Pointer to basic structure: */
  struct   cropparams *p;

  /* About input image. */
  size_t              in_ind;  /* Index of this image in the input names.  */
  fitsfile           *infits;  /* Pointer to the input FITS image.         */
  long        fpixel[MAXDIM];  /* Position of first pixel in input image.  */
  long        lpixel[MAXDIM];  /* Position of last pixel in input image.   */
  double           *ipolygon;  /* Input image based polygon vertices.      */

  /* Output (cropped) image. */
  size_t             out_ind;  /* Index of this crop in the output list.   */
  double       world[MAXDIM];  /* World coordinates of crop center.        */
  double       sized[MAXDIM];  /* Width and height of image in degrees.    */
  double         corners[24];  /* RA and Dec of this crop's corners.       */
  double      equatorcorr[2];  /* Crop crosses the equator, see wcsmode.c. */
  fitsfile          *outfits;  /* Pointer to the output FITS image.        */

  /* For log */
  char                 *name;  /* Filename of crop.                        */
  size_t              numimg;  /* Number of images used to make this crop. */
  unsigned char centerfilled;  /* ==1 if the center is filled.             */

  /* Thread parameters. */
  size_t             *indexs;  /* Indexs to be used in this thread.        */
  pthread_barrier_t       *b;  /* pthread barrier to keep threads waiting. */
};

void
onecrop_name(struct onecropparams *crp);

void
onecrop(struct onecropparams *crp);

int
onecrop_center_filled(struct onecropparams *crp);

void
crop_print_log(struct onecropparams *p);

#endif
