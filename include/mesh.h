/*********************************************************************
meshgrid -- Create a mesh grid ready for multithreaded analysis.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef MESHGRID_H
#define MESHGRID_H

struct meshparams
{
  /* Image: */
  size_t           s0;	 /* Height of input image.                      */
  size_t           s1;	 /* Width of input image.                       */

  /* Threads: */
  size_t   numthreads;   /* Number of CPU threads.                      */
  size_t      *indexs;   /* 2D array of mesh indexs for each thread.    */
  size_t     thrdcols;   /* The number of columns in indexs array.      */

  /* Channel(s): */
  size_t          nch;   /* Total number of channels.                   */
  size_t         nch1;   /* Number of channels along first FITS axis.   */
  size_t         nch2;   /* Number of channels along first FITS axis.   */

  /* Meshs: */
  float  lastmeshfrac;   /* Last mesh fraction of size, to add new.     */
  size_t     meshsize;   /* Size of each mesh.                          */
  size_t       nmeshc;   /* Number of meshes in each channel.           */
  size_t       nmeshi;   /* Number of meshes in all image.              */
  size_t       *start;   /* Starting pixel for each mesh.               */
  size_t       *types;   /* The type of each mesh.                      */
  size_t     *chindex;   /* The index of each mesh in its channel.      */
  size_t    *imgindex;   /* The index of each mesh in the image.        */
  size_t          gs0;   /* Number of meshes on axis 0 in each channel. */
  size_t          gs1;   /* Number of meshes on axis 1 in each channel. */

  /* Mesh types and information: */
  size_t       ts0[4];   /* Size (along first FITS axis) of mesh types. */
  size_t       ts1[4];   /* Size (along second FITS axis) of mesh types.*/
};

void
makemesh(struct meshparams *mp);

void
checkmesh(float *img, struct meshparams *mp, long **out);

void
freemesh(struct meshparams *mp);

#endif
