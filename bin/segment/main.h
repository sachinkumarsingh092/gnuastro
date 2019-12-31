/*********************************************************************
Segment - Segment initial labels based on signal structure.
Segment is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2018-2020, Free Software Foundation, Inc.

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
#ifndef MAIN_H
#define MAIN_H

/* Include necessary headers */
#include <gnuastro/data.h>

#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "Segment"    /* Program full name.       */
#define PROGRAM_EXEC   "astsegment" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION


/* Macros */
#define DETECTION_ALL  "all"






/* Main program parameters structure */
struct segmentparams
{
  /* From command-line */
  struct gal_options_common_params  cp; /* Common parameters.             */
  struct gal_tile_two_layer_params ltl; /* Large tessellation.            */
  char             *inputname;  /* Input filename.                        */
  char            *kernelname;  /* Input kernel filename.                 */
  char                  *khdu;  /* Kernel HDU.                            */
  char         *convolvedname;  /* Convolved image (to avoid convolution).*/
  char                  *chdu;  /* HDU of convolved image.                */
  char         *detectionname;  /* Detection image file name.             */
  char                  *dhdu;  /* Detection image file name.             */
  char               *skyname;  /* Filename of Sky image.                 */
  char                *skyhdu;  /* Filename of Sky image.                 */
  char               *stdname;  /* File name of Standard deviation image. */
  char                *stdhdu;  /* HDU of Stanard deviation image.        */
  uint8_t            variance;  /* The input STD is actually variance.    */
  uint8_t           rawoutput;  /* Output only object and clump labels.   */

  float            minskyfrac;  /* Undetected area min. frac. in tile.    */
  uint8_t              minima;  /* Build clumps from their minima, maxima.*/
  size_t            snminarea;  /* Minimum area for segmentation.         */
  uint8_t             checksn;  /* Save the clump S/N values to a file.   */
  size_t          minnumfalse;  /* Min No. of det/seg for true quantile.  */
  float               snquant;  /* Quantile of clumps in sky for true S/N.*/
  uint8_t    keepmaxnearriver;  /* Keep clumps with a peak near a river.  */
  float         clumpsnthresh;  /* Clump S/N threshold.                   */
  uint8_t          onlyclumps;  /* Finish after finding true clumps.      */
  float               gthresh;  /* Multiple of STD to stop growing clumps.*/
  size_t       minriverlength;  /* Min, len of good grown clump rivers.   */
  float           objbordersn;  /* Minimum S/N for grown clumps to be one.*/
  uint8_t         grownclumps;  /* Save grown clumps instead of original. */
  uint8_t  continueaftercheck;  /* Don't abort after the check steps.     */
  uint8_t   checksegmentation;  /* Save the segmentation steps in file.   */

  /* Internal. */
  char        *clumpsn_s_name;  /* Sky clump S/N name.                    */
  char        *clumpsn_d_name;  /* Detection clumps S/N name.             */
  char      *segmentationname;  /* Name of segmentation steps file.       */

  gal_data_t           *input;  /* Input dataset.                         */
  gal_data_t          *kernel;  /* Given kernel for convolution.          */
  gal_data_t            *conv;  /* Convolved dataset.                     */
  gal_data_t          *binary;  /* For binary operations.                 */
  gal_data_t          *olabel;  /* Object labels.                         */
  gal_data_t          *clabel;  /* Clumps labels.                         */
  gal_data_t             *std;  /* STD of undetected pixels, per tile.    */
  gal_data_t       *clumpvals;  /* Values to build clumps (avoid bugs).   */

  float               cpscorr;  /* Counts/second correction.              */
  size_t        numdetections;  /* Number of final detections.            */
  size_t            numclumps;  /* Number of true clumps.                 */
  size_t           numobjects;  /* Number of objects.                     */

  char     *useddetectionname;  /* Name of file USED for detection image. */
  char           *usedstdname;  /* Name of file USED for sky STD image.   */

  float                medstd;  /* For output STD image: median STD.      */
  float                minstd;  /* For output STD image: median STD.      */
  float                maxstd;  /* For output STD image: median STD.      */

  /* Output: */
  time_t              rawtime;  /* Starting time of the program.          */
};

#endif
