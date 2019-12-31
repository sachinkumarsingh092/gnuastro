/*********************************************************************
NoiseChisel - Detect signal in a noisy dataset.
NoiseChisel is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef MAIN_H
#define MAIN_H

/* Include necessary headers */
#include <gnuastro/data.h>

#include <gnuastro-internal/options.h>

/* Progarm names.  */
#define PROGRAM_NAME   "NoiseChisel"    /* Program full name.       */
#define PROGRAM_EXEC   "astnoisechisel" /* Program executable name. */
#define PROGRAM_STRING PROGRAM_NAME" (" PACKAGE_NAME ") " PACKAGE_VERSION





/* Main program parameters structure */
struct noisechiselparams
{
  /* From command-line */
  struct gal_options_common_params cp; /* Common parameters.              */
  struct gal_tile_two_layer_params ltl;/* Large tessellation.             */
  char             *inputname;  /* Input filename.                        */
  char            *kernelname;  /* Input kernel filename.                 */
  char                  *khdu;  /* Kernel HDU.                            */
  char         *convolvedname;  /* Convolved image (to avoid convolution).*/
  char                  *chdu;  /* HDU of convolved image.                */
  char        *widekernelname;  /* Name of wider kernel to be used.       */
  char                  *whdu;  /* Wide kernel HDU.                       */

  uint8_t  continueaftercheck;  /* Don't abort after the check steps.     */
  uint8_t  ignoreblankintiles;  /* Ignore input's blank values.           */
  uint8_t           rawoutput;  /* Only detection & 1 elem/tile output.   */
  uint8_t               label;  /* Label detections that are connected.   */

  float          meanmedqdiff;  /* Difference between mode and median.    */
  float               qthresh;  /* Quantile threshold on convolved image. */
  float          outliersigma;  /* Multiple of sigma to define outlier.   */
  double      outliersclip[2];  /* Outlier Sigma-clipping params.         */
  size_t          smoothwidth;  /* Interpolation: flat kernel to smooth.  */
  uint8_t        checkqthresh;  /* Save the quantile threhsold steps.     */
  uint8_t   blankasforeground;  /* Blank as foreg. in erosion and opening.*/
  size_t                erode;  /* Number of erosions after thresholding. */
  size_t             erodengb;  /* Connectivity for erosion.              */
  float          noerodequant;  /* Quantile for no erosion.               */
  size_t              opening;  /* Depth of opening after erosion.        */
  size_t           openingngb;  /* Connectivity to use for opening.       */
  uint8_t      skyfracnoblank;  /* No blanks in estimating non-det frac.  */
  float            minskyfrac;  /* Undetected area min. frac. in tile.    */
  double         sigmaclip[2];  /* Sigma-clipping parameters.             */
  uint8_t         checkdetsky;  /* Check pseudo-detection sky value.      */
  float               dthresh;  /* Sigma threshold for Pseudo-detections. */
  size_t             dopening;  /* Depth of opening after dthresh.        */
  size_t          dopeningngb;  /* Connectivity for opening after dthresh.*/
  size_t              holengb;  /* Connectivity for defining a hole.      */
  size_t        pseudoconcomp;  /* Connectivity for connected components. */
  size_t            snminarea;  /* Minimum pseudo-detection area for S/N. */
  uint8_t             checksn;  /* Save pseudo-detection S/N values.      */
  size_t          minnumfalse;  /* Min No. of det/seg for true quantile.  */
  float               snquant;  /* True detection quantile.               */
  float              snthresh;  /* Manually input S/N value.              */
  float          detgrowquant;  /* Quantile to grow true detections.      */
  size_t   detgrowmaxholesize;  /* Max. size of holes to fill in growth.  */
  uint8_t       cleangrowndet;  /* Remove grown objects with small S/N.   */
  uint8_t      checkdetection;  /* Save all detection steps to a file.    */
  uint8_t            checksky;  /* Check the Sky value estimation.        */

  /* Internal. */
  char           *qthreshname;  /* Name of Quantile threshold check image.*/
  char            *detskyname;  /* Name of Initial det sky check image.   */
  char          *detsn_s_name;  /* Sky pseudo-detections S/N name.        */
  char          *detsn_d_name;  /* Detection pseudo-detections S/N name.  */
  char          *detsn_D_name;  /* Final detection S/N name.              */
  char         *detectionname;  /* Name of detection steps file.          */
  char               *skyname;  /* Name of Sky estimation steps file.     */

  gal_data_t           *input;  /* Input image.                           */
  gal_data_t          *kernel;  /* Sharper kernel.                        */
  gal_data_t      *widekernel;  /* Wider kernel.                          */
  gal_data_t            *conv;  /* Convolved wth sharper kernel.          */
  gal_data_t           *wconv;  /* Convolved with wider kernel.           */
  gal_data_t          *binary;  /* For binary operations.                 */
  gal_data_t          *olabel;  /* Labels of objects in the detection.    */
  gal_data_t   *expand_thresh;  /* Quantile threshold to expand per tile. */
  gal_data_t *exp_thresh_full;  /* Full array containing growth thresh.   */
  gal_data_t             *sky;  /* Mean of undetected pixels, per tile.   */
  gal_data_t             *std;  /* STD of undetected pixels, per tile.    */
  size_t           maxtcontig;  /* Maximum contiguous space for a tile.   */
  size_t          maxltcontig;  /* Maximum contiguous space for a tile.   */
  size_t            *maxtsize;  /* Maximum size of a single small tile.   */
  size_t           *maxltsize;  /* Maximum size of a single large tile.   */
  size_t            numexpand;  /* Initial number of pixels to expand.    */
  time_t              rawtime;  /* Starting time of the program.          */

  float                medstd;  /* Median STD before interpolation.       */
  float                minstd;  /* Minimum STD before interpolation.      */
  float                maxstd;  /* Maximum STD before interpolation.      */
  float               cpscorr;  /* Counts/second correction.              */

  size_t       numinitialdets;  /* Number of initial detections.          */
  size_t        numdetections;  /* Number of final detections.            */
  float           detsnthresh;  /* Pseudo-detection S/N threshold.        */
};

#endif
