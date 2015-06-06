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
#ifndef MAIN_H
#define MAIN_H

#include "mesh.h"
#include "fitsarrayvv.h"
#include "commonparams.h"

/* Progarm name macros: */
#define SPACK_VERSION   "0.1"
#define SPACK           "astnoisechisel" /* Subpackage executable name. */
#define SPACK_NAME      "NoiseChisel"    /* Subpackage full name.       */
#define SPACK_STRING    SPACK_NAME" ("PACKAGE_STRING") "SPACK_VERSION




struct uiparams
{
  char          *inputname;  /* Name of input file.                 */
  char           *maskname;  /* Name of mask image file.            */
  char               *mhdu;  /* Name of mask image header name.     */
  char         *kernelname;  /* Name of kernel image file.          */
  char               *khdu;  /* Name of kernel header name.         */

  int   fullconvolutionset;
  int fullinterpolationset;
  int        fullsmoothset;
  int          masknameset;
  int              mhduset;
  int        kernelnameset;
  int              khduset;
  int     skysubtractedset;

  int         smeshsizeset;
  int         lmeshsizeset;
  int              nch1set;
  int              nch2set;
  int      lastmeshfracset;
  int        numnearestset;
  int       smoothwidthset;
  int        mirrordistset;
  int          minmodeqset;

  int           qthreshset;
  int        numerosionset;
  int          erodengbset;
  int           openingset;
  int        openingngbset;
  int          minbfracset;
  int     sigclipmultipset;
  int  sigcliptoleranceset;
  int           dthreshset;
  int      detsnminareaset;
  int       minnumfalseset;
  int          detquantset;
  int    detsnhistnbinsset;
  int            dilateset;
};





struct noisechiselparams
{
  /* Other structures: */
  struct uiparams      up;  /* User interface parameters.                  */
  struct commonparams  cp;  /* Common parameters.                          */
  struct meshparams   smp;  /* Smaller mesh grid of input image.           */
  struct meshparams   lmp;  /* Larger mesh grid of input image.            */

  /* Input: */
  int                nwcs;  /* Number of WCS structures.                   */
  struct wcsprm      *wcs;  /* Pointer to WCS structures.                  */
  int              bitpix;  /* Input image bitpix value.                   */
  size_t         numblank;  /* Number of blank pixels in image.            */
  int       skysubtracted;  /* ==1: the input image is already sky subted. */

  /* output: */
  char          *meshname;  /* Name of --checkmesh output.                 */
  char        *threshname;  /* !=NULL: Name of threshold steps.            */
  char     *detectionname;  /* !=NULL: Name of detection steps.            */
  char  *detectionskyname;  /* !=NULL: Name of detection sky steps.        */
  char   *detectionsnname;  /* !=NULL: Name of detection S/N on meshs.     */
  char           *skyname;  /* !=NULL: Name of image showing sky and STD.  */
  char  *segmentationname;  /* !=NULL: Name of segmentation steps.         */

  /* Detection: */
  float             *conv;  /* Convolved image.                            */
  float           qthresh;  /* Quantile threshold on convolved img.        */
  size_t       numerosion;  /* Number of times to erode thresholded image. */
  int            erodengb;  /* Use 4 or 8 connectivity in erosion.         */
  size_t          opening;  /* Depth of opening to apply to eroded image.  */
  int          openingngb;  /* Use 4 or 8 connectivity in opening.         */
  float          minbfrac;  /* Minimum fraction of undetected pixels.      */
  float     sigclipmultip;  /* Multiple of standard deviation, sigma clip. */
  float  sigcliptolerance;  /* Tolerance in sigma clip.                    */
  float           dthresh;  /* Threshold to remove false detections.       */
  size_t     detsnminarea;  /* Minimum area to calculate S/N for detection.*/
  size_t      minnumfalse;  /* Min No. false detections/segments for quant.*/
  float          detquant;  /* False detection S/N quantile to find true.  */
  size_t   detsnhistnbins;  /* ==1: Save false detection S/Ns histogram.   */
  size_t           dilate;  /* Number of times to dilate true detections.  */

  /* Operating mode: */

  /* Internal: */
  time_t          rawtime;  /* Starting time of the program.               */
  long              *olab;  /* Object labels.                              */
  long              *clab;  /* Clump labels.                               */
  unsigned char      *byt;  /* Array of single bytes for binary operations.*/
  unsigned char     *dbyt;  /* False detection removal thresholded array.  */
  float           cpscorr;  /* Correction for counts per second data.      */
  unsigned char      b0f1;  /* ==1: we are now working on data, not noise. */
  float              *img;  /* Input image.                                */
  float              *std;  /* The standard deviation of each pixel.       */
  int             stepnum;  /* Number of step if user wants to see steps.  */
  size_t       numobjects;  /* Total number of objects (detections).       */
  size_t        numclumps;  /* Total number of clumps.                     */
  size_t         *topinds;  /* Indexs of the top flux pixels in each clump.*/
  size_t        relngb[8];  /* Indexs of neighbors.                        */
};

#endif
