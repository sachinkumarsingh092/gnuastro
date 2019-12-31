/*********************************************************************
Statistics - Statistical analysis on input dataset.
Statistics is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2019-2020, Free Software Foundation, Inc.

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
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro-internal/checkset.h>

#include <gnuastro/wcs.h>
#include <gnuastro/binary.h>
#include <gnuastro/arithmetic.h>

#include "main.h"
#include "contour.h"

/* Pixels containing contour */
static gal_data_t *
contour_pixels(gal_data_t *input, double level, size_t minmapsize,
               int quietmmap)
{
  size_t one=1;
  uint8_t *b, *a, *af;
  int flags=GAL_ARITHMETIC_NUMOK;
  gal_data_t *number, *thresh, *eroded;

  /* Allocate the single-element dataset to use in arithmetic.*/
  number=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &one, NULL, 1,
                        -1, 1, NULL, NULL, NULL);
  *(double *)(number->array)=level;

  /* Only keep the pixels above the requested level, we are using the
     arithmetic library to not have to worry about the type of the
     input. */
  thresh=gal_arithmetic(GAL_ARITHMETIC_OP_GT, 1, flags, input, number);

  /* Erode the thresholded image by one. */
  eroded=gal_binary_erode(thresh, 1, 1, 0);

  /* Only keep the outer pixels. */
  b=eroded->array;
  af=(a=thresh->array)+thresh->size;
  do if(*b++ == 1) *a=0; while(++a<af);

  /* Clean up and return. */
  gal_data_free(number);
  gal_data_free(eroded);
  return thresh;
}





/* Given the indexs of the contours, write them in the proper format. */
static void
contour_pgfplots(gal_data_t *edgeindexs, gal_data_t *input, float level,
                 FILE *fp)
{
  double *xa, *ya;
  gal_data_t *x, *y, *tmp;
  size_t *s, *sf, w=input->dsize[1];

  /* Go through each connected edge and add the contour positions. */
  for(tmp=edgeindexs; tmp!=NULL; tmp=tmp->next)
    if(tmp->size>10)
      {
        if(input->wcs)
          {
            {
              /* Allocate the coordinate arrays. */
              x=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &tmp->size,
                               NULL, 0, -1, 1, NULL, NULL, NULL);
              y=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &tmp->size,
                               NULL, 0, -1, 1, NULL, NULL, NULL);

              /* Fill in the coordinates. */
              xa=x->array;
              ya=y->array;
              sf=(s=tmp->array)+tmp->size;
              do {*xa++=*s%w+1; *ya++=*s/w+1;} while(++s<sf);

              /* Convert the pixel positions to WCS. */
              x->next=y;
              gal_wcs_img_to_world(x, input->wcs, 1);

              /* Write them. */
              xa=x->array;
              ya=y->array;
              sf=(s=tmp->array)+tmp->size;
              do
                {fprintf(fp, "%.10f  %.10f  %f\n", *xa++, *ya++, level);}
              while(++s<sf);

              /* Clean up. */
              gal_data_free(x);
              gal_data_free(y);
            }
          }
        else
          {
            sf=(s=tmp->array)+tmp->size;
            do
              {fprintf(fp, "%zu  %zu  %f\n", *s%w+1, *s/w+1, level);}
            while(++s<sf);
          }

        /* To separate the different connected regions, we'll need to put
           an empty line between them. */
        fprintf(fp, "\n");
      }
}





/* Contour for each level. */
static void
contour_level(gal_data_t *input, double level, FILE *fp,
              size_t minmapsize, int quietmmap)
{
  gal_data_t *edge, *edgeindexs;

  /* Find the edge pixels given this threshold. */
  edge=contour_pixels(input, level, minmapsize, quietmmap);

  /* Indexs of the edges (separated by groups of connected edges). */
  edgeindexs=gal_binary_connected_indexs(edge, 2);

  /* Make the PGFPlots contours. */
  contour_pgfplots(edgeindexs, input, level, fp);

  /* Clean up and return. */
  gal_list_data_free(edgeindexs);
  gal_data_free(edge);
}





void
contour(struct statisticsparams *p)
{
  FILE *fp;
  char *outname;
  double *d, *df;
  uint8_t keepinputdir=p->cp.keepinputdir;

  /* Make sure the dataset is 2D. */
  if(p->input->ndim!=2)
    error(EXIT_FAILURE, 0, "contours are currently only supported for "
          "2D datasets (images). The input dataset has %zu dimensions",
          p->input->ndim);

  /* Make sure the output doesn't exist. */
  p->cp.keepinputdir = p->cp.output ? 1 : keepinputdir;
  outname=gal_checkset_automatic_output(&p->cp,
                                        ( p->cp.output
                                          ? p->cp.output
                                          : p->inputname ), "_contour.txt");
  p->cp.keepinputdir=keepinputdir;

  /* Open the output file. */
  fp=fopen(outname, "w+");

  /* Print basic information about the columns. */
  fprintf(fp, "# %s Contour positions\n", PACKAGE_STRING);
  fprintf(fp, "# Column 1: Coord_1 [position,f64] Position in first axis.\n");
  fprintf(fp, "# Column 2: Coord_2 [position,f64] Position in second axis.\n");
  fprintf(fp, "# Column 3: Level   [value,   f32] Contour level.\n");
  fprintf(fp, "#\n");
  fprintf(fp, "# Each connected contour is separated by an empty line.\n");
  fprintf(fp, "# This format is recognized in PGFPlots (package of LaTeX).\n");

  /* Estimate the contours for each level. */
  df=(d=p->contour->array)+p->contour->size;
  do
    contour_level(p->input, *d, fp, p->cp.minmapsize,
                  p->cp.quietmmap);
  while(++d<df);

  /* Clean up and cose the file. */
  free(outname);
  fclose(fp);
}
