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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/mesh.h>
#include <gnuastro/mode.h>
#include <gnuastro/qsort.h>
#include <gnuastro/array.h>
#include <gnuastro/timing.h>
#include <gnuastro/statistics.h>

#include "main.h"


void *
avestdonthread(void *inparam)
{
  struct gal_mesh_thread_params *mtp=(struct gal_mesh_thread_params *)inparam;
  struct gal_mesh_params *mp=mtp->mp;
  struct subtractskyparams *p=(struct subtractskyparams *)mp->params;

  float *mponeforall=mp->oneforall;
  float *oneforall=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  float *f, *cofa, *c=NULL;       /* cofa: convolved-oneforall */
  size_t modeindex=(size_t)(-1);
  size_t s0, s1, ind, row, num, start, is1=mp->s1;
  size_t i, *indexs=&mp->indexs[mtp->id*mp->thrdcols];
  float mirrordist=mp->mirrordist, minmodeq=mp->minmodeq;
  float ave, med, std, *sorted, sigclipmultip=p->sigclipmultip;
  float *imgend, *inimg=mp->img, *inconv=p->conv, modesym=0.0f;
  float *conv=NULL, *img, sigcliptolerance=p->sigcliptolerance;

  /* Allocate the oneforall array for the convolved image, since we
     only want to use those meshs whose convolved image mode is larger
     than minmodeq. */
  if(p->conv!=mp->img)
    {
      errno=0; cofa=malloc(mp->maxs0*mp->maxs1*sizeof *oneforall);
      if(cofa==NULL)
        error(EXIT_FAILURE, errno, "unable to allocate %lu bytes for"
              "cofa in avestdonthread of subtractsky.c",
              mp->maxs0*mp->maxs1*sizeof *mp->img);
    }
  else cofa=NULL;

  /* Start this thread's work: */
  for(i=0;indexs[i]!=GAL_THREADS_NON_THRD_INDEX;++i)
    {
      /* Prepare the values: */
      num=row=0;
      f=oneforall;
      ind=indexs[i];
      start=mp->start[ind];
      s0=mp->ts0[mp->types[ind]];
      s1=mp->ts1[mp->types[ind]];
      if(cofa) c=sorted=cofa; else sorted=oneforall;

      /* Copy all the non-NaN pixels images pixels of this mesh into
         the mesh array. Note that currently, the spatial positioning
         of the pixels is irrelevant, so we only keep those that are
         non-NaN. Recall that both the convolved an unconvolved image
         have the same NaN pixels.*/
      do
        {
          if(cofa) conv = inconv + start + row * is1;
          imgend=(img = inimg + start + row++ * is1 ) + s1;
          do
            {
              if(!isnan(*img))
                {
                  ++num;
                  *f++ = *img;
                  if(cofa) *c++=*conv;
                }
              if(cofa) ++conv;
            }
          while(++img<imgend);
        }
      while(row<s0);

      /* Do the desired operation on the mesh: */
      qsort(sorted, num, sizeof *oneforall, gal_qsort_float_increasing);
      gal_statistics_mode_index_in_sorted(sorted, num, mirrordist,
                                          &modeindex, &modesym);
      if( modesym>GAL_STATISTICS_MODE_SYM_GOOD
          && (float)modeindex/(float)num>minmodeq )
        {
          /* If cofa was defined, then oneforall was not sorted. */
          if(cofa)
            qsort(oneforall, num, sizeof *oneforall,
                  gal_qsort_float_increasing);

          /* Do sigma-clipping and save the result if it is
             accurate. Note that all meshs were initialized to NaN, so
             if they don't fit the criteria, they can simply be
             ignored. */
          if(gal_statistics_sigma_clip_converge(oneforall, 1, num,
                                                sigclipmultip, sigcliptolerance,
                                                &ave, &med, &std, 0))
            {
              mp->cgarray1[ind]=ave;
              if(mp->ngarrays==2) mp->cgarray2[ind]=std;
            }
        }
    }

  /* Free any allocated space and if multiple threads were used, wait
     until all other threads finish. */
  if(p->conv!=mp->img) free(cofa);
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}





void
subtractsky(struct subtractskyparams *p)
{
  struct gal_mesh_params *mp=&p->mp;

  long *meshindexs;
  struct timeval t1;
  int checkstd=p->checkstd;
  size_t s0=mp->s0, s1=mp->s1;
  float *sky, *std, *skysubtracted;


  /* Prepare the mesh array. */
  gettimeofday(&t1, NULL);
  gal_mesh_make_mesh(mp);
  if(p->meshname)
    {
      gal_mesh_check_mesh_id(mp, &meshindexs);
      gal_fits_array_to_file(p->meshname, "Input", FLOAT_IMG,
                             p->mp.img, s0, s1, p->anyblank, p->wcs,
                             NULL, SPACK_STRING);
      gal_fits_array_to_file(p->meshname, "MeshIndexs", LONG_IMG,
                             meshindexs, s0, s1, 0, p->wcs,
                             NULL, SPACK_STRING);
      free(meshindexs);
    }
  if(p->cp.verb) gal_timing_report(&t1, "Mesh grid ready.", 1);



  /* Convolve the image if the user has asked for it: */
  if(p->up.kernelnameset)
    {
      gal_mesh_spatial_convolve_on_mesh(mp, &p->conv);
      if(p->convname)
        {
          gal_fits_array_to_file(p->convname, "Input", FLOAT_IMG,
                                 p->mp.img, s0, s1, p->anyblank,
                                 p->wcs, NULL, SPACK_STRING);
          gal_fits_array_to_file(p->convname, "Input", FLOAT_IMG,
                                 p->conv, s0, s1, p->anyblank, p->wcs,
                                 NULL, SPACK_STRING);
        }
    }
  else p->conv=p->mp.img;
  if(p->cp.verb)
    gal_timing_report(&t1, "Input image convolved with kernel.", 1);



  /* Find the sky value and its standard deviation on each mesh. */
  gal_mesh_operate_on_mesh(mp, avestdonthread, sizeof(float), checkstd, 1);
  if(p->skyname)
    gal_mesh_value_file(mp, p->skyname, "Sky value", "Sky STD", p->wcs,
                        SPACK_STRING);
  if(p->cp.verb)
    gal_timing_report(&t1, "Sky and its STD found on some meshes.", 1);



  /* Interpolate over the meshs to fill all the blank ones in both the
     sky and the standard deviation arrays: */
  gal_mesh_interpolate(mp, "Interpolating the sky and its standard deviation");
  if(p->skyname)
    gal_mesh_value_file(mp, p->skyname, "Sky Interpolated",
                  "Sky STD interpolated", p->wcs, SPACK_STRING);
  if(p->cp.verb)
    gal_timing_report(&t1, "All blank meshs filled (interplated).", 1);



  /* Smooth the interpolated array:  */
  if(mp->smoothwidth>1)
    {
      gal_mesh_smooth(mp);
      if(p->cp.verb)
        gal_timing_report(&t1, "Mesh grid smoothed.", 1);
    }


  /* Make the sky array and save it if the user has asked for it: */
  gal_mesh_check_garray(mp, &sky, &std);
  if(p->skyname)
    gal_mesh_value_file(mp, p->skyname, "Sky Smoothed", "Sky STD smoothed",
                        p->wcs, SPACK_STRING);


  /* Subtract the sky value */
  gal_array_fmultip_const(sky, s0*s1, -1.0f);
  skysubtracted=gal_array_fsum_arrays(mp->img, sky, s0*s1);
  gal_fits_array_to_file(p->cp.output ,"SkySubtracted", FLOAT_IMG,
                         skysubtracted, s0, s1, p->anyblank, p->wcs,
                         NULL, SPACK_STRING);


  /* Clean up: */
  free(sky);
  gal_mesh_free_mesh(mp);
  free(skysubtracted);
  if(checkstd) free(std);
  if(p->up.kernelnameset) free(p->conv);
}
