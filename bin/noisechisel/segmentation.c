/*********************************************************************
NoiseChisel - Detect and segment signal in a noisy dataset.
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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/threads.h>

#include <gnuastro-internal/timing.h>

#include "main.h"

#include "ui.h"
#include "clumps.h"
#include "segmentation.h"





/***********************************************************************/
/*****************            Over detections          *****************/
/***********************************************************************/
/* Find the true clumps over each detection. */
static void *
segmentation_on_threads(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct clumps_params *clprm=(struct clumps_params *)(tprm->params);
  struct noisechiselparams *p=clprm->p;

  size_t i;
  gal_data_t *topinds;
  struct clumps_thread_params cltprm;

  /* Initialize the general parameters for this thread. */
  cltprm.clprm   = clprm;

  /* Go over all the detections given to this thread (counting from zero.) */
  for(i=0; tprm->indexs[i] != GAL_THREADS_NON_THRD_INDEX; ++i)
    {

      /* Set the ID of this detection, note that for the threads, we
         counted from zero, but the IDs start from 1, so we'll add a 1 to
         the ID given to this thread. */
      cltprm.id     = tprm->indexs[i]+1;
      cltprm.indexs = &clprm->labindexs[ cltprm.id ];


      /* The `topinds' array is only necessary when the user wants to
         ignore true clumps with a peak touching a river. */
      if(p->keepmaxnearriver==0)
        {
          /* Allocate the list of local maxima. For each clump there is
             going to be one local maxima. But we don't know the number of
             clumps a-priori, so we'll just allocate the number of pixels
             given to this detected region. */
          topinds=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1,
                                 cltprm.indexs->dsize, NULL, 0,
                                 p->cp.minmapsize, NULL, NULL, NULL);
          cltprm.topinds=topinds->array;
        }
      else { cltprm.topinds=NULL; topinds=NULL; }


      /* Find the clumps over this region. */
      clumps_oversegment(&cltprm);

      /* Make the clump S/N table. This table is made before (possibly)
         stopping the process (if a check is requested). This is because if
         the user has also asked for a check image, we can break out of the
         loop at that point.*/
      clumps_make_sn_table(&cltprm);

      /* If the user wanted to check the segmentation steps or the clump
         S/N values in a table, then we have to stop the process at this
         point. */
      if(clprm->step==1 || p->checkclumpsn)
        { gal_data_free(topinds); continue; }



      /* Clean up. */
      gal_data_free(topinds);
    }

  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* If the user wanted to see the S/N table in a file, this function will be
   called and will do the job. */
static void
segmentation_save_sn_table(struct clumps_params *clprm)
{
  char *msg;
  float *sarr;
  uint32_t *oiarr, *cioarr;
  size_t i, j, c=0, totclumps=0;
  gal_data_t *sn, *objind, *clumpinobj;
  struct noisechiselparams *p=clprm->p;
  struct gal_linkedlist_stll *comments=NULL;


  /* Find the total number of clumps in all the initial detections. Recall
     that the `size' values were one more than the actual number because
     the labelings start from 1. */
  for(i=1;i<p->numinitdets+1;++i) totclumps += clprm->sn[i].size-1;


  /* Allocate the columns for the table. */
  sn=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &totclumps, NULL, 0,
                    p->cp.minmapsize, "CLUMP_S/N", "ratio",
                    "Signal-to-noise ratio.");
  objind=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, &totclumps, NULL, 0,
                        p->cp.minmapsize, "HOST_DET_ID", "counter",
                        "ID of detection hosting this clump.");
  clumpinobj=gal_data_alloc(NULL, GAL_TYPE_UINT32, 1, &totclumps, NULL, 0,
                            p->cp.minmapsize, "CLUMP_ID_IN_OBJ", "counter",
                            "ID of clump in host detection.");


  /* Fill in the columns. */
  sarr=sn->array;
  oiarr=objind->array;
  cioarr=clumpinobj->array;
  for(i=1;i<p->numinitdets+1;++i)
    for(j=1;j<clprm->sn[i].size;++j)
      {
        oiarr[c]  = i;
        cioarr[c] = j;
        sarr[c]   = ((float *)(clprm->sn[i].array))[j];
        ++c;
      }


  /* Write the comments. */
  gal_linkedlist_add_to_stll(&comments, "See also: `CLUMPS_ALL_DET' HDU of "
                             "output with `--checksegmentation'.", 1);
  asprintf(&msg, "S/N values of `nan': clumps smaller than `--segsnminarea' "
           "of %zu.", p->segsnminarea);
  gal_linkedlist_add_to_stll(&comments, msg, 1);
  free(msg);
  gal_linkedlist_add_to_stll(&comments, "S/N of clumps over detected "
                             "regions.", 1);
  gal_table_comments_add_intro(&comments, PROGRAM_STRING, &p->rawtime);


  /* Set the column pointers and write them into a table.. */
  clumpinobj->next=sn;
  objind->next=clumpinobj;
  gal_table_write(objind, comments, p->cp.tableformat, p->clumpsn_d_name, 1);


  /* Clean up. */
  gal_data_free(sn);
  gal_data_free(objind);
  gal_data_free(clumpinobj);
  gal_linkedlist_free_stll(comments, 1);


  /* Abort NoiseChisel if necessary. */
  if(!p->continueaftercheck)
    ui_abort_after_check(p, p->clumpsn_s_name, p->clumpsn_d_name,
                         "showing all clump S/N values");
}





/* Find true clumps over the detected regions. */
void
segmentation_detections(struct noisechiselparams *p)
{
  struct clumps_params clprm;
  gal_data_t *labindexs, *claborig;


  /* Get the indexs of all the pixels in each label. */
  labindexs=clumps_label_indexs(p);


  /* Initialize the necessary thread parameters. Note that since the object
     labels begin from one, the `sn' array will have one extra element.*/
  clprm.p=p;
  clprm.sky0_det1=1;
  clprm.snind = NULL;
  clprm.labindexs=labindexs;
  clprm.sn=gal_data_array_calloc(p->numinitdets+1);


  /* Spin off the threads to start the work. Note that several steps are
     done on each tile within a thread. So if the user wants to check
     steps, we need to break out of the processing get an over-all output,
     then reset the input and call it again. So it will be slower, but its
     is natural, since the user is testing to find the correct combination
     of parameters for later use. */
  if(p->segmentationname)
    {
      /* Necessary initializations. */
      clprm.step=1;
      claborig=p->clabel;
      p->clabel=gal_data_copy(claborig);

      /* Do each step. */
      while(clprm.step<3)
        {
          /* Reset the temporary copy of clabel back to its original. */
          if(clprm.step>1)
            memcpy(p->clabel->array, claborig->array,
                   claborig->size*gal_type_sizeof(claborig->type));

          /* (Re-)do everything until this step. */
          gal_threads_spin_off(segmentation_on_threads, &clprm,
                               p->numinitdets, p->cp.numthreads);

          /* Set the extension name. */
          switch(clprm.step)
            {
            case 1: p->clabel->name = "CLUMPS_ALL_DET";      break;
            case 2: p->clabel->name = "TRUE_CLUMPS";         break;
            default:
              error(EXIT_FAILURE, 0, "a bug! the value %d is not recognized "
                    "in `segmentation_detections'. Please contact us at "
                    "%s so we can address the issue", clprm.step,
                    PACKAGE_BUGREPORT);
            }

          /* Write the temporary array into the check image. */
          gal_fits_img_write(p->clabel, p->segmentationname, NULL,
                             PROGRAM_STRING);

          /* If the user wanted to check the clump S/N values, then break
             out of the loop, we don't need the rest of the process any
             more. */
          if( clprm.step==1
              && ( p->checkclumpsn && !p->continueaftercheck ) ) break;

          /* Increment the step counter. */
          ++clprm.step;
        }

      /* Clean up (we don't need the original any more). */
      gal_data_free(claborig);
      p->clabel->name=NULL;
    }
  else
    {
      clprm.step=0;
      gal_threads_spin_off(segmentation_on_threads, &clprm, p->numinitdets,
                           p->cp.numthreads);
    }


  /* If the user wanted to see the S/N table, then make the S/N table. */
  if(p->checkclumpsn) segmentation_save_sn_table(&clprm);


  /* Free the array of indexs. */
  gal_data_array_free(clprm.sn, p->numinitdets, 1);
  gal_data_array_free(labindexs, p->numinitdets, 1);
}




















/***********************************************************************/
/*****************         High level function         *****************/
/***********************************************************************/
void
segmentation(struct noisechiselparams *p)
{
  float *f;
  uint32_t *l, *lf;

  /* To keep the user up to date. */
  if(!p->cp.quiet)
    gal_timing_report(NULL, "Starting over-segmentation (finding clumps).",
                      1);


  /* If a check segmentation image was requested, then put in the
     inputs. */
  if(p->segmentationname)
    {
      gal_fits_img_write(p->input, p->segmentationname, NULL, PROGRAM_STRING);
      gal_fits_img_write(p->conv, p->segmentationname, NULL, PROGRAM_STRING);
      gal_fits_img_write(p->olabel, p->segmentationname, NULL,
                         PROGRAM_STRING);
    }


  /* Allocate the clump labels image. */
  p->clabel=gal_data_alloc(NULL, p->olabel->type, p->olabel->ndim,
                           p->olabel->dsize, p->olabel->wcs, 1,
                           p->cp.minmapsize, NULL, NULL, NULL);


  /* Set any possibly existing NaN values to blank. */
  f=p->input->array; lf=(l=p->clabel->array)+p->clabel->size;
  do if(isnan(*f++)) *l=GAL_BLANK_UINT32; while(++l<lf);


  /* Find the clump S/N threshold. */
  clumps_true_find_sn_thresh(p);


  /* Reset the clabel array to find true clumps in objects. */
  f=p->input->array; lf=(l=p->clabel->array)+p->clabel->size;
  do *l = isnan(*f++) ? GAL_BLANK_UINT32 : 0; while(++l<lf);


  /* Find true clumps over the detected regions. */
  segmentation_detections(p);


  /* If the user wanted to check the segmentation and hasn't called
     `continueaftercheck', then stop NoiseChisel. */
  if(p->segmentationname && !p->continueaftercheck)
    ui_abort_after_check(p, p->segmentationname, NULL,
                         "showing all segmentation steps");
}
