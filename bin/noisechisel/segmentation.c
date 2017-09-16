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
#include <gnuastro/binary.h>
#include <gnuastro/threads.h>
#include <gnuastro/dimension.h>

#include <gnuastro-internal/timing.h>

#include "main.h"

#include "ui.h"
#include "clumps.h"
#include "segmentation.h"




/***********************************************************************/
/*****************      Relabeling (grown) clumps      *****************/
/***********************************************************************/
/* Correct the label of an detection when it doesn't need segmentation (it
   is fully one object). The final labels for the object(s) with a detected
   region will be set later (don't forget that we have detections that are
   composed of multiple objects). So the labels within each detection start
   from 1.*/
static void
segmentation_relab_noseg(struct clumps_thread_params *cltprm)
{
  int32_t *olabel=cltprm->clprm->p->olabel->array;
  size_t *s=cltprm->indexs->array, *sf=s+cltprm->indexs->size;
  do olabel[ *s ] = 1; while(++s<sf);
}





/* Find the adjacency matrixs (number, sum and signal to noise) for the
   rivers between potentially separate objects in a detection region. They
   have to be allocated prior to entering this funciton.

   The way to find connected objects is through an adjacency matrix. It is
   a square matrix with a side equal to numobjs. So to see if regions `a`
   and `b` are connected. All we have to do is to look at element
   a*numobjs+b or b*numobjs+a and get the answer. Since the number of
   objects in a given region will not be too high, this is efficient. */
static void
segmentation_relab_to_objects(struct clumps_thread_params *cltprm)
{
  size_t amwidth=cltprm->numtrueclumps+1;
  struct noisechiselparams *p=cltprm->clprm->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  size_t mdsize[2]={amwidth, amwidth};
  gal_data_t *nums_d=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 2, mdsize, NULL,
                                    1, p->cp.minmapsize, NULL, NULL, NULL);
  gal_data_t *sums_d=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 2, mdsize, NULL,
                                    1, p->cp.minmapsize, NULL, NULL, NULL);
  gal_data_t *adjacency_d=gal_data_alloc(NULL, GAL_TYPE_UINT8, 2, mdsize,
                                         NULL, 1, p->cp.minmapsize, NULL,
                                         NULL, NULL);
  float *imgss=p->input->array;
  uint8_t *adjacency=adjacency_d->array;
  size_t nngb=gal_dimension_num_neighbors(ndim);
  int32_t *clumptoobj, *olabel=p->olabel->array;
  size_t *dinc=gal_dimension_increment(ndim, dsize);
  size_t *s, *sf, i, j, ii, rpnum, *nums=nums_d->array;
  double ave, rpsum, c=sqrt(1/p->cpscorr), *sums=sums_d->array;
  double err=cltprm->std*cltprm->std*(p->skysubtracted?1.0f:2.0f);
  int32_t *ngblabs=gal_data_malloc_array(GAL_TYPE_UINT32, nngb, __func__,
                                         "ngblabs");


  /* Go over all the still-unlabeled pixels (if they exist) and see which
     labels they touch. In the process, get the average value of the
     river-pixel values and put them in the respective adjacency
     matrix. Note that at this point, the rivers are also part of the
     "diffuse" regions. So we don't need to go over all the indexs of this
     object, only its diffuse indexs. */
  if(cltprm->diffuseindexs->size)
    {
      sf=(s=cltprm->diffuseindexs->array)+cltprm->diffuseindexs->size;
      do
        /* We only want to work on pixels that have already been identified
           as touching more than one label: river pixels. */
        if( olabel[ *s ]==CLUMPS_RIVER )
          {
            /* Initialize the values. */
            i=ii=0;
            rpnum=1;              /* River-pixel number of points used. */
            rpsum=imgss[*s];      /* River-pixel sum of values used.    */
            memset(ngblabs, 0, nngb*sizeof *ngblabs);

            /* Check all the fully-connected neighbors of this pixel and
               see if it touches a label or not */
            GAL_DIMENSION_NEIGHBOR_OP(*s, ndim, dsize, ndim, dinc, {
                if( olabel[nind] > 0 )
                  {
                    /* Add this neighbor's value and increment the number. */
                    if( !isnan(imgss[nind]) ) { ++rpnum; rpsum+=imgss[nind]; }

                    /* Go over the already found neighbors and see if this
                       grown clump has already been considered or not. */
                    for(i=0;i<ii;++i) if(ngblabs[i]==olabel[nind]) break;

                    /* This is the first time we are getting to this
                       neighbor: */
                    if(i==ii) ngblabs[ ii++ ] = olabel[nind];
                  }
              } );

            /* For a check:
            if(ii>0)
              {
                printf("%zu, %zu:\n", *s%dsize[1]+1, *s/dsize[1]+1);
                for(i=0;i<ii;++i) printf("\t%u\n", ngblabs[i]);
              }
            */

            /* If more than one neighboring label was found, fill in the
               'sums' and 'nums' adjacency matrixs with the values for this
               pixel. Recall that ii is the number of neighboring labels to
               this river pixel. */
            if(ii>i)
              for(i=0;i<ii;++i)
                for(j=0;j<ii;++j)
                  if(i!=j)
                    {
                      /* For safety, we will fill both sides of the
                         diagonal. */
                      ++nums[ ngblabs[i] * amwidth + ngblabs[j] ];
                      ++nums[ ngblabs[j] * amwidth + ngblabs[i] ];
                      sums[ ngblabs[i] * amwidth + ngblabs[j] ] +=
                        rpsum/rpnum;
                      sums[ ngblabs[j] * amwidth + ngblabs[i] ] +=
                        rpsum/rpnum;
                    }
          }
      while(++s<sf);


      /* We now have the average values and number of all rivers between
         the grown clumps. We now want to finalize their connection (given
         the user's criteria). */
      for(i=1;i<amwidth;++i)
        for(j=1;j<i;++j)
          {
            ii = i * amwidth + j;
            if(nums[ii]>p->minriverlength)       /* There is a connection. */
              {
                /* For easy reading. */
                ave=sums[ii]/nums[ii];

                /* In case the average is negative (only possible if `sums'
                   is negative), don't change the adjacency: it is already
                   initialized to zero. Note that even an area of 1 is
                   acceptable, and we put no area criteria here, because
                   the fact that a river exists between two clumps is
                   important. */
                if( ave>0.0f && ( c * ave / sqrt(ave+err) ) > p->objbordersn )
                  {
                    adjacency[ii]=1;   /* We want to set both sides of the */
                    adjacency[ j * amwidth + i ] = 1; /* Symmetric matrix. */
                  }
              }
          }


      /* For a check:
      if(cltprm->id==XXX)
        {
          printf("=====================\n");
          printf("%zu:\n--------\n", cltprm->id);
          for(i=1;i<amwidth;++i)
            {
              printf(" %zu...\n", i);
              for(j=1;j<amwidth;++j)
                {
                  ii=i*amwidth+j;
                  if(nums[ii])
                    {
                      ave=sums[ii]/nums[ii];
                      printf("    ...%zu: N:%-4zu S:%-10.2f S/N: %-10.2f "
                             "--> %u\n", j, nums[ii], sums[ii],
                             c*ave/sqrt(ave+err), adjacency[ii]);
                    }
                }
              printf("\n");
            }
        }
      */


      /* Calculate the new labels for each grown clump. */
      cltprm->clumptoobj = gal_binary_connected_adjacency_matrix(adjacency_d,
                                                         &cltprm->numobjects);
      clumptoobj = cltprm->clumptoobj->array;
    }

  /* There was no list of diffuse pixels, this happens when the user sets a
     very high `gthresh' threshold and wants to make sure that each clump
     is a separate object. So we need to define the number of objects and
     `clumptoobj' manually. */
  else
    {
      /* Allocate the `clumptoobj' array. */
      cltprm->clumptoobj = gal_data_alloc(NULL, GAL_TYPE_INT32, 1, &amwidth,
                                          NULL, 1, p->cp.minmapsize, NULL,
                                          NULL, NULL);
      clumptoobj = cltprm->clumptoobj->array;

      /* Fill in the `clumptoobj' array with the indexs of the objects. */
      for(i=0;i<amwidth;++i) clumptoobj[i]=i;

      /* Set the number of objects. */
      cltprm->numobjects = cltprm->numtrueclumps;
    }


  /* For a check
  if(cltprm->id==XXXX)
    {
      printf("NUMTRUECLUMPS: %zu\n----------\n", cltprm->numtrueclumps);
      for(i=0;i<cltprm->numtrueclumps+1;++i)
        printf("\t%zu --> %d\n", i, clumptoobj[i]);
      printf("=== numobjects: %zu====\n", cltprm->numobjects);
      exit(0);
    }
  */


  /* Correct all the labels. */
  sf=(s=cltprm->indexs->array)+cltprm->indexs->size;
  do
    if( olabel[*s] > 0 )
      olabel[*s] = clumptoobj[ olabel[*s] ];
  while(++s<sf);


  /* Clean up and return. */
  free(dinc);
  free(ngblabs);
  gal_data_free(nums_d);
  gal_data_free(sums_d);
  gal_data_free(adjacency_d);
}




/* The correspondance between the clumps and objects has been found. With
   this function, we want to correct the clump labels so the clump IDs in
   each object start from 1 and are contiguous. */
static void
segmentation_relab_clumps_in_objects(struct clumps_thread_params *cltprm)
{
  size_t numobjects=cltprm->numobjects, numtrueclumps=cltprm->numtrueclumps;

  int32_t *clumptoobj=cltprm->clumptoobj->array;
  int32_t *clabel=cltprm->clprm->p->clabel->array;
  size_t i, *s=cltprm->indexs->array, *sf=s+cltprm->indexs->size;
  size_t *nclumpsinobj=gal_data_calloc_array(GAL_TYPE_SIZE_T, numobjects+1,
                                             __func__, "nclumpsinobj");
  int32_t *newlabs=gal_data_calloc_array(GAL_TYPE_UINT32, numtrueclumps+1,
                                         __func__, "newlabs");

  /* Fill both arrays. */
  for(i=1;i<numtrueclumps+1;++i)
    newlabs[i] = ++nclumpsinobj[ clumptoobj[i] ];

  /* Reset the clump labels over the detection region. */
  do if(clabel[*s]>0) clabel[*s] = newlabs[ clabel[*s] ]; while(++s<sf);

  /* Clean up. */
  free(newlabs);
  free(nclumpsinobj);
}





/* Prior to this function, the objects have labels that are unique and
   contiguous (the labels are contiguous, not the objects!) within each
   detection and start from 1. However, for the final output, it is
   necessary that each object over the whole dataset have a unique
   ID. Since multiple threads are working on separate objects at every
   instance, this function will use a mutex to limit the reading and
   writing to the variable keeping the total number of objects counter. */
static void
segmentation_relab_overall(struct clumps_thread_params *cltprm)
{
  struct clumps_params *clprm=cltprm->clprm;
  int32_t startinglab, *olabel=clprm->p->olabel->array;
  size_t *s=cltprm->indexs->array, *sf=s+cltprm->indexs->size;

  /* Lock the mutex if we are working on more than one thread. NOTE: it is
     very important to keep the number of operations within the mutex to a
     minimum so other threads don't get delayed. */
  if(clprm->p->cp.numthreads>1)
    pthread_mutex_lock(&clprm->labmutex);

  /* Save the total number of objects found so far into `startinglab', then
     increment the total number of objects and clumps. */
  startinglab        = clprm->totobjects;
  clprm->totobjects += cltprm->numobjects;
  clprm->totclumps  += cltprm->numtrueclumps;

  /* Unlock the mutex (if it was locked). */
  if(clprm->p->cp.numthreads>1)
    pthread_mutex_unlock(&clprm->labmutex);

  /* Increase all the object labels by `startinglab'. */
  do olabel[*s] += startinglab; while(++s<sf);
}




















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

  size_t i, *s, *sf;
  gal_data_t *topinds;
  struct clumps_thread_params cltprm;
  int32_t *clabel=p->clabel->array, *olabel=p->olabel->array;

  /* Initialize the general parameters for this thread. */
  cltprm.clprm   = clprm;

  /* Go over all the detections given to this thread (counting from zero.) */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
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
         loop at that point.

         Note that the array of `gal_data_t' that keeps the S/N table for
         each detection is allocated before threading starts. However, when
         the user wants to inspect the steps, this function is called
         multiple times. So we need to avoid over-writing the allocations. */
      if( clprm->sn[ cltprm.id ].dsize==NULL )
        clumps_make_sn_table(&cltprm);
      else cltprm.sn=&clprm->sn[ cltprm.id ];


      /* If the user wanted to check the segmentation steps or the clump
         S/N values in a table, then we have to stop the process at this
         point. */
      if(clprm->step==1 || p->checkclumpsn)
        { gal_data_free(topinds); continue; }

      /* Only keep true clumps. */
      clumps_det_keep_true_relabel(&cltprm);
      gal_data_free(topinds);
      if(clprm->step==2) continue;


      /* Set the internal (with the detection) clump and object
         labels. Segmenting a detection into multiple objects is only
         defined when there is more than one true clump over the
         detection. When there is only one true clump
         (cltprm->numtrueclumps==1) or none (p->numtrueclumps==0), then
         just set the required preliminaries to make the next steps be
         generic for all cases. */
      if(cltprm.numtrueclumps<=1)
        {
          /* Set the basics. */
          cltprm.numobjects=1;
          segmentation_relab_noseg(&cltprm);


          /* If the user wanted a check image, this object doesn't
             change. */
          if( clprm->step >= 3 && clprm->step <= 6) continue;


          /* If the user has asked for grown clumps in the clumps image
             instead of the raw clumps, then replace the indexs in the
             `clabel' array is well. In this case, there will always be one
             "clump". */
          if(p->grownclumps)
            {
              sf=(s=cltprm.indexs->array)+cltprm.indexs->size;
              do clabel[ *s++ ] = 1; while(s<sf);
              cltprm.numtrueclumps=1;
            }
        }
      else
        {
          /* Grow the true clumps over the detection. */
          clumps_grow_prepare_initial(&cltprm);
          if(cltprm.diffuseindexs->size)
            clumps_grow(p->olabel, cltprm.diffuseindexs, 1);
          if(clprm->step==3)
            { gal_data_free(cltprm.diffuseindexs); continue; }


          /* If grown clumps are desired instead of the raw clumps, then
             replace all the grown clumps with those in clabel. */
          if(p->grownclumps)
            {
              sf=(s=cltprm.indexs->array)+cltprm.indexs->size;
              do
                if(olabel[*s]>0) clabel[*s]=olabel[*s];
              while(++s<sf);
            }


          /* Identify the objects in this detection using the grown clumps
             and correct the grown clump labels into new object labels. */
          segmentation_relab_to_objects(&cltprm);
          if(clprm->step==4)
            {
              gal_data_free(cltprm.clumptoobj);
              gal_data_free(cltprm.diffuseindexs);
              continue;
            }


          /* Continue the growth and cover the whole area, we don't need
             the diffuse indexs any more, so after filling the detected
             region, free the indexs. */
          if( cltprm.numobjects == 1 )
            segmentation_relab_noseg(&cltprm);
          else
            {
              /* Correct the labels so every non-labeled pixel can be
                 grown. */
              clumps_grow_prepare_final(&cltprm);

              /* Cover the whole area. */
              clumps_grow(p->olabel, cltprm.diffuseindexs, 0);
            }
          gal_data_free(cltprm.diffuseindexs);
          if(clprm->step==5) { gal_data_free(cltprm.clumptoobj); continue; }


          /* Correct the clump labels. Note that this is only necessary
             when there is more than object over the detection or when
             there were multiple clumps over the detection. */
          if(cltprm.numobjects>1)
            segmentation_relab_clumps_in_objects(&cltprm);
          gal_data_free(cltprm.clumptoobj);
          if(clprm->step==6) {continue;}
        }

      /* Convert the object labels to their final value */
      segmentation_relab_overall(&cltprm);
    }

  /* Wait until all the threads finish then return. */
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
  int32_t *oiarr, *cioarr;
  gal_list_str_t *comments=NULL;
  size_t i, j, c=0, totclumps=0;
  gal_data_t *sn, *objind, *clumpinobj;
  struct noisechiselparams *p=clprm->p;


  /* Find the total number of clumps in all the initial detections. Recall
     that the `size' values were one more than the actual number because
     the labelings start from 1. */
  for(i=1;i<p->numdetections+1;++i) totclumps += clprm->sn[i].size-1;


  /* Allocate the columns for the table. */
  sn=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &totclumps, NULL, 0,
                    p->cp.minmapsize, "CLUMP_S/N", "ratio",
                    "Signal-to-noise ratio.");
  objind=gal_data_alloc(NULL, GAL_TYPE_INT32, 1, &totclumps, NULL, 0,
                        p->cp.minmapsize, "HOST_DET_ID", "counter",
                        "ID of detection hosting this clump.");
  clumpinobj=gal_data_alloc(NULL, GAL_TYPE_INT32, 1, &totclumps, NULL, 0,
                            p->cp.minmapsize, "CLUMP_ID_IN_OBJ", "counter",
                            "ID of clump in host detection.");


  /* Fill in the columns. */
  sarr=sn->array;
  oiarr=objind->array;
  cioarr=clumpinobj->array;
  for(i=1;i<p->numdetections+1;++i)
    for(j=1;j<clprm->sn[i].size;++j)
      {
        oiarr[c]  = i;
        cioarr[c] = j;
        sarr[c]   = ((float *)(clprm->sn[i].array))[j];
        ++c;
      }


  /* Write the comments. */
  gal_list_str_add(&comments, "See also: `CLUMPS_ALL_DET' HDU of "
                   "output with `--checksegmentation'.", 1);
  asprintf(&msg, "S/N values of `nan': clumps smaller than `--segsnminarea' "
           "of %zu.", p->segsnminarea);
  gal_list_str_add(&comments, msg, 0);
  gal_list_str_add(&comments, "S/N of clumps over detected regions.", 1);
  gal_table_comments_add_intro(&comments, PROGRAM_STRING, &p->rawtime);


  /* Set the column pointers and write them into a table.. */
  clumpinobj->next=sn;
  objind->next=clumpinobj;
  gal_table_write(objind, comments, p->cp.tableformat, p->clumpsn_d_name, 1);


  /* Clean up. */
  gal_data_free(sn);
  gal_data_free(objind);
  gal_data_free(clumpinobj);
  gal_list_str_free(comments, 1);


  /* Abort NoiseChisel if necessary. */
  if(!p->continueaftercheck)
    ui_abort_after_check(p, p->clumpsn_s_name, p->clumpsn_d_name,
                         "showing all clump S/N values");
}





/* Find true clumps over the detected regions. */
static void
segmentation_detections(struct noisechiselparams *p)
{
  char *msg;
  struct clumps_params clprm;
  gal_data_t *labindexs, *claborig, *demo=NULL;


  /* Get the indexs of all the pixels in each label. */
  labindexs=clumps_det_label_indexs(p);


  /* Initialize the necessary thread parameters. Note that since the object
     labels begin from one, the `sn' array will have one extra element.*/
  clprm.p=p;
  clprm.sky0_det1=1;
  clprm.totclumps=0;
  clprm.totobjects=0;
  clprm.snind = NULL;
  clprm.labindexs=labindexs;
  clprm.sn=gal_data_array_calloc(p->numdetections+1);


  /* When more than one thread is to be used, initialize the mutex. */
  if( p->cp.numthreads > 1 ) pthread_mutex_init(&clprm.labmutex, NULL);


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
      while(clprm.step<8)
        {
          /* Reset the temporary copy of clabel back to its original. */
          if(clprm.step>1)
            memcpy(p->clabel->array, claborig->array,
                   claborig->size*gal_type_sizeof(claborig->type));

          /* (Re-)do everything until this step. */
          gal_threads_spin_off(segmentation_on_threads, &clprm,
                               p->numdetections, p->cp.numthreads);

          /* Set the extension name. */
          switch(clprm.step)
            {
            case 1:
              demo=p->clabel;
              demo->name = "DET_CLUMPS_ALL";
              if(!p->cp.quiet)
                {
                  asprintf(&msg, "Identified clumps over detections  "
                           "(HDU: `%s').", demo->name);
                  gal_timing_report(NULL, msg, 2);
                  free(msg);
                }
              break;

            case 2:
              demo=p->clabel;
              demo->name = "DET_CLUMPS_TRUE";
              if(!p->cp.quiet)
                {
                  asprintf(&msg, "True clumps found                  "
                           "(HDU: `%s').", demo->name);
                  gal_timing_report(NULL, msg, 2);
                  free(msg);
                }
              break;

            case 3:
              demo=p->olabel;
              demo->name = "DET_CLUMPS_GROWN";
              if(!p->cp.quiet)
                {
                  gal_timing_report(NULL, "Starting to identify objects.", 1);
                  asprintf(&msg, "True clumps grown                  "
                           "(HDU: `%s').", demo->name);
                  gal_timing_report(NULL, msg, 2);
                  free(msg);
                }
              break;

            case 4:
              demo=p->olabel;
              demo->name = "DET_OBJ_IDENTIFIED";
              if(!p->cp.quiet)
                {
                  asprintf(&msg, "Identified objects over detections "
                           "(HDU: `%s').", demo->name);
                  gal_timing_report(NULL, msg, 2);
                  free(msg);
                }
              break;

            case 5:
              demo=p->olabel;
              demo->name = "DET_OBJECTS_FULL";
              if(!p->cp.quiet)
                {
                  asprintf(&msg, "Objects grown to cover full area   "
                           "(HDU: `%s').", demo->name);
                  gal_timing_report(NULL, msg, 2);
                  free(msg);
                }
              break;

            case 6:
              demo=p->clabel;
              demo->name = "CLUMPS_FINAL";
              if(!p->cp.quiet)
                {
                  asprintf(&msg, "Clumps given their final label     "
                           "(HDU: `%s').", demo->name);
                  gal_timing_report(NULL, msg, 2);
                  free(msg);
                }
              break;

            case 7:
              demo=p->olabel;
              demo->name = "OBJECTS_FINAL";
              if(!p->cp.quiet)
                {
                  asprintf(&msg, "Objects given their final label    "
                           "(HDU: `%s').", demo->name);
                  gal_timing_report(NULL, msg, 2);
                  free(msg);
                }
              break;

            default:
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so "
                    "we can address the issue. The value %d is not "
                    "recognized for clprm.step", __func__, PACKAGE_BUGREPORT,
                    clprm.step);
            }

          /* Write the demonstration array into the check image. The
             default values are hard to view, so we'll make a copy of the
             demo, set all Sky regions to blank and all clump macro values
             to zero. */
          gal_fits_img_write(demo, p->segmentationname, NULL, PROGRAM_NAME);

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
      p->olabel->name = p->clabel->name = NULL;
    }
  else
    {
      clprm.step=0;
      gal_threads_spin_off(segmentation_on_threads, &clprm, p->numdetections,
                           p->cp.numthreads);
    }


  /* Save the final number of objects and clumps. */
  p->numclumps=clprm.totclumps;
  p->numobjects=clprm.totobjects;

  /* If the user wanted to see the S/N table, then make the S/N table. */
  if(p->checkclumpsn) segmentation_save_sn_table(&clprm);


  /* Clean up allocated structures and destroy the mutex. */
  gal_data_array_free(clprm.sn, p->numdetections+1, 1);
  gal_data_array_free(labindexs, p->numdetections+1, 1);
  if( p->cp.numthreads>1 ) pthread_mutex_destroy(&clprm.labmutex);
}




















/***********************************************************************/
/*****************         High level function         *****************/
/***********************************************************************/
void
segmentation(struct noisechiselparams *p)
{
  float *f;
  char *msg;
  int32_t *l, *lf;
  struct timeval t1;

  /* To keep the user up to date. */
  if(!p->cp.quiet)
    {
      if(!p->cp.quiet) gettimeofday(&t1, NULL);
      gal_timing_report(NULL, "Starting segmentation.",
                        1);
    }


  /* If a check segmentation image was requested, then put in the
     inputs. */
  if(p->segmentationname)
    {
      gal_fits_img_write(p->input, p->segmentationname, NULL, PROGRAM_NAME);
      gal_fits_img_write(p->conv, p->segmentationname, NULL, PROGRAM_NAME);
      p->olabel->name="DETECTION_LABELS";
      gal_fits_img_write(p->olabel, p->segmentationname, NULL,
                         PROGRAM_NAME);
      p->olabel->name=NULL;
    }


  /* Allocate the clump labels image. */
  p->clabel=gal_data_alloc(NULL, p->olabel->type, p->olabel->ndim,
                           p->olabel->dsize, p->olabel->wcs, 1,
                           p->cp.minmapsize, NULL, NULL, NULL);
  p->clabel->flag=p->input->flag;


  /* Set any possibly existing NaN values to blank. */
  f=p->input->array; lf=(l=p->clabel->array)+p->clabel->size;
  do if(isnan(*f++)) *l=GAL_BLANK_INT32; while(++l<lf);


  /* Find the clump S/N threshold. */
  clumps_true_find_sn_thresh(p);


  /* Reset the clabel array to find true clumps in objects. */
  f=p->input->array; lf=(l=p->clabel->array)+p->clabel->size;
  do *l = isnan(*f++) ? GAL_BLANK_INT32 : 0; while(++l<lf);


  /* Find true clumps over the detected regions. */
  segmentation_detections(p);


  /* Report the results and timing to the user. */
  if(!p->cp.quiet)
    {
      asprintf(&msg, "%zu object%s""containing %zu clump%sfound.",
               p->numobjects, p->numobjects==1 ? " " : "s ",
               p->numclumps,  p->numclumps ==1 ? " " : "s ");
      gal_timing_report(&t1, msg, 1);
      free(msg);
    }


  /* If the user wanted to check the segmentation and hasn't called
     `continueaftercheck', then stop NoiseChisel. */
  if(p->segmentationname && !p->continueaftercheck)
    ui_abort_after_check(p, p->segmentationname, NULL,
                         "showing all segmentation steps");
}
