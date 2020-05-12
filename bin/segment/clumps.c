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
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/fits.h>
#include <gnuastro/blank.h>
#include <gnuastro/label.h>
#include <gnuastro/threads.h>
#include <gnuastro/pointer.h>
#include <gnuastro/dimension.h>
#include <gnuastro/statistics.h>

#include <gnuastro-internal/timing.h>

#include "main.h"

#include "ui.h"
#include "clumps.h"










/**********************************************************************/
/*****************              Grow clumps           *****************/
/**********************************************************************/
/* Make the preparations for the intiial growing the clumps to identify
   objects: a single standard deviation for the whole object and preparing
   the labels (because the growth is going to happen on the 'olabel'
   image. */
void
clumps_grow_prepare_initial(struct clumps_thread_params *cltprm)
{
  gal_data_t *indexs=cltprm->indexs;
  gal_data_t *input=cltprm->clprm->p->input;
  struct segmentparams *p=cltprm->clprm->p;

  size_t ndiffuse=0, coord[2], *dindexs;
  double wcoord[2]={0.0f,0.0f}, sum=0.0f;
  size_t *s, *sf, *dsize=input->dsize, ndim=input->ndim;
  float glimit, *imgss=input->array, *std=p->std->array;
  int32_t *olabel=p->olabel->array, *clabel=p->clabel->array;


  /* Find the flux weighted center (meaningful only for positive valued
     pixels). */
  sf=(s=indexs->array)+indexs->size;
  do
    if( imgss[ *s ] > 0.0f )
      {
        sum        += imgss[ *s ];
        wcoord[0]  += imgss[ *s ] * (*s/dsize[1]);
        wcoord[1]  += imgss[ *s ] * (*s%dsize[1]);
      }
  while(++s<sf);


  /* Calculate the center, if no pixels were positive, use the
     geometric center (irrespective of flux). */
  if(sum==0.0f)
    {
      sf=(s=indexs->array)+indexs->size;
      do
        {
          wcoord[0] += *s / dsize[1];
          wcoord[1] += *s % dsize[1];
        }
      while(++s<sf);
      sum = indexs->size;
    }


  /* Convert floating point coordinates to integers. */
  coord[0] = GAL_DIMENSION_FLT_TO_INT(wcoord[0]/sum);
  coord[1] = GAL_DIMENSION_FLT_TO_INT(wcoord[1]/sum);


  /* Find the growth limit. Note that the STD may be a value, or a dataset
     (which may be a full sized image or a tessellation). If its not a
     single value, we'll check through the number of elements to see what
     kind of dataset it is (if its a tile or full image). */
  cltprm->std = ( p->std->size>1
                  ? ( p->std->size==p->input->size
                      ? std[gal_dimension_coord_to_index(ndim, dsize, coord)]
                      : std[gal_tile_full_id_from_coord(&p->cp.tl, coord)] )
                  : std[0] );
  if(p->variance) cltprm->std = sqrt(cltprm->std);


  /* From the standard deviation, find the growth limit. */
  glimit = p->gthresh * cltprm->std;


  /* Allocate space to keep the diffuse indexs over this detection. We need
     to keep the actual indexs since it is our only connection to the
     object at this stage: we are also going to re-label the pixels to
     grow. For most astronomical objects, the major part of the detection
     area is going to be diffuse flux, so we will just allocate the same
     size as 'indexs' array (the 'dsize' will be corrected after getting
     the exact number.

     Also note that since 'indexs' is already sorted, therefore
     'diffuseindexs' will also be already sorted. */
  cltprm->diffuseindexs=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1,
                                       cltprm->indexs->dsize, NULL, 0,
                                       p->cp.minmapsize, p->cp.quietmmap,
                                       NULL, NULL, NULL);
  dindexs=cltprm->diffuseindexs->array;
  sf=(s=indexs->array)+indexs->size;
  do
    {
      olabel[*s] = clabel[*s];
      if( clabel[*s]==GAL_LABEL_INIT )
        if( imgss[*s]>glimit ) dindexs[ ndiffuse++ ] = *s;
    }
  while(++s<sf);


  /* Correct the sizes of the 'diffuseindexs' data structure. */
  cltprm->diffuseindexs->size = cltprm->diffuseindexs->dsize[0] = ndiffuse;
}





/* Add all the remaining pixels in the detection (below the growth
   threshold, or those that were not touching). Note that initially
   'diffuseindexs' was filled with the pixels that are above the growth
   threshold. That was necessary for identifying the objects. Now that we
   have identified the objects and labeled them, we want to add the
   remaining diffuse pixels to it too before doing the final growth.

   Note that the most efficient way is just to re-fill the 'diffuseindexs'
   array instead of adding the pixels below the threshold and sorting them
   afterwards.*/
void
clumps_grow_prepare_final(struct clumps_thread_params *cltprm)
{
  size_t ndiffuse=0;
  size_t *dindexs=cltprm->diffuseindexs->array;
  int32_t *olabel=cltprm->clprm->p->olabel->array;
  size_t *s=cltprm->indexs->array, *sf=s+cltprm->indexs->size;

  /* Recall that we initially allocated 'diffuseindexs' to have the same
     size as the indexs. So there is no problem if there are more pixels in
     this final round compared to the initial round. */
  do
    if( olabel[*s] < 0 )
      dindexs[ ndiffuse++ ] = *s;
  while(++s<sf);

  /* Correct the sizes of the 'diffuseindexs' data structure. */
  cltprm->diffuseindexs->size = cltprm->diffuseindexs->dsize[0] = ndiffuse;
}
























/**********************************************************************/
/*****************             S/N threshold          *****************/
/**********************************************************************/
/* Correct the labels of the clumps that will be used in determining the
   S/N threshold for true clumps.   */
static void
clumps_correct_sky_labels_for_check(struct clumps_thread_params *cltprm,
                                    gal_data_t *tile)
{
  gal_data_t *newinds;
  int32_t *ninds, curlab, *l, *lf;
  size_t len=cltprm->numinitclumps+1;
  struct segmentparams *p=cltprm->clprm->p;

  /* If any of the clumps must be kept ('cltprm->snind->size!=0'), then
     re-label them for the check image. Otherwise, remove all clumps. */
  if(cltprm->snind->size)
    {
      /* A small sanity check. */
      if(gal_tile_block(tile)!=p->clabel)
        error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to "
              "address the problem. 'tile->block' must point to the "
              "'clabel' dataset", __func__, PACKAGE_BUGREPORT);


      /* Allocate a dataset with the new indexs, note that it will need to
         have one element for each initial label (the excluded clumps need
         to be set to zero). So we also need to clear the allocated
         space. */
      newinds=gal_data_alloc(NULL, p->clabel->type, 1, &len, NULL, 0,
                             p->cp.minmapsize, p->cp.quietmmap,
                             NULL, NULL, NULL);


      /* Get the next available label for these clumps. If more than one
         thread was used, we are first going to lock the mutex (so no other
         thread changes these values), we will then read the shared number
         for this thread to use, then update the shared number and finally,
         unlock the mutex so other threads can do the same when they get to
         this point. */
      if(p->cp.numthreads>1) pthread_mutex_lock(&cltprm->clprm->labmutex);
      curlab        = p->numclumps+1; /* Note that counting begins from 1. */
      p->numclumps += cltprm->snind->size;
      if(p->cp.numthreads>1) pthread_mutex_unlock(&cltprm->clprm->labmutex);


      /* Initialize the newinds array to GAL_LABEL_INIT (which be used as a
         new label for all the clumps that must be removed. */
      lf = (l=newinds->array) + newinds->size;
      do *l++=GAL_LABEL_INIT; while(l<lf);


      /* The new indexs array has been initialized to zero. So we just need
         to go over the labels in 'cltprm->sninds' and give them a value of
         'curlab++'. */
      ninds=newinds->array;
      lf = (l=cltprm->snind->array) + cltprm->snind->size;
      do { ninds[*l]=curlab++; *l=ninds[*l]; } while(++l<lf);


      /* Go over this tile and correct the values. */
      GAL_TILE_PARSE_OPERATE( tile, NULL, 0, 1,
                              {if(*i>0) *i=ninds[ *(int32_t *)i ];} );

      /* Clean up. */
      gal_data_free(newinds);
    }
  else
    /* There were no usable clumps in this tile, so just set all the pixels
       larger than zero (a clump) to 'GAL_LABEL_INIT'. */
    GAL_TILE_PARSE_OPERATE( tile, NULL, 0, 1, {*i=*i>0?GAL_LABEL_INIT:*i;} );
}





static void *
clumps_find_make_sn_table(void *in_prm)
{
  struct gal_threads_params *tprm=(struct gal_threads_params *)in_prm;
  struct clumps_params *clprm=(struct clumps_params *)(tprm->params);
  struct segmentparams *p=clprm->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  void *tarray;
  double numdet;
  gal_data_t *tile, *tblock, *tmp;
  uint8_t *binary=p->binary->array;
  struct clumps_thread_params cltprm;
  size_t i, c, ind, tind, num, numsky, *indarr;
  size_t *scoord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                       "scoord");
  size_t *icoord=gal_pointer_allocate(GAL_TYPE_SIZE_T, ndim, 0, __func__,
                                       "icoord");


  /* Initialize the parameters for this thread. */
  cltprm.clprm   = clprm;
  cltprm.topinds = NULL;


  /* Go over all the tiles/detections given to this thread. */
  for(i=0; tprm->indexs[i] != GAL_BLANK_SIZE_T; ++i)
    {
      /* IDs. */
      cltprm.id = tind  = tprm->indexs[i];
      tile = &p->ltl.tiles[tind];


      /* Change the tile's pointers to the binary image (which has 1 for
         detected pixels and 0 for un-detected regions). */
      tarray=tile->array;
      tblock=tile->block;
      tile->array = gal_tile_block_relative_to_other(tile, p->binary);
      tile->block = p->binary;


      /* Get the number of usable elements in this tile (note that tiles
         can have blank pixels), so we can't simply use 'tile->size'. */
      if(p->input->flag & GAL_DATA_FLAG_HASBLANK)
        {
          tmp=gal_statistics_number(tile);
          num=*((size_t *)(tmp->array));
          gal_data_free(tmp);
        }
      else num=tile->size;


      /* Find the number of detected pixels over this tile. Since this is
         the binary image, this is just the sum of all the pixels.

         Note that 'numdet' can be 'nan' when the whole tile is blank and
         so there was no values to sum. Recall that in summing, when there
         is not input, the output is 'nan'. */
      tmp=gal_statistics_sum(tile);
      numdet=*((double *)(tmp->array));
      gal_data_free(tmp);


      /* See if this tile should be used or not (has enough undetected
         pixels). Note that it might happen that some tiles are fully
         blank. In such cases, it is important to first check the number of
         detected pixels. */
      numsky=num-numdet;
      if( num && (float)numsky/(float)num > p->minskyfrac )
        {
          /* Add the indexs of all undetected pixels in this tile into an
             array. */
          cltprm.indexs=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &numsky,
                                       NULL, 0, p->cp.minmapsize,
                                       p->cp.quietmmap, NULL, NULL, NULL);


          /* Change the tile's block to the clump labels dataset (because
             we'll need to set the labels of the rivers on the edge of the
             tile here). */
          tile->array = gal_tile_block_relative_to_other(tile, p->clabel);
          tile->block = p->clabel;


          /* We need to set all the pixels on the edge of the tile to
             rivers and not include them in the list of indexs to set
             clumps. To do that, we need this tile's starting
             coordinates. */
          gal_dimension_index_to_coord(gal_pointer_num_between(
                           p->clabel->array, tile->array, p->clabel->type),
                                       ndim, dsize, scoord);


          /* Add the index of every sky element to the array of
             indexs. Note that since we know the array is always of type
             'int32_t', we can call the 'GAL_TILE_PO_OISET' macro to avoid
             having to deal with multiple possible types in
             'GAL_TILE_PARSE_OPERATE'. Since the OUT macro-variable is
             NULL, the 'int' is just a place-holder, it will not be
             used. */
          c=0;
          indarr=cltprm.indexs->array;
          GAL_TILE_PO_OISET(int32_t, int, tile, NULL, 0, 1, {
              /* This pixel's index over all the image. */
              ind = (int32_t *)i - (int32_t *)(p->clabel->array);
              gal_dimension_index_to_coord(ind, ndim, dsize, icoord);

              /* If the pixel is on the tile edge, set it as river and
                 don't include it in the indexs. */
              if( icoord[0]==scoord[0]
                  || icoord[0]==scoord[0]+tile->dsize[0]-1
                  || icoord[1]==scoord[1]
                  || icoord[1]==scoord[1]+tile->dsize[1]-1 )
                *(int32_t *)i=GAL_LABEL_RIVER;

              /* This pixel is not on the edge, check if it had a value of
                 '0' in the binary image (is not detected) then add it to
                 the list of indexs (note that the binary image also
                 contains the blank pixels, so only sky regions have a
                 value of 0 in the binary image). */
              else if( binary[ind]==0 )
                {
                  /*
                  if(c!=cltprm.indexs->size)
                    {
                      if(cltprm.id==282) *i+=2;
                  */
                      indarr[c++]=gal_pointer_num_between(p->clabel->array,
                                                          i, p->clabel->type);
                  /*
                    }
                  else
                    if(cltprm.id==282)
                      {
                        int32_t *clabel=p->clabel->array;
                        size_t kjd=gal_data_num_between(p->clabel->array, i,
                                                        p->clabel->type);
                        printf("%zu, %zu: %u\n", kjd%dsize[1]+1,
                               kjd/dsize[1]+1, clabel[kjd]);
                      }
                  */
                }
            });


          /* Correct the number of indexs. */
          cltprm.indexs->size=cltprm.indexs->dsize[0]=c;


          /* Generate the clumps over this region. */
          cltprm.numinitclumps=gal_label_watershed(p->conv, cltprm.indexs,
                                                   p->clabel,
                                                   cltprm.topinds,
                                                   !p->minima);


          /* Set all river pixels to GAL_LABEL_INIT (to be distinguishable
             from the detected regions). */
          GAL_TILE_PO_OISET( int32_t, int, tile, NULL, 0, 1,
                             {if(*i==GAL_LABEL_RIVER) *i=GAL_LABEL_INIT;} );


          /* For a check, the step variable will be set. */
          if(clprm->step==1)
            { gal_data_free(cltprm.indexs); continue; }


          /* Make the clump S/N table. */
          cltprm.sn    = &cltprm.clprm->sn[cltprm.id];
          cltprm.snind = ( cltprm.clprm->snind
                           ? &cltprm.clprm->snind[cltprm.id]
                           : NULL );
          gal_label_clump_significance(p->clumpvals, p->std, p->clabel,
                                       cltprm.indexs, &p->cp.tl,
                                       cltprm.numinitclumps, p->snminarea,
                                       p->variance, clprm->sky0_det1,
                                       cltprm.sn, cltprm.snind);


          /* If the user wanted to check the steps, remove the clumps that
             weren't used from the 'clabel' image (they have been already
             excluded from the table). */
          if(cltprm.snind)
            clumps_correct_sky_labels_for_check(&cltprm, tile);


          /* If there were no clumps, then just set the S/N table to
             NULL. This must be done after the check image creation (if
             necessary), because we use 'cltprm.snind' as a proxy for the
             check image.*/
          if( cltprm.clprm->sn[ cltprm.id ].size==0 )
            cltprm.snind=cltprm.sn=NULL;


          /* Clean up. */
          gal_data_free(cltprm.indexs);
        }

      /* Reset the tile's pointers back to what they were. */
      tile->array=tarray;
      tile->block=tblock;
    }

  /* Clean up. */
  free(scoord);
  free(icoord);

  /* Wait for the all the threads to finish and return. */
  if(tprm->b) pthread_barrier_wait(tprm->b);
  return NULL;
}





/* Write the S/N table. */
static void
clumps_write_sn_table(struct segmentparams *p, gal_data_t *insn,
                      gal_data_t *inind, char *filename,
                      gal_list_str_t *comments)
{
  gal_data_t *sn, *ind, *cols;

  /* Remove all blank elements. The index and sn values must have the same
     set of blank elements, but checking on the integer array is faster. */
  if( gal_blank_present(inind, 1) )
    {
      /* Remove blank elements. */
      ind=gal_data_copy(inind);
      sn=gal_data_copy(insn);
      gal_blank_remove(ind);
      gal_blank_remove(sn);

      /* A small sanity check. */
      if(ind->size==0 || sn->size==0)
        error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
              "the problem. For some reason, all the elements in 'ind' or "
              "'sn' are blank", __func__, PACKAGE_BUGREPORT);
    }
  else
    {
      sn  = insn;
      ind = inind;
    }

  /* Set the columns. */
  cols       = ind;
  cols->next = sn;

  /* Prepare the comments. */
  gal_table_comments_add_intro(&comments, PROGRAM_STRING, &p->rawtime);

  /* write the table. */
  gal_table_write(cols, comments, p->cp.tableformat, filename,
                  "SKY_CLUMP_SN", 0);

  /* Clean up (if necessary). */
  if(sn!=insn) gal_data_free(sn);
  if(ind==inind) ind->next=NULL; else gal_data_free(ind);
}





/* Find the true clump signal to noise value from the clumps in the sky
   region.

   Each thread will find the useful signal to noise values for the tiles
   that have been assigned to it. It will then store the pointer to the S/N
   table into the sntablearr array (with the size of the number of
   meshs). If no clumps could be found in a mesh, then
   sntablearr[i]=NULL. Otherwise, it points to an array of the useful S/N
   values in that clump. Note that we don't care about the order of S/N
   values any more! There is also an accompanying array to keep the number
   of elements in the final S/N array of each mesh: numclumpsarr.

   Using these two arrays, after all the threads are finished, we can
   concatenate all the S/N values into one array and send it to the main
   findsnthresh function in thresh.c. */
void
clumps_true_find_sn_thresh(struct segmentparams *p)
{
  char *msg;
  struct timeval t1;
  size_t i, j, c, numsn=0;
  struct clumps_params clprm;
  gal_list_str_t *comments=NULL;
  gal_data_t *sn, *snind, *quant, *claborig;

  /* Get starting time for later reporting if necessary. */
  if(!p->cp.quiet) gettimeofday(&t1, NULL);


  /* Initialize/allocate the clump parameters structure, Note that the S/N
     indexs are also needed when we want to check the segmentation steps
     (they are used to correct the indexs in the final output). */
  clprm.p=p;
  clprm.sky0_det1=0;
  clprm.sn=gal_data_array_calloc(p->ltl.tottiles);
  clprm.snind = ( p->checksegmentation || p->checksn
                  ? gal_data_array_calloc(p->ltl.tottiles) : NULL );


  /* If the user wants to check the steps of get an S/N table, then we need
     a unique label for each clump. But in each region, the labels start
     from 1. So we need a central place to keep the next available
     label. Since 'p->numclumps' is not used yet, we will use it here. When
     multiple threads are used, we will need a mutex to make sure that only
     one thread can change this central variable at every one moment. */
  if(p->checksegmentation || p->checksn)
    {
      p->numclumps=0;
      if( p->cp.numthreads > 1 ) pthread_mutex_init(&clprm.labmutex, NULL);
    }


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

          /* Do this step. */
          gal_threads_spin_off(clumps_find_make_sn_table, &clprm,
                               p->ltl.tottiles, p->cp.numthreads);

          /* Set the extension name. */
          switch(clprm.step)
            {
            case 1: p->clabel->name = "SKY_CLUMPS_ALL";    break;
            case 2: p->clabel->name = "SKY_CLUMPS_FOR_SN"; break;
            default:
              error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s so "
                    "we can address the issue. The value %d is not valid for "
                    "clprm.step", __func__, PACKAGE_BUGREPORT, clprm.step);
            }

          /* Write the demonstration array into the check image. The
             default values are hard to view, so we'll make a copy of the
             demo, set all Sky regions to blank and all clump macro values
             to zero. */
          gal_fits_img_write(p->clabel, p->segmentationname, NULL,
                             PROGRAM_NAME);

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
      gal_threads_spin_off(clumps_find_make_sn_table, &clprm,
                           p->ltl.tottiles, p->cp.numthreads);
    }


  /* Destroy the mutex if it was initialized. */
  if( p->cp.numthreads>1 && (p->checksegmentation || p->checksn) )
    pthread_mutex_destroy(&clprm.labmutex);


  /* Find the total number of S/N values we have for all the clumps. */
  for(i=0;i<p->ltl.tottiles;++i)
    if(clprm.sn[i].ndim)  /* Only on tiles were an S/N was calculated. */
      numsn+=clprm.sn[i].size;
  if( numsn < p->minnumfalse )
    error(EXIT_FAILURE, 0, "%zu usable clumps found in the undetected "
          "regions. This is smaller than the requested minimum number of "
          "false/reference clumps (%zu, value to the '--minnumfalse' "
          "option).\n\n"
          "There are several ways to address the problem. The best and most "
          "highly recommended is to use a larger input if possible (when the "
          "input is a crop from a larger dataset). If that is not the case, "
          "or it doesn't solve the problem, you need to loosen the "
          "parameters (and therefore cause more scatter/bias in the final "
          "result). Thus don't loosen them too much. Recall that you can "
          "see all the option values to Gnuastro's programs by appending "
          "'-P' to the end of your command.\n\n"
          "  * Slightly decrease '--largetilesize' to have more tiles.\n"
          "  * Decrease '--minskyfrac' (currently %g) to look into more "
          "tiles.\n"
          "  * Slightly decrease '--snminarea' (currently %zu) to "
          "measure more clumps.\n"
          "  * If Segment already works on a dataset with similar noise "
          "properties, you can directly pass the 'true' clump "
          "signal-to-noise ratio found there to '--clumpsnthresh' and "
          "avoid having to study the undetected regions any more.\n\n"
          "Append your previous command with '--checksegmentation' to see "
          "the steps and get a better feeling of the cause/solution. Note "
          "that the output is a multi-extension FITS file).\n\n"
          "To better understand the segmentation process and options, "
          "please run the following command (press 'SPACE'/arrow-keys to "
          "navigate and 'Q' to return back to the command-line):\n\n"
          "    $ info gnuastro \"Segmentation options\"\n",
          numsn, p->minnumfalse, p->minskyfrac, p->snminarea);


  /* Allocate the space to keep all the S/N values. */
  sn=gal_data_alloc(NULL, GAL_TYPE_FLOAT32, 1, &numsn, NULL, 0,
                    p->cp.minmapsize, p->cp.quietmmap, "CLUMP_S/N", "ratio",
                    "Signal-to-noise ratio");
  snind = ( p->checksn
            ? gal_data_alloc(NULL, GAL_TYPE_INT32, 1, &numsn, NULL, 0,
                             p->cp.minmapsize, p->cp.quietmmap, "CLUMP_ID",
                             "counter", "Unique ID for this clump.")
            : NULL );


  /* Copy the S/N values of all the clumps into the unified array. */
  c=0;
  for(i=0;i<p->ltl.tottiles;++i)
    if(clprm.sn[i].ndim)
      for(j=0;j<clprm.sn[i].size;++j)
        {
          ((float *)(sn->array))[c] = ((float *)(clprm.sn[i].array))[j];
          if(snind)
            ((int32_t *)(snind->array))[c] =
              ((int32_t *)(clprm.snind[i].array))[j];
          ++c;
        }


  /* The S/N array of sky clumps is desiged to have no blank values, so set
     the flags accordingly to avoid a redundant blank search. */
  sn->flag  =  GAL_DATA_FLAG_BLANK_CH;
  sn->flag &= ~GAL_DATA_FLAG_HASBLANK;


  /* If the user wanted to see the S/N table, then save it. */
  if(p->checksn)
    {
      /* Make the comments, then write the table and free the comments. */
      if(p->cp.numthreads>1)
        gal_list_str_add(&comments, "NOTE: In multi-threaded mode, clump "
                         "IDs differ in each run and are not sorted.", 1);
      gal_list_str_add(&comments, "See also: 'SKY_CLUMPS_FOR_SN' HDU of "
                       "output with '--checksegmentation'.", 1);
      gal_list_str_add(&comments, "S/N of clumps over undetected regions.",
                       1);
      clumps_write_sn_table(p, sn, snind, p->clumpsn_s_name, comments);
      gal_list_str_free(comments, 1);
    }


  /* Find the desired quantile from the full S/N distribution. */
  quant = gal_statistics_quantile(sn, p->snquant, 1);
  p->clumpsnthresh = *((float *)(quant->array));
  if(!p->cp.quiet)
    {
      if( asprintf(&msg, "Clump peak S/N: %g (%.3f quant of %zu).",
                   p->clumpsnthresh, p->snquant, sn->size)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(&t1, msg, 2);
      free(msg);
    }


  /* Clean up. */
  gal_data_free(sn);
  gal_data_free(snind);
  gal_data_free(quant);
  gal_data_array_free(clprm.sn, p->ltl.tottiles, 1);
  gal_data_array_free(clprm.snind, p->ltl.tottiles, 1);
}


















/***********************************************************************/
/*****************           Over detections           *****************/
/***********************************************************************/
/* Only keep true clumps over detections. */
void
clumps_det_keep_true_relabel(struct clumps_thread_params *cltprm)
{
  struct segmentparams *p=cltprm->clprm->p;
  size_t ndim=p->input->ndim, *dsize=p->input->dsize;

  int istouching;
  size_t i, *s, *sf, *dinc;
  float *sn = cltprm->sn ? cltprm->sn->array : NULL;
  int32_t *l, *lf, *newlabs, curlab=1, *clabel=p->clabel->array;

  /* If there were no clumps over the detection, then just set the number
     of true clumps to zero, otherwise, see which ones should be
     removed. */
  if(cltprm->sn)
    {
      /* Allocate the necessary arrays. */
      newlabs=gal_pointer_allocate(GAL_TYPE_INT32,
                                   cltprm->numinitclumps+1, 0, __func__,
                                   "newlabs");
      dinc=gal_dimension_increment(ndim, dsize);

      /* Initialize the new labels with GAL_LABEL_INIT (so the diffuse area
         can be distinguished from the clumps). */
      lf=(l=newlabs)+cltprm->numinitclumps+1;
      do *l++=GAL_LABEL_INIT; while(l<lf);

      /* Set the new labels. Here we will also be removing clumps with a peak
         that touches a river pixel. */
      if(p->keepmaxnearriver)
        {
          for(i=1;i<cltprm->numinitclumps+1;++i)
            if( sn[i] > p->clumpsnthresh ) newlabs[i]=curlab++;
        }
      else
        {
          for(i=1;i<cltprm->numinitclumps+1;++i)
            {
              /* Check if all the neighbors of this top element are
                 touching a river or not. */
              istouching=0;
              GAL_DIMENSION_NEIGHBOR_OP(cltprm->topinds[i], ndim, dsize,
                                        ndim, dinc,
                                        {
                                          if(clabel[nind]==0)
                                            istouching=1;
                                        });

              /* If the peak isn't touching a river, then check its S/N and
                 if that is also good, give it a new label. */
              if( !istouching && sn[i] > p->clumpsnthresh )
                newlabs[i]=curlab++;
            }
        }

      /* Correct the clump labels. Note that the non-clumpy regions over
         the detections (rivers) have already been initialized to
         GAL_LABEL_INIT (which is negative). So we'll just need to correct
         the ones with a value larger than 0. */
      sf=(s=cltprm->indexs->array)+cltprm->indexs->size;
      do if(clabel[*s]>0) clabel[*s] = newlabs[ clabel[*s] ]; while(++s<sf);

      /* Save the total number of true clumps in this detection. */
      cltprm->numtrueclumps=curlab-1;

      /* Clean up. */
      free(dinc);
      free(newlabs);
    }
  else cltprm->numtrueclumps=0;
}
