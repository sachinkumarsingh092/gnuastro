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
#include <config.h>

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <stdlib.h>

#include "forqsort.h"
#include "neighbors.h"
#include "arraymanip.h"
#include "linkedlist.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "clumps.h"
#include "binary.h"









/****************************************************************
 *****************   Over segmentation       ********************
 ****************************************************************/
/* This is very similar to the immersion method of Vincent &
   Soille(1991), but I will not separate the image into layers, and
   work based on the ordered flux values. If a certain pixel (at a
   certain level) has no neighbors it is a local maximum and will be
   assigned a new label. If it has a labeled neighbor, it will take
   that label and if there is more than one neighboring labeled region
   that pixel will be a `river` pixel.

   x0, y0, x1 and y1:

     - For the noise (mesh boxs), they specify the region in the image
       where this over-segmentation is taking place.

     - For the detections, they should be set to 0, 0, is0 and is1,
       where is0 and is1 are the height and width of the input image.
*/
void
oversegment(struct noisechiselparams *p, size_t *inds, size_t area,
            size_t x0, size_t y0, size_t x1, size_t y1,
            size_t *numclumps)
{
  /* pix is not actually used its self, however, the pointer to it
     will be extensively used (ind) for the macro FILL_NGB_8_REGION.
     This macro works on the pointer to the index, not the index its
     self. `pix` is filled with different index values, so the pointer
     to it doesn't change. */
  size_t pix;

  float *arr=p->conv;
  size_t *topinds=p->topinds;
  struct sll *Q=NULL, *cleanup=NULL;
  long n1, rlab, nlab, curlab=1, *clab=p->clab;
  size_t *n, *nf, *indf, *pind, *ind=&pix, is1=p->lmp.s1;
  size_t ng, *rn, *rnf, numngb, ngb[8], *relngb=p->relngb;

  /* Initialize the region you want to over-segment. */
  indf=(pind=inds)+area;
  do clab[*pind]=SEGMENTINIT; while(++pind<indf);

  /* In the case where a connected region with the same flux or masked
     regions exists, some later indices might already be labeled. Note
     that in the convolved image that is being used here, the masked
     pixels have the smallest possible float value.  */
  indf=(pind=inds)+area;
  do
    /* When regions of a constant flux or masked regions exist, some
       later indexs (although they have same flux) will be filled
       before hand. If they are done, there is no need to do them
       again. */
    if(clab[*pind]==SEGMENTINIT)
      {
        /* View the over-segmentation pixel by pixel. Make sure that
           `--viewsegtf` is off, or else several copies will be made.

        if(p->data0_noise1==0 && p->mid==81 && pind!=p->indexs)
          {
            long *tmp;
            float *ftmp;
            static size_t counter=0;
            size_t tmpx0=495, tmpy0=499, tmpx1=510, tmpy1=508;

            if(counter==0)
              {
                floatshrinkarraytonew(arr, np->lmesh.s0, np->lmesh.s1,
                                      tmpx0, tmpy0, tmpx1, tmpy1, &ftmp);
                array_to_fits("overseg.fits", NULL, "1", FLOAT_IMG,
                              ftmp, tmpx1-tmpx0, tmpy1-tmpy0, NULL);
                free(ftmp);
                counter++;
              }

            if( *(pind-1)/is1>=tmpx0 && *(pind-1)/is1<tmpx1
                && *(pind-1)%is1>=tmpy0 && *(pind-1)%is1<tmpy1
                && counter<=(tmpx1-tmpx0)*(tmpy1-tmpy0) )
              {
                longshrinkarraytonew(clab, np->lmesh.s0, np->lmesh.s1,
                                     tmpx0, tmpy0, tmpx1, tmpy1, &tmp);
                array_to_fits("overseg.fits", NULL, "1", LONG_IMG,
                              tmp, tmpx1-tmpx0, tmpy1-tmpy0, NULL);
                free(tmp);
                counter++;
              }
          }
        */
        /* Some cases might happen where one or multiple regions of
           the pixels under study have the same flux. In particular
           note that masked pixels were all given a value of
           FLT_MAX. We have sorted the pixels by flux. So two equal
           valued pixels of two separate (but equal flux) regions
           might fall immediately after each other (For example two
           nearby stars whose centers are masked and are initially
           detected as one object because their wings touch above the
           noise).

           Therefore, if we see that the next pixel in the index list
           has the same flux as this one, it does not guarantee that
           it should be given the same label. Similar to the breadth
           first search algorithm for finding connected components, we
           will search all the neighbours and the neighbours of those
           neighbours that have the same flux of this pixel to see if
           they touch any label or not and to finally give them all
           the same label. */
        if( pind+1<indf && arr[*pind]==arr[*(pind+1)] )
          {
            n1=0;
            if(Q!=NULL || cleanup!=NULL)
              error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we "
                    "can fix this problem. In oversegment (clumps.c) Q and "
                    "cleanup should be NULL but while checking the equal "
                    "flux regions they aren't.", PACKAGE_BUGREPORT);
            add_to_sll(&Q, *pind);
            add_to_sll(&cleanup, *pind);
            clab[*pind]=SEGMENTTMPCHECK;

            /* Find all the pixels that have the same flux and are
               connected. */
            while(Q!=NULL)
              {
                /* Pop an element from the queue. */
                pop_from_sll(&Q, ind);

                /* Check the vicinity of this pixel that was just
                   popped to see if it can find any already labeled
                   neighbour or not. */
                FILL_NGB_8_REGION;

                /* If the pixel is on the side of the region, set it
                   as a river, no more need to look around it. */
                if(numngb<8)
                  clab[*ind]=SEGMENTRIVER;
                else
                  {
                    /* Begin looking into the neighbours of this pixel. */
                    nf=(n=ngb)+numngb;
                    do
                      {
                        /* If this neighbour has not been labeled yet
                           and has an equal flux, add it to the queue
                           to expand the studied region.*/
                        nlab=clab[*n];
                        if( nlab == SEGMENTINIT && arr[*n] == arr[*pind] )
                          {
                            clab[*n]=SEGMENTTMPCHECK;
                            add_to_sll(&cleanup, *n);
                            add_to_sll(&Q, *n);
                          }
                        /* If this neighbour has a positive nlab, it
                           means that it belongs to another object, so
                           if n1 has not been set for the whole region
                           put this label equal to n1. If n1 has been
                           set and is different from `nlab' then this
                           whole equal flux region should be a wide
                           river because it is connecting two
                           connected regions.*/
                        else if (nlab>0)
                          {
                            if(n1==0)
                              n1=nlab;
                            else if(nlab!=n1)
                              n1=SEGMENTRIVER;
                          }
                        /* If this neigbour has a label of zero, then
                           we are on the edge of the region. When
                           over-segmenting the noise and the
                           detections, clab is zero for the parts of
                           the image that we are not interested in
                           (detections and noise respectively). */
                        else if(nlab==0)
                          clab[*ind]=SEGMENTRIVER;
                      }
                    while(++n<nf);
                  }
              }

            /* Set the label that is to be given to this equal flux
               region. If n1 was set to any value, then that label
               should be used for the whole region. Otherwise, this is
               a new label. */
            if(n1) rlab = n1;
            else   rlab = curlab++;

            /* Give the same label to the whole connected equal flux
               region, except those that might have been on the side
               of the image and were a river pixel. */
            while(cleanup!=NULL)
              {
                pop_from_sll(&cleanup, &pix);
                /* If it was on the sides of the image, it has been
                   changed to a river pixel. */
                if(clab[pix]==SEGMENTTMPCHECK)
                  clab[pix]=rlab;
              }
          }
        /* The flux of this pixel is not the same as the next sorted
           flux, so simply find the label for this object. */
        else
          {
            /* Check if the pixel is on the side of the image (for
               detections) or mesh box (for noise).*/
            if(*pind/is1==x0 || *pind%is1==y0
               || *pind/is1==x1-1 || *pind%is1==y1-1)
              n1=SEGMENTRIVER;
            else
              {
                /* Set the first neighbour's label to zero. */
                n1=0;
                if(Q!=NULL || cleanup!=NULL)
                  error(EXIT_FAILURE, 0, "A bug! Please contact us at %s "
                        "so we can fix this problem. In oversegment "
                        "(clumps.c) Q and cleanup should be NULL but while "
                        "checking the equal flux regions they aren't.",
                        PACKAGE_BUGREPORT);

                /* Go over all the 8 neighbors of this pixel and see
                   if all the neighbors that have a non-negative value
                   belong to one label or not. If the pixel is
                   neighboured by more than one label, set it as a
                   river pixel. Also if it is touching a zero valued
                   pixel (which does not belong to this object), set
                   it as a river pixel.

                   relngb was defined in ui.c: it keeps the relative
                   indexs of the neighbors of a pixel.*/
                rnf=(rn=relngb)+8;
                do
                  {
                    ng=*pind+*rn;
                    nlab=clab[ng];
                    if(nlab>0)
                      {
                        if(n1==0) n1=nlab;
                        else if(nlab!=n1)
                          {
                            n1=SEGMENTRIVER;
                            break;
                          }
                      }
                    else if(nlab==0)
                      {
                        n1=SEGMENTRIVER;
                        break;
                      }
                  }
                while(++rn<rnf);
              }

            /* Either assign a new label to this pixel, or give it the
               one of its neighbors. If n1 equals zero, then it is a
               new peak, and a new label should be created.  But if
               n1!=0, it is either a river pixel (has more than one
               labeled neighbor) or all its neighbors have the same
               label. */
            if(n1) rlab = n1;
            else
              {
                rlab = curlab++;
                if(topinds) topinds[rlab]=*pind;
              }

            /* Put the found label in the pixel. */
            clab[ *pind ] = rlab;
          }
      }
  while(++pind<indf);

  *numclumps=curlab;
}



















/******************************************************************/
/*************     Find/Apply S/N threshold      ******************/
/******************************************************************/
void *
clumpsnthreshonmesh(void *inparams)
{
  struct meshthreadparams *mtp=(struct meshthreadparams *)inparams;
  struct meshparams *mp=mtp->mp;
  struct noisechiselparams *p=(struct noisechiselparams *)mp->params;

  size_t *mponeforall=mp->oneforall;
  size_t *inds=&mponeforall[mtp->id*mp->maxs0*mp->maxs1];

  size_t numclumps;
  long *clab=p->clab;
  float minbfrac=p->minbfrac;
  char *segmentationname=p->segmentationname;
  size_t i, *indexs=&mp->indexs[mtp->id*mp->thrdcols];
  size_t s0, s1, nf, nb, ind, size, startind, is1=mp->s1;

  for(i=0;indexs[i]!=NONTHRDINDEX;++i)
    {
      /* Set the necesary parameters: */
      ind=indexs[i];
      s0=mp->ts0[mp->types[ind]];
      s1=mp->ts1[mp->types[ind]];
      startind=mp->start[ind];
      size=s0*s1;


      /* Check to see if we have enough blank area for getting the
	 background noise statistics. */
      count_f_b_onregion(p->byt, startind, s0, s1, is1, &nf, &nb);
      if( (float)nb < (float)(size)*minbfrac )
        {
          if(segmentationname)
            longinitonregion(clab, 0, startind, s0, s1, is1);
          continue;
        }


      /* We want to find the clumps on the noise, not the signal, the
         number of noise pixels in this mesh were calculated above
         (nb). Here we want to pull out the indexs of those pixels
         which is necessary for the over-segmentation. */
      index_f_b_onregion(p->byt, startind, s0, s1, is1, inds, 0);


      /* Sort the indexs based on the flux within them. */
      qsort(inds, nb, sizeof(size_t), indexfloatdecreasing);


      /* Do the over-segmentation: */
      oversegment(p, inds, nb, startind/is1, startind%is1,
                  startind/is1+s0, startind%is1+s1, &numclumps);
    }

  /* Free any allocated space and if multiple threads were used, wait
     until all other threads finish. */
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}





void
clumpsngrid(struct noisechiselparams *p)
{
  struct meshparams *lmp=&p->lmp;

  /* Set the convolved image as the basis for sorting the indexs and
     finding clumps. */
  forqsortindexarr=p->conv;

  /* For finding the noise signal to noise, we know there is no
     objects, so there is no need to identify the index with the top
     flux. Later on, when dealing with the detections, they become
     important. */
  p->topinds=NULL;

  /* Find the clump signal to noise on each mesh: */
  operateonmesh(lmp, clumpsnthreshonmesh, sizeof(size_t), 0, 1);

  arraytofitsimg(p->segmentationname, "NoiseOversegmentaion",
                 LONG_IMG, p->clab, p->smp.s0, p->smp.s1, 0, p->wcs,
                 NULL, SPACK_STRING);
}
