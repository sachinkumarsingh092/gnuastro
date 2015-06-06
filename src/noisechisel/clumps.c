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
#include <string.h>
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
/*************             Clump S/N             ******************/
/******************************************************************/
/* In this function we want to find the general information for each
   clump in an over-segmented labeled array. The signal in each clump
   is the average signal inside it subtracted by the average signal in
   the river pixels around it. So this function will go over all the
   pixels in the object (already found in deblendclumps()) and add
   them appropriately.

   The output is an array of size numseg*INFOTABCOLS. INFOTABCOLS=4. The
   columns are:
   0: Average signal in clump.
   1: Number of pixels in clump.
   2: Average signal around clump.
   3: Number of pixels around clump.
   4: Standard deviation on flux weighted center of clump.
*/
void
getclumpinfo(struct noisechiselparams *p, size_t *inds, size_t area,
             size_t numclumps, size_t x0, size_t y0, size_t x1,
             size_t y1, double **outclumpinfo)
{
  struct meshparams *smp=&p->smp;

  double *clumpinfo;
  float *img=p->img, *smpstd=smp->garray2;
  long *clab=p->clab, wngb[WNGBSIZE], ngblab;
  size_t *xys, index, lab, is0=p->lmp.s0, is1=p->lmp.s1;
  size_t i=0, ii=0, row, *n, *nf, ngb[8], *ind, *indf, numngb;

  /* Just make sure that the box size is not only around one pixel! */
  if(x1-x0<=1 || y1-y0<=1)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can find "
          "and fix the problem in clumpinfo (clumps.c). For some reason, "
          "the specified input region is %lu by %lu wide.",
          PACKAGE_BUGREPORT, y1-y0, x1-x0);

  /* Allocate the clump information array. */
  errno=0;
  clumpinfo=*outclumpinfo=calloc(numclumps*INFOTABCOLS, sizeof *clumpinfo);
  if(clumpinfo==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for clumpinfo (clumps.c)",
          numclumps*INFOTABCOLS*sizeof *clumpinfo);
  errno=0; xys=calloc(2*numclumps, sizeof *xys);
  if(xys==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for xys (clumps.c)",
          2*numclumps*sizeof *xys);

  /* Go over all the pixels in this set of pixels and fill in the
     proper information for each clump. */
  indf = (ind=inds) + area;
  do
    if(!isnan(img[*ind]))
      {
        if(clab[*ind]==SEGMENTRIVER)
          {
            /* Fill in the neighbours arrary for this pixel. If we are
               working on the mesh grid (the noise), then we only want
               the neighbors within a region. Otherwise (when working
               on the detections) we want the neighbors on the full
               image.*/
            if(p->b0f1) {FILL_NGB_8_ALLIMG}
            else        {FILL_NGB_8_REGION}


            /* We are on a river pixel. So its value has to be added
               to the borders of any object it touches. But since it
               might touch a labeled region more than once, we use
               wngb to keep track of which label we have already added
               its value to. `ii` is the number of different labels
               this river pixel has already been added to. wngb will
               keep the labels. */
            ii=0;
            memset(wngb, 0, sizeof(wngb));

            /* Look into the 8-connected neighbors: */
            nf=(n=ngb)+numngb;
            do
              if( (ngblab=clab[ *n ]) >0)
                {
                  /* Go over wngb to see if this river pixel's value
                     has been added to a segment or not. */
                  for(i=0;i<ii;++i)
                    if(wngb[i]==ngblab)
                      /* It is already present. break out. */
                      break;

                  if(i==ii) /* This label was not added yet. */
                    {
                      clumpinfo[ ngblab*INFOTABCOLS+2 ]+=img[*ind];;
                      ++clumpinfo[ ngblab*INFOTABCOLS+3 ];
                      wngb[ii]=ngblab;
                      ++ii;
                    }
                }
            while(++n<nf);
          }
        else
          {
            lab=clab[*ind];     /* The label of this clump. */
            xys[ 2 * lab     ] = *ind/is1;
            xys[ 2 * lab + 1 ] = *ind%is1;
            clumpinfo[ lab * INFOTABCOLS     ]+=img[*ind];
            clumpinfo[ lab * INFOTABCOLS + 1 ]+=1;
          }
      }
  while(++ind<indf);

  /* Do the final preparations. All the calculations are only
     necessary for the clumps that satisfy the minimum area. So there
     is no need to waste time on the smaller ones. */
  for(lab=1;lab<numclumps;++lab)
    {
      row=lab*INFOTABCOLS;
      if(clumpinfo[row+1]>p->segsnminarea)
        {
          /* Convert sum to average: */
          clumpinfo[row  ] /= clumpinfo[row+1];
          clumpinfo[row+2] /= clumpinfo[row+3];

          /* Find the index of the flux weighted center and use it to
             find the standard deviation for this clump. Note that
             this is only needed if the input image was already sky
             subtracted. If it wasn't, then we are not subtracting the
             sky to worry about its error! The error in the total flux
             is simply its square root. */
          if(p->skysubtracted)
            {
              index = ( (xys[2*lab] / clumpinfo[row+1]) * is1
                        + xys[2*lab+1] / clumpinfo[row+1] );
              clumpinfo[row+4]  = smpstd[imgindextomeshid(smp, index)];
            }
        }
    }

  /* Clean up: */
  free(xys);
}





void
clumpsntable(struct noisechiselparams *p, size_t *inds, size_t area,
             size_t *numclumps, size_t x0, size_t y0, size_t x1,
             size_t y1, float **sntable)
{
  float *sntab;
  double *clumpinfo, err;
  size_t i, ind, counter=0;
  double I, O, Ni, cpscorr=p->cpscorr;


  /* Get the information for all the segments. */
  getclumpinfo(p, inds, area, *numclumps, x0, y0, x1, y1, &clumpinfo);


  /* Allocate the signal to noise table. */
  errno=0;
  sntab = *sntable = malloc( *numclumps * sizeof *sntab );
  if(sntab==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for sntab (clumps.c)",
          *numclumps*sizeof *sntab);


  /* Start calculating the Signal to noise ratios. */
  sntab[0]=0;
  for(i=1;i<*numclumps;++i)
    {
      /* These variables used for easy readability. */
      I = clumpinfo[ i * INFOTABCOLS     ];
      O = clumpinfo[ i * INFOTABCOLS + 2 ];

      /* If the inner flux is smaller than the outer flux (happens
	 only in noise cases) or the area is smaller than the minimum
	 area to calculate signal-to-noise, then set the S/N of this
	 segment to zero. */
      if( (Ni=clumpinfo[i*INFOTABCOLS+1]) > p->segsnminarea
          && I>O )              /* This is O, not 0 (zero). */
	{
          /* If the sky was subtracted then put in the standard
             deviation. */
          err  = ( p->skysubtracted
                   ? ( 2.0f * clumpinfo[i*INFOTABCOLS+4]
                       * clumpinfo[i*INFOTABCOLS+4] )
                   : 0.0f );

          /* Calculate the Signal to noise ratio. */
          ind = p->b0f1 ? i : counter++;
	  sntab[ind]=( sqrt((float)(Ni)/cpscorr)*(I-O)
		       / sqrt( (I>0?I:-1*I) + (O>0?O:-1*O) + err ) );
	}
      else
	sntab[i]=0;
    }

  /* If we are dealing with noise, replace the number of clumps with
     the number of those with a sufficient area and inner flux. */
  if(p->b0f1==0)
    *numclumps=counter;
}


















/******************************************************************/
/*************       S/N threshold on grid       ******************/
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
  float *sntable, minbfrac=p->minbfrac;
  char *segmentationname=p->segmentationname;
  size_t i, *indexs=&mp->indexs[mtp->id*mp->thrdcols];
  size_t s0, s1, nf, area, ind, size, startind, is1=mp->s1;

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
      count_f_b_onregion(p->byt, startind, s0, s1, is1, &nf, &area);
      if( (float)area < (float)(size)*minbfrac )
        {
          if(segmentationname)
            longinitonregion(clab, 0, startind, s0, s1, is1);
          continue;
        }


      /* We want to find the clumps on the noise, not the signal, the
         number of noise pixels in this mesh were calculated above
         (area). Here we want to pull out the indexs of those pixels
         which is necessary for the over-segmentation. */
      index_f_b_onregion(p->byt, startind, s0, s1, is1, inds, 0);


      /* Sort the indexs based on the flux within them. */
      qsort(inds, area, sizeof(size_t), indexfloatdecreasing);


      /* Do the over-segmentation: */
      oversegment(p, inds, area, startind/is1, startind%is1,
                  startind/is1+s0, startind%is1+s1, &numclumps);
      if(numclumps<p->minnumfalse) continue;

      /* Find the signal to noise of all the clumps. */
      clumpsntable(p, inds, area, &numclumps, startind/is1, startind%is1,
                   startind/is1+s0, startind%is1+s1, &sntable);


      /* Cleanup: */
      free(sntable);
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

  if(p->segmentationname)
    arraytofitsimg(p->segmentationname, "NoiseOversegmentaion",
                   LONG_IMG, p->clab, p->smp.s0, p->smp.s1, 0, p->wcs,
                   NULL, SPACK_STRING);
}
