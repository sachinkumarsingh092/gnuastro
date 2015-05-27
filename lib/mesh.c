/*********************************************************************
meshgrid -- Create a mesh grid ready for multithreaded analysis.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include "mesh.h"
#include "mode.h"
#include "forqsort.h"
#include "neighbors.h"
#include "linkedlist.h"
#include "statistics.h"
#include "fitsarrayvv.h"
#include "astrthreads.h"
#include "spatialconvolve.h"



/*
   Basic idea
   ==========

   The meshes in an image are completely independent. There is no
   overlap. So when we want to do any sort of processing on separate
   meshes, they can be done completely independently. Whether they are
   positioned in the same channel or not. Usually, after the initial
   processing is complete, and some meshes pass the test and some
   don't, we need to do interpolation to find values for all the
   meshes and finally smooth the result.

   Maybe in the future all these arrays can be removed and calculated
   immediately for each mesh during the processing. This will increase
   the speed and decrease the amount of necessary memory. However, at
   this early stage, I am doing this job separately here for now so
   the main functions can be cleaner and more easily
   understandable. Once they have become mature enough, all these
   calculations should be done during the processing of each mesh.
*/









/*********************************************************************/
/********************     Preparing for Save      ********************/
/*********************************************************************/
/* Save the meshid of each pixel into an array the size of the image. */
void
checkmeshid(struct meshparams *mp, long **out)
{
  long i, *l, *lp;
  size_t row, start, *types=mp->types;
  size_t s0, s1, is1=mp->s1, *ts0=mp->ts0, *ts1=mp->ts1;

  /* Allocate the array to keep the mesh indexs. Calloc is used so we
     can add all the indexes to the existing value to make sure that
     there is no overlap. */
  errno=0;
  *out=calloc(mp->s0*mp->s1, sizeof **out);
  if(*out==NULL)
    error(EXIT_FAILURE, errno, "The array to show mesh labels in checkmesh.");

  /* Fill the indexs: */
  for(i=0;i<mp->nmeshi;++i)
    {
      row=0;
      s0=ts0[types[i]];
      s1=ts1[types[i]];
      start=mp->start[i];
      do
        {
          lp=(l = *out + start + row++ * is1 ) + s1;
          do *l++ += i; while(l<lp);
        }
      while(row<s0);
    }
}





/* Put the values of the check array(s) into an array the size of the
   input image. Note that the check arrays are only the size of the
   number of meshs, not the actual input image size. */
void
checkgarray(struct meshparams *mp, float **out1, float **out2)
{
  size_t gs0=mp->gs0, gs1=mp->gs1;
  float *f, *fp, *ff, garray1, garray2;
  size_t i, row, start, meshid, *types=mp->types;
  size_t f0, f1, fs1=mp->gs1*mp->nch1, chid, inchid;
  size_t s0, s1, is1=mp->s1, *ts0=mp->ts0, *ts1=mp->ts1;

  /* Allocate the array to keep the mesh indexs. Calloc is used so we
     can add all the indexes to the existing value to make sure that
     there is no overlap. */
  errno=0; *out1=malloc(mp->s0*mp->s1*sizeof **out1);
  if(*out1==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for out1 in checkgarray "
          "(mesh.c)", mp->s0*mp->s1*sizeof **out1);
  if(mp->cgarray2)
    {
      errno=0; *out2=malloc(mp->s0*mp->s1*sizeof **out2);
      if(*out2==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for out2 in checkgarray "
              "(mesh.c)", mp->s0*mp->s1*sizeof **out2);
    }

  /* Fill the array: */
  for(i=0;i<mp->nmeshi;++i)
    {
      /* Find the proper index for this mesh. garrays can either be
         indexed by each channel (so the indexs in each channel are
         contiguous) or by the image (so the mesh indexs are
         contiguous over the image). If mp->fgarray1==mp->garray1,
         then it is the latter. The last operation that worked on
         garrays, set garray1 either to cgarray1 or fgarray1. */
      if(mp->garray1==mp->fgarray1)
        {
          /* In this case, `i' is the index of a mesh in the full
             image, not the channel separated image. So we have to
             correct it to the channel independent indexing system. */
          f0=i/fs1;
          f1=i%fs1;
          inchid = (f0%gs0)*gs1      + f1%gs1;
          chid   = (f0/gs0)*mp->nch1 + f1/gs1;
          meshid = chid * mp->nmeshc + inchid;
        }
      else meshid=i;

      /* Fill the output array with the value in this mesh: */
      row=0;
      s0=ts0[types[meshid]];
      s1=ts1[types[meshid]];
      start=mp->start[meshid];
      garray1 = mp->garray1[i];
      if(mp->cgarray2)
        garray2 = mp->garray2[i];
      do
        {
          fp= ( f = *out1 + start + row * is1 ) + s1;
          if(mp->cgarray2) ff= *out2 + start + row * is1;
          do
            {
              *f++ = garray1;
              if(mp->cgarray2) *ff++ = garray2;
            }
          while(f<fp);
          ++row;
        }
      while(row<s0);
    }
}





/* By default, the garray1 and garray2 arrays keep the meshes of each
   channel contiguously. So in practice, each channel is like a small
   independent image. This will cause problems when we want to work on
   the meshs irrespective of which channel they belong to. This
   function allocates and fills in the fgarray1 and fgarray2 arrays.

   The explanation above is for the case when reverse==0. If it is set
   equal to 1 (or any non-zero number), then
*/
void
fullgarray(struct meshparams *mp, int reverse)
{
  size_t nch1=mp->nch1;
  size_t ind, gs1=mp->gs1, gs0=mp->gs0;
  float *fgarray1=NULL, *fgarray2=NULL;
  float *cgarray1=mp->cgarray1, *cgarray2=mp->cgarray2;
  size_t g0, g1, f0, f1, fmeshind, chid, is1=mp->nch1*mp->gs1;

  /* First allocate the fgarrays if they were not already full. */
  if(mp->fgarray1==NULL)
    {
      /* A simple sanity check */
      if(reverse)
        error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can "
              "fix this problem. For some reason, fullgarray has been called "
              "with the `reverse' flag set to true while fgarray is not "
              "allocated! This should not happen.", PACKAGE_BUGREPORT);

      /* Allocate the fgarrays */
      errno=0; fgarray1=mp->fgarray1=malloc(mp->nmeshi*sizeof *fgarray1);
      if(fgarray1==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for fgarray1",
              mp->nmeshi*sizeof *fgarray1);
      if(cgarray2)
        {
          errno=0;
          fgarray2=mp->fgarray2=malloc(mp->nmeshi*sizeof *fgarray2);
          if(fgarray2==NULL)
            error(EXIT_FAILURE, errno, "%lu bytes for fgarray2",
                  mp->nmeshi*sizeof *fgarray2);
        }
    }
  else
    {
      fgarray1=mp->fgarray1;
      fgarray2=mp->fgarray1;
    }

  /* Fill the fgarrays with the proper values: */
  for(chid=0;chid<mp->nch;++chid)
    {
      /* The first pixel ID in this channel is: chid*mp->nmeshi. */
      f0=(chid/nch1)*gs0;      /* Position of first channel mesh */
      f1=(chid%nch1)*gs1;      /* in full image.                 */
      fmeshind=chid*mp->nmeshc;

      /* Go over all the meshes in this channel and put them in their
         righful place. */
      for(ind=0;ind<mp->nmeshc;++ind)
        {
          g0=ind/gs1;    /* Position of this mesh in this channel. */
          g1=ind%gs1;
          if(reverse)
            cgarray1[ind+fmeshind]=fgarray1[(g0+f0)*is1+g1+f1];
          else
            fgarray1[(g0+f0)*is1+g1+f1]=cgarray1[ind+fmeshind];
          if(cgarray2)
            {
              if(reverse)
                cgarray2[ind+fmeshind]=fgarray2[(g0+f0)*is1+g1+f1];
              else
                fgarray2[(g0+f0)*is1+g1+f1]=cgarray2[ind+fmeshind];
            }
        }
    }

  /* Just for a check:
  arraytofitsimg("nochannels.fits", "fgarray1", FLOAT_IMG, fgarray1,
                 mp->nch2*mp->gs0, mp->nch1*mp->gs1, 1, NULL, NULL,
                 "mesh");
  if(cgarray2)
    arraytofitsimg("nochannels.fits", "fgarray2", FLOAT_IMG, fgarray2,
                   mp->nch2*mp->gs0, mp->nch1*mp->gs1, 1, NULL, NULL,
                   "mesh");
  */
}


















/*********************************************************************/
/********************      Creating the mesh      ********************/
/*********************************************************************/
/* In each channel we can have four types (sizes of meshs):

   0. If meshsize is chosen nicely in relation to the channel size,
      the majority (or ideally all) of the mesh tiles are going to
      have an area of meshsize*meshsize. These are built from pixel
      (0,0) until ((s0/mesh_size)*mesh_size,
      (s1/mesh_size)*mesh_size), note that these are integer
      divisions, not floating point, so 5/2=2.

   1. The top row of meshes. If s0%meshsize is larger than a certain
      fraction of meshsize, they will have s0%meshsize pixels in each
      row. If not, they will have meshsize+s0%meshsize rows. In any
      case, all but the last (type 3) will have meshsize columns. The
      number of columns of the last mesh in the top row, will depend
      on the third type of meshes.

   2. Just like type 2, but for the width or the number of columns in
      the image.

   3. The size of the last (top left) mesh in the grid will get its
      height (number of rows) from type 2s and its width from type 3s.

   Note that the numbers for the types used here are only for when all
   four types exist, when there is only one or two types, the first
   (in the order above) gets type id of zero and the second gets type
   id of 1. In the actual function, the indexs are given in the order
   the types are found.
 */
void
fillmeshinfo(struct meshparams *mp, size_t chs0, size_t chs1,
             size_t lasts0, size_t lasts1)
{
  size_t i, j, chi, chj, gs0=mp->gs0, gs1=mp->gs1;
  size_t meshsize=mp->meshsize, *chindex=mp->chindex;
  size_t numtypes=0, totalmeshcount=0, chid, s1=mp->s1;
  size_t meshid, typeind, lasti, lastj, nmeshc=mp->nmeshc;
  size_t nch1=mp->nch1, nch2=mp->nch2, *ts0=mp->ts0, *ts1=mp->ts1;
  size_t *types=mp->types, *start=mp->start, *imgindex=mp->imgindex;

  /* Initialize the maximum number of rows and columns: */
  mp->maxs1=mp->maxs0=0;

  /* Main meshs (type 0) in each channel. */
  if(gs0>1 && gs1>1)
    {
      typeind=numtypes++;
      ts1[typeind]=ts0[typeind]=meshsize;
      lasti = lasts0==meshsize ? gs0 : gs0-1;
      lastj = lasts1==meshsize ? gs1 : gs1-1;
      if(ts0[typeind]>mp->maxs0) mp->maxs0=ts0[typeind];
      if(ts1[typeind]>mp->maxs1) mp->maxs1=ts1[typeind];

      for(chi=0;chi<nch2;++chi)   /* Don't forget that C and FITS */
	for(chj=0;chj<nch1;++chj) /* axises have different order. */
          {
            chid=chi*nch1+chj;
            for(i=0;i<lasti;++i)
              for(j=0;j<lastj;++j)
                {
                  ++totalmeshcount;
                  meshid=chid*nmeshc + i*gs1+j;
                  types[meshid]=typeind;
                  chindex[meshid]=i*gs1+j;
                  imgindex[meshid]=( (chi*gs0+i) * nch1*gs1 + chj*gs1+j );
                  start[meshid]=( (chi*chs0+i*meshsize) * s1
                                  + (chj*chs1+j*meshsize) );
                }
          }
    }

  /* Top row of meshes (type 1) in each channel. */
  if(gs1>1 && lasts0!=meshsize)
    {
      typeind=numtypes++;
      ts0[typeind] = lasts0;
      ts1[typeind] = meshsize;
      lastj = lasts1==meshsize ? gs1 : gs1-1;
      if(ts0[typeind]>mp->maxs0) mp->maxs0=ts0[typeind];
      if(ts1[typeind]>mp->maxs1) mp->maxs1=ts1[typeind];

      for(chi=0;chi<nch2;++chi)
	for(chj=0;chj<nch1;++chj)
          {
            chid=chi*nch1+chj;
            for(j=0;j<lastj;++j)
              {
                ++totalmeshcount;
                meshid=chid*nmeshc + (gs0-1)*gs1+j;
                types[meshid]=typeind;
                chindex[meshid]=(gs0-1)*gs1+j;
                imgindex[meshid]=(chi*gs0+gs0-1) * nch1*gs1 + (chj*gs1+j);
                start[meshid] = ( (chi*chs0+(gs0-1)*meshsize) * s1
                                  + (chj*chs1+j*meshsize) );
              }
          }
    }

  /* Left column of meshes (type 2) in each channel. */
  if(gs0>1 && lasts1!=meshsize)
    {
      typeind=numtypes++;
      ts0[typeind] = meshsize;
      ts1[typeind] = lasts1;
      lasti = lasts0==meshsize ? gs0 : gs0-1;
      if(ts0[typeind]>mp->maxs0) mp->maxs0=ts0[typeind];
      if(ts1[typeind]>mp->maxs1) mp->maxs1=ts1[typeind];

      for(chi=0;chi<nch2;++chi)
	for(chj=0;chj<nch1;++chj)
          {
            chid=chi*nch1+chj;
            for(i=0;i<lasti;++i)
              {
                ++totalmeshcount;
                meshid=chid*nmeshc + i*gs1+(gs1-1);
                types[meshid]=typeind;
                chindex[meshid]=i*gs1+(gs1-1);
                imgindex[meshid]=( (chi*gs0+i) * nch1*gs1
                                   + (chj*gs1+(gs1-1)) );
                start[meshid] = ( (chi*chs0+i*meshsize) * s1
                                  + (chj*chs1+(gs1-1)*meshsize) );
              }
          }
    }

  /* Top left mesh or only mesh (type 3) in each channel. Note that
     this might only happen for one mesh in each channel (if there are
     4 types of meshes). */
  if(mp->nmeshi-totalmeshcount==mp->nch)
    {
      typeind=numtypes++;
      ts0[typeind] = lasts0;
      ts1[typeind] = lasts1;
      if(ts0[typeind]>mp->maxs0) mp->maxs0=ts0[typeind];
      if(ts1[typeind]>mp->maxs1) mp->maxs1=ts1[typeind];

      for(chi=0;chi<nch2;++chi)
	for(chj=0;chj<nch1;++chj)
	  {
            ++totalmeshcount;
	    meshid=(chi*nch1+chj)*nmeshc + gs0*gs1-1;
	    types[meshid]=typeind;
            chindex[meshid]=gs0*gs1-1;
	    imgindex[meshid]=( (chi*gs0+(gs0-1)) * nch1*gs1
                               + (chj*gs1+(gs1-1)) );
	    start[meshid]=( (chi*chs0+(gs0-1)*meshsize) * s1
				 + (chj*chs1+(gs1-1)*meshsize) );
	  }
    }

  /* Just for a check: */
  if(totalmeshcount!=mp->nmeshi)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix it."
          "The basic information for some meshes has not been found (in "
          "fillmeshinfo of mesh.c).", PACKAGE_BUGREPORT);
}





/* In the explanation below, all the parameters named are from the
   `struct meshparams` of `mesh.h`.

   An image contians `s0*s1` pixels in s0 rows and s1 columns. A mesh
   is a box of pixels within the image. A channel is defined as large
   regions of the image that should contain the meshes discussed
   above. So each channel contains a certain number of meshs and each
   mesh contains a certain number of pixels.

   In general, for any channel size, there can be four types of
   meshes. See the explanations above fillmeshinfo() for a more
   detailed discussion. So when all the channels have the same size,
   in the whole image, there can only be four types of meshes.

   The idea behind this classification is that CCDs often have more
   than one reading channel and essentially all the properties, vary
   from channel to channel. Therefore it is very important to do
   convolution, thresholding, detection and finding the sky background
   and noise properties differently for each channel.

   Each mesh is treated independently in its thread, but when
   interpolation over the meshes is desired, only meshes in each
   channel should be interpolated over.
*/
void
makemesh(struct meshparams *mp)
{
  size_t meshsize=mp->meshsize, lasts0, lasts1;
  size_t i, chs0=mp->s0/mp->nch2, chs1=mp->s1/mp->nch1;

  /* Set the number of channels. */
  mp->nch=mp->nch1*mp->nch2;

  /* Incase meshsize is larger than the channel sizes, make it equal
     to the smaller of the channel axis sizes. */
  if(meshsize>chs0 || meshsize>chs1)
    meshsize = (chs0<chs1
                ? (chs0%2 ? chs0-1 : chs0)
                : (chs1%2 ? chs1-1 : chs1));

  /* Initialize the necessary arrays: */
  for(i=0;i<4;++i) mp->ts0[i]=mp->ts1[i]=0;

  /* Set the mesh grid height (gs0, number of rows) and width (gs1,
     number of columns) on each channel. Note that the channels are
     assumed to be identical in size, so these values are the same for
     all the channels. In case the remainder of the image side size
     and meshsize is larger than lastmeshfrac, add a new mesh. If it
     is smaller than that fraction, then the remainder should be
     merged with the last mesh in the row or column. */
  if( (float)(chs0%meshsize) > mp->lastmeshfrac*(float)(meshsize) )
    {
      mp->gs0 = chs0/meshsize + 1;
      lasts0 = chs0%meshsize;
    }
  else
    {
      mp->gs0 = chs0/meshsize;
      lasts0 = meshsize+chs0%meshsize;
    }
  if( (float)(chs1%meshsize) > mp->lastmeshfrac*(float)(meshsize) )
    {
      mp->gs1 = chs1/meshsize + 1;
      lasts1 = chs1%meshsize;
    }
  else
    {
      mp->gs1 = chs1/meshsize;
      lasts1 = meshsize+chs1%meshsize;
    }
  mp->nmeshc=mp->gs0*mp->gs1;
  mp->nmeshi=mp->nmeshc*mp->nch;

  /* Allocate the arrays to keep all the mesh starting points and
     types. Irrespective of which channel they lie in. */
  mp->cgarray1=mp->cgarray2=mp->fgarray1=mp->fgarray2=NULL;
  errno=0; mp->start=malloc(mp->nmeshi*sizeof *mp->start);
  if(mp->start==NULL) error(EXIT_FAILURE, errno, "Mesh starting points");
  errno=0; mp->types=malloc(mp->nmeshi*sizeof *mp->types);
  if(mp->types==NULL) error(EXIT_FAILURE, errno, "Mesh types");
  errno=0; mp->chindex=malloc(mp->nmeshi*sizeof *mp->chindex);
  if(mp->chindex==NULL) error(EXIT_FAILURE, errno, "Mesh in channel index");
  errno=0; mp->imgindex=malloc(mp->nmeshi*sizeof *mp->imgindex);
  if(mp->imgindex==NULL) error(EXIT_FAILURE, errno, "Mesh in image index");

  /* Fill in the information for each mesh and each type. */
  fillmeshinfo(mp, chs0, chs1, lasts0, lasts1);
}





void
freemesh(struct meshparams *mp)
{
  free(mp->start);
  free(mp->types);
  free(mp->chindex);
  free(mp->cgarray1);
  free(mp->cgarray2);
  free(mp->fgarray1);
  free(mp->fgarray2);
  free(mp->imgindex);
}




















/*********************************************************************/
/********************       Mesh operations       ********************/
/*********************************************************************/
/* The purpose of this function is to prepare the (possibly)
   multi-threaded environemnt, spin-off each thread and wait for them
   to finish for any operation that is to be done on the mesh
   grid.

   The arguments are:

   1. A pointer to the meshparams structure that keeps all the
      information.

   2. A pointer to a function that returns and gets a `void *' as its
      only argument. This function will be directly given to
      pthread_create. Through this function, you can any function that
      you wish to operate on the mesh grid with.

   3. The size of each element to copy the mesh grid into, this has to
      be type size of the same type that constitutes `img' in
      meshparams. If the value to this argument is non-zero, an array
      will be allocated that can contain all the pixels in all the
      meshs and can be used by threads to manipute the pixels (for
      example sort them) in each mesh.

   4. If the value of this argument is 1, then a second garray will be
      allocated in case your operation needs one.
*/
void
operateonmesh(struct meshparams *mp, void *(*meshfunc)(void *),
              size_t oneforallsize, int makegarray2)
{
  int err;
  size_t i, nb;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  pthread_attr_t attr;
  struct meshthreadparams *mtp;
  size_t numthreads=mp->numthreads;

  /* Allocate the arrays to keep the thread and parameters for each
     thread. */
  errno=0; mtp=malloc(numthreads*sizeof *mtp);
  if(mtp==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes in fillmesh (mesh.c) for mtp",
          numthreads*sizeof *mtp);

  /* Allocate the array to keep the values for each  */
  errno=0; mp->cgarray1=malloc(mp->nmeshi*sizeof *mp->cgarray1);
  if(mp->cgarray1==NULL) error(EXIT_FAILURE, errno, "mp->cgarray1");
  if(makegarray2)
    {
      errno=0; mp->cgarray2=malloc(mp->nmeshi*sizeof *mp->cgarray2);
      if(mp->cgarray2==NULL) error(EXIT_FAILURE, errno, "mp->cgarray2");
    }
  else mp->cgarray2=NULL;
  mp->garray1=mp->cgarray1;
  mp->garray2=mp->cgarray2;

  /* `oneforall' is an array with the sides equal to the maximum side
     of the meshes in the image. The purpose is to enable manipulating
     the mesh pixels (for example sorting them and so on.). One such
     array is allocated for each thread in this one full
     allocation. Each thread can then use its own portion of this
     array through the following declaration:

     float *oneforall=&mp->oneforall[mtp->id*mp->maxs0*mp->maxs1];

     In meshparams, `oneforall' is defined as a `void *', so the
     caller function, can cast it to any type it wants. The size of
     each type is given to `fillmesh' through the `oneforallsize'
     argument.
  */
  if(oneforallsize)
    {
      errno=0;
      mp->oneforall=malloc(numthreads*mp->maxs0*mp->maxs1*oneforallsize);
      if(mp->oneforall==NULL)
        error(EXIT_FAILURE, errno, "Unable to allocate %lu bytes for"
              "mtp->oneforall in fillmesh of mesh.c.",
              numthreads*mp->maxs0*mp->maxs1*oneforallsize);
    }

  /* Distribute the meshes in all the threads. */
  distinthreads(mp->nmeshi, mp->numthreads, &mp->indexs, &mp->thrdcols);

  /* Spin off the threads: */
  if(numthreads==1)
    {
      mtp[0].id=0;
      mtp[0].mp=mp;
      (*meshfunc)(&mtp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
	 (that spinns off the nt threads) is also a thread, so the
	 number the barrier should be one more than the number of
	 threads spinned off. */
      if(mp->nmeshi<numthreads) nb=mp->nmeshi+1;
      else                      nb=numthreads+1;
      attrbarrierinit(&attr, &mp->b, nb);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
	if(mp->indexs[i*mp->thrdcols]!=NONTHRDINDEX)
	  {
            mtp[i].id=i;
	    mtp[i].mp=mp;
	    err=pthread_create(&t, &attr, meshfunc, &mtp[i]);
	    if(err) error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
	  }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&mp->b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&mp->b);
    }

  free(mtp);
  free(mp->indexs);
  if(oneforallsize) free(mp->oneforall);
}




















/*********************************************************************/
/********************         Interpolate         ********************/
/*********************************************************************/
/* Ideally, this should be a radial distance using the square root of
   the differences. But here we are not dealing with subpixel
   distances, so manhattan distance is fine. It also avoids the square
   root function which is much slower. */
float
manhattandistance(long ind, long xc, long yc, long s1)
{
  return labs(ind%s1 - xc) + labs(ind/s1 - yc);
}




/* Some of the variables have different names than the meshparams
   structure because they are to be fed into the FILL_NGB_4_ALLIMG
   macro. */
void *
meshinterponthread(void *inparams)
{
  struct meshthreadparams *mtp=(struct meshthreadparams *)inparams;
  struct meshparams *mp=mtp->mp;

  /* Basic variables used in other definitions: */
  size_t numnearest=mp->numnearest;
  size_t is0=mp->fullinterpolation ? mp->gs0*mp->nch2 : mp->gs0;
  size_t is1=mp->fullinterpolation ? mp->gs1*mp->nch1 : mp->gs1;

  /* Variables for this function: */
  struct tosll *lQ, *sQ;
  size_t xc, yc, *n, *nf, currentnum, thisind;
  unsigned char *byt=&mp->byt[mtp->id*is0*is1];
  float *nearest1=&mp->nearest1[mtp->id*numnearest];
  size_t i, *indexs=&mp->indexs[mtp->id*mp->thrdcols];
  float mdist, *garray1=mp->garray1, *garray2=mp->garray2;
  size_t fmeshid, position, *ind=&position, numngb, ngb[4];
  float *outgarray1=mp->outgarray1, *outgarray2=mp->outgarray2;
  float *nearest2=garray2 ? &mp->nearest2[mtp->id*numnearest] : NULL;

  /* Go over all the meshes for this thread. */
  for(i=0;indexs[i]!=NONTHRDINDEX;++i)
    {
      /* Reset the byt array: */
      memset(byt, 0, is0*is1);

      /* Get the index of this NaN mesh: */
      thisind=mp->naninds[indexs[i]];
      /*printf("Started mesh %lu\n", thisind);*/

      /*
         `fmeshid' (first-mesh-id) keeps the index of the first mesh
         in this channel. We are going to consider each channel as a
         sparate image for the neighbor finding job. This is why the
         pixels in each channel were designed to be contiguous. When
         we want to look at the value on each mesh, we simply sum its
         channel index with fmeshid. From this point on, *ind
         corresponds to the index of the mesh in the channel.

         We will also only be concerened with the portion of `mp->byt'
         array that is within this channel. So the `byt' pointer is
         set to point to the first mesh index in this channel.
      */
      lQ=sQ=NULL;
      fmeshid = ( mp->fullinterpolation
                  ? 0
                  : (thisind/mp->nmeshc) * mp->nmeshc);

      /* Initialize the necessary parameters. */
      *ind=thisind-fmeshid;
      xc=*ind%is1;
      yc=*ind/is1;
      byt[*ind]=1;
      currentnum=0;
      add_to_tosll_end( &lQ, &sQ, *ind, 0 );

      /* Start finding the nearest filled pixels. */
      while(sQ)
        {
          /* Pop out a pixel index (p) from the queue: */
          pop_from_tosll_start(&lQ, &sQ, ind, &mdist);

          /* If it isn't a NaN, then put it in the `nearest1' and
             `nearest2' arrays. */
          if(!isnan(garray1[*ind+fmeshid]))
            {
              nearest1[currentnum]=garray1[*ind+fmeshid];
              if(garray2) nearest2[currentnum]=garray2[*ind+fmeshid];
              if(++currentnum>=numnearest) break;
            }

          /* Check the four neighbors and if they have not already
             been checked, put them into the queue. */
          FILL_NGB_4_ALLIMG;
          nf=(n=ngb)+numngb;
          do
            if(byt[*n]==0)
              {
                byt[*n]=1;
                add_to_tosll_end( &lQ, &sQ, *n,
                                  manhattandistance(*n, xc, yc, is1) );
              }
          while(++n<nf);

          /* If there are no more meshes to add to the queue, then
             this shows, there were not enough points for
             interpolation. Normally, this loop should only be exited
             through the `currentnum>=numnearest' check above. */
          if(sQ==NULL)
            error(EXIT_FAILURE, 0, "Only %lu mesh(s) are filled for "
                  "interpolation in channel %lu. Either set less "
                  "restrictive requirements to get more interpolation "
                  "points or decrease the number of nearest points to "
                  "use for interpolation. Problem encountered on thread "
                  "%lu, for pixel %lu.\n", currentnum, thisind/mp->nmeshc,
                  mtp->id, thisind);
        }
      tosll_free(lQ);       /* The rest of the queue not needed. */


      /* Find the median of the nearest neighbors and put it in: */
      qsort(nearest1, numnearest, sizeof *nearest1, floatincreasing);
      outgarray1[thisind] = ( numnearest%2 ?
                              nearest1[numnearest/2] : /* Odd.  */
                              (nearest1[numnearest/2]  /* Even. */
                               +nearest1[numnearest/2-1])/2 );
      /*
      if(outgarray1[thisind]<4000 || outgarray1[thisind]>5000)
        {
          printf("%lu", thisind);
          for(i=0;i<numnearest;++i)
            printf("\t%f\n", nearest1[i]);
          printf("Result: %f\n\n\n", outgarray1[thisind]);
        }
      */
      if(garray2)
        {
          qsort(nearest2, numnearest, sizeof *nearest2, floatincreasing);
          outgarray2[thisind] = ( numnearest%2 ?
                                  nearest2[numnearest/2] : /* Odd.  */
                                  (nearest2[numnearest/2]  /* Even. */
                                   +nearest2[numnearest/2-1])/2 );
        }
    }

  /* If there is more than one thread, wait until the others
     finish. */
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}




void
preparemeshinterparrays(struct meshparams *mp)
{
  float *garray1, *garray2;
  size_t numthreads=mp->numthreads;
  size_t numnan=0, bs0=mp->gs0, bs1=mp->gs1;
  float *f, *ff, *fff=NULL, *ffff=NULL, *fp;


  /* If full-interpolation is to be done (ignoring channels) we will
     need two arrays to temporarily keep the actual positions of the
     meshs. Note that by default, the meshs in each channels are
     placed contiguously. In these arrays, the meshs are placed as
     they would in the image (channels are ignored).*/
  if(mp->fullinterpolation)
    {
      if(mp->fgarray1==NULL)
        fullgarray(mp, 0);
      bs0=mp->nch2*mp->gs0;
      bs1=mp->nch1*mp->gs1;
      mp->garray1=mp->fgarray1;
      mp->garray2=mp->fgarray2;
    }
  else
    {
      mp->garray1=mp->cgarray1;
      mp->garray2=mp->cgarray2;
    }
  garray1=mp->garray1;
  garray2=mp->garray2;


  /* Allocate the output arrays to keep the final values. Note that we
     cannot just do the interpolaton on the same input grid, because
     the newly filled interpolated pixels will affect the later
     ones. In the end, these copies are going to replace the
     garrays. */
  errno=0; ff=mp->outgarray1=malloc(mp->nmeshi*sizeof *mp->outgarray1);
  if(mp->outgarray1==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for outgarray1",
          mp->nmeshi*sizeof *mp->outgarray1);
  if(garray2)
    {
      fff=garray2;
      errno=0; ffff=mp->outgarray2=malloc(mp->nmeshi*sizeof *mp->outgarray2);
      if(mp->outgarray2==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for outgarray2",
              mp->nmeshi*sizeof *mp->outgarray2);
    }



  /* Find how many blank pixels there are and keep all their indexs in
     the naninds array. Note that all garrays are going to have the
     same number of NaN meshs. */
  fp=(f=garray1)+mp->nmeshi;
  do if(isnan(*f++)) ++numnan; while(f<fp);
  errno=0; mp->naninds=malloc(numnan*sizeof *mp->naninds);
  if(mp->naninds==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for blank pixel indexs in "
          "meshinterpolate (mesh.c).", numnan*sizeof *mp->naninds);
  numnan=0; fp=(f=garray1)+mp->nmeshi;
  do
    {
      if(isnan(*f))
        mp->naninds[numnan++]=f-garray1;
      else
        {
          *ff=*f;
          if(garray2)
            *ffff=*fff;
        }
      if(garray2)
        {
          ++fff;
          ++ffff;
        }
      ++ff;
    }
  while(++f<fp);
  mp->numnan=numnan;


  /* Allocate the array for the values of the nearest pixels. For each
     pixel in each thread, the nearest pixels will be put here and
     sorted to find the median value. If garray2 is not NULL, then it
     should be interpolated. Note that both arrays have blank pixels
     on the same places.*/
  errno=0;
  mp->nearest1=malloc(numthreads*mp->numnearest*sizeof *mp->nearest1);
  if(mp->nearest1==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for the array to keep the "
          "nearest1 values for interpolation (mesh.c).",
          numthreads*mp->numnearest*sizeof *mp->nearest1);
  if(garray2)
    {
      errno=0;
      mp->nearest2=malloc(numthreads*mp->numnearest*sizeof *mp->nearest2);
      if(mp->nearest2==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for the array to keep the "
              "nearest2 values for interpolation (mesh.c).",
              numthreads*mp->numnearest*sizeof *mp->nearest2);
    }


  /* Allocate space for the byte array keeping a record of which
     pixels were used to search. Note that each thread needs one byt
     array. */
  errno=0; mp->byt=malloc(numthreads*bs0*bs1*sizeof *mp->byt);
  if(mp->byt==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for mp->byt array in "
          "mesh.c", numthreads*bs0*bs1*sizeof *mp->byt);
}





void
meshinterpolate(struct meshparams *mp)
{
  int err;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  size_t i, nb;
  pthread_attr_t attr;
  struct meshthreadparams *mtp;
  size_t numthreads=mp->numthreads;

  /* Prepare all the meshparams arrays: */
  preparemeshinterparrays(mp);

  /* Allocate the arrays to keep the thread and parameters for each
     thread. */
  errno=0; mtp=malloc(numthreads*sizeof *mtp);
  if(mtp==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes in fillmesh (mesh.c) for mtp",
          numthreads*sizeof *mtp);

  /* Distribute the blank pixels between the threads: */
  distinthreads(mp->numnan, numthreads, &mp->indexs, &mp->thrdcols);

  /* Spin off the threads: */
  if(numthreads==1)
    {
      mtp[0].id=0;
      mtp[0].mp=mp;
      meshinterponthread(&mtp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
	 (that spinns off the nt threads) is also a thread, so the
	 number the barrier should be one more than the number of
	 threads spinned off. */
      if(mp->numnan<numthreads) nb=mp->numnan+1;
      else                      nb=numthreads+1;
      attrbarrierinit(&attr, &mp->b, nb);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
	if(mp->indexs[i*mp->thrdcols]!=NONTHRDINDEX)
	  {
            mtp[i].id=i;
	    mtp[i].mp=mp;
	    err=pthread_create(&t, &attr, meshinterponthread, &mtp[i]);
	    if(err) error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
	  }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&mp->b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&mp->b);
    }

  /* Replace garray1 and garray2 with outgarray1 and outgarray2 */
  free(mp->garray1);
  if(mp->fullinterpolation) mp->garray1=mp->fgarray1=mp->outgarray1;
  else                      mp->garray1=mp->cgarray1=mp->outgarray1;
  if(mp->cgarray2)
    {
      free(mp->garray2);
      if(mp->fullinterpolation) mp->garray2=mp->fgarray2=mp->outgarray2;
      else                      mp->garray2=mp->cgarray2=mp->outgarray2;
    }
  mp->outgarray1=mp->outgarray2=NULL;

  /* For a check
  system("rm test.fits");
  arraytofitsimg("test.fits", "garray1", FLOAT_IMG, mp->garray1,
                 mp->nch2*mp->gs0, mp->nch1*mp->gs1, 1, NULL, NULL,
                 "mesh");
  arraytofitsimg("test.fits", "garray2", FLOAT_IMG, mp->garray2,
                 mp->nch2*mp->gs0, mp->nch1*mp->gs1, 1, NULL, NULL,
                 "mesh");
  */

  /* Clean up. */
  free(mtp);
  free(mp->byt);
  free(mp->indexs);
  free(mp->naninds);
  free(mp->nearest1);
  free(mp->nearest2);
}




















/*********************************************************************/
/********************           Smooth            ********************/
/*********************************************************************/
void
meshsmooth(struct meshparams *mp)
{
  float *charray;
  float *f, *o, *fp, *tmp, *kernel, *sgarray1, *sgarray2;
  size_t numthreads=mp->numthreads, gs0=mp->gs0, gs1=mp->gs1;
  size_t chid, nmeshc=mp->nmeshc, smoothwidth=mp->smoothwidth;

  /* Make the smoothing kernel and set all its elements to 1. */
  errno=0;
  kernel=malloc(smoothwidth*smoothwidth*sizeof *kernel);
  if(kernel==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for kernel in mesh.c",
          mp->smoothwidth*sizeof *kernel);
  fp=(f=kernel)+smoothwidth*smoothwidth; do *f++=1; while(f<fp);


  /* Smooth the garrays. */
  if(mp->fullsmooth)
    {
      if(mp->fgarray1==NULL || mp->fgarray1!=mp->garray1)
        fullgarray(mp, 0);

      spatialconvolve(mp->fgarray1, gs0*mp->nch2, gs1*mp->nch1, kernel,
                      smoothwidth, smoothwidth, numthreads, 1, &sgarray1);
      free(mp->fgarray1);
      mp->garray1=mp->fgarray1=sgarray1;
      if(mp->cgarray2)
        {
          spatialconvolve(mp->fgarray2, gs0*mp->nch2, gs1*mp->nch1, kernel,
                          smoothwidth, smoothwidth, mp->numthreads, 1,
                          &sgarray2);
          free(mp->fgarray2);
          mp->garray2=mp->fgarray2=sgarray2;
        }
    }
  else
    for(chid=0;chid<mp->nch;++chid) /* Note that in this mode, we do not */
      {                             /* Allocate anything :-).            */

        /* If the last operation was done on the full image, then
           mp->garray1==mp->fgarray1. Therefore, the mesh boxes in
           each channel will not be congituous. So we have to update
           cgarray and set mp->garray1=mp->cgarray1. */
        if(mp->garray1==mp->fgarray1)
          fullgarray(mp, 1);
        mp->garray1=mp->cgarray1;
        mp->garray2=mp->cgarray2;

        charray=&mp->cgarray1[chid*nmeshc];
        spatialconvolve(charray, gs0, gs1, kernel, smoothwidth, smoothwidth,
                        numthreads, 1, &tmp);
        o=tmp; fp=(f=charray)+gs0*gs1; do *f=*o++; while(++f<fp);
        free(tmp);
        if(mp->cgarray2)
          {
            charray=&mp->cgarray2[chid*nmeshc];
            spatialconvolve(charray, gs0, gs1, kernel, smoothwidth,
                            smoothwidth, numthreads, 1, &tmp);
            o=tmp; fp=(f=charray)+gs0*gs1; do *f=*o++; while(++f<fp);
            free(tmp);
          }
      }

  /* Clean up: */
  free(kernel);
}




















/*********************************************************************/
/********************      Convolve in mesh       ********************/
/*********************************************************************/
void*
meshspatialconvonthreads(void *inparam)
{
  struct meshthreadparams *mtp = (struct meshthreadparams *) inparam;
  struct meshparams *mp = mtp->mp;
  const size_t ks0=mp->ks0, ks1=mp->ks1, is1=mp->s1;

  double sum, ksum;
  const size_t *types=mp->types;
  size_t chid, *chbrd=mtp->chbrd, ystart;
  const size_t *ts0=mp->ts0, *ts1=mp->ts1;
  size_t i, meshid, ind, *start=mp->start;
  size_t x, y, fx, fy, a, b, k0h=ks0/2, k1h=ks1/2;
  size_t xmin, xmax, ymin, ymax, nmeshc=mp->nmeshc;
  size_t *indexs=&mp->indexs[mtp->id*mp->thrdcols];
  float *img=mp->img, *conv=mtp->conv, *kernel=mp->kernel;

  /* For each mesh in this thread, convolve the mesh */
  for(i=0; (meshid=indexs[i])!=NONTHRDINDEX; ++i)
    {
      /* Get the channel ID and the relevant information for this mesh. */
      chid=meshid/nmeshc;
      xmin=chbrd[chid*4]+k0h;    /* Within these four points, the      */
      xmax=chbrd[chid*4+2]-k0h;  /* convolution does not have to worry */
      ymin=chbrd[chid*4+1]+k1h;  /* about the edges. */
      ymax=chbrd[chid*4+3]-k1h;

      /* Get the starting and ending positions of this mesh. NOTE:
         ystart is needed because for every `x', the value of y needs
         to change. */
      fx = ( x      = start[meshid]/is1 ) + ts0[types[meshid]];
      fy = ( ystart = start[meshid]%is1 ) + ts1[types[meshid]];

      /* If the mesh is on the edge of the channel it should be
         treated differently compared to when it is not. */
      if( x>=xmin && ystart>=ymin && fx<xmax && fy<ymax )
	{                       /* All pixels in this mesh are distant  */
          for(;x<fx;++x)        /* enough from the edge of the channel. */
            for(y=ystart;y<fy;++y)
              {
                if(isnan(img[x*is1+y]))
                  conv[x*is1+y]=NAN;
                else
                  {
                    ksum=sum=0.0f;
                    for(a=0;a<ks0;++a)
                      for(b=0;b<ks1;++b)
                        if( !isnan( img[ ind = (x+a-k0h) * is1 + y+b-k1h ] ) )
                          {
                            ksum+=kernel[a*ks1+b];
                            sum+=img[ind] * kernel[a*ks1+b];
                          }
                    conv[x*is1+y] = sum/ksum;
                  }
              }
	}
      else
	{                       /* Some pixels in this mesh are too  */
	  for(;x<fx;++x)        /* close to the edge.                */
	    for(y=ystart;y<fy;++y)
	      {
                if( isnan(img[x*is1+y]) )
                    conv[x*is1+y]=NAN;
                else
                  {
                    ksum=sum=0.0f;
                    if(x>=xmin && y>=ymin && fx<xmax && fy<ymax)
                      for(a=0;a<ks0;++a)
                        for(b=0;b<ks1;++b)
                          {     /* {} needed, because of next `else'. */
                            if( !isnan(img[ ind=(x+a-k0h) * is1 + y+b-k1h ]))
                              {
                                ksum+=kernel[a*ks1+b];
                                sum+=img[ind] * kernel[a*ks1+b];
                              }
                          }
                    else
                      {
                        for(a=0;a<ks0;++a)
                          if(x+a>=xmin && x+a<chbrd[chid*4+2]+k0h)
                            for(b=0;b<ks1;++b)
                              if(y+b>=ymin && y+b<chbrd[chid*4+3]+k1h)
                                {
                                  if( !isnan( img[ind=(x+a-k0h)
                                                  * is1
                                                  + y+b-k1h]) )
                                    {
                                      ksum+=kernel[a*ks1+b];
                                      sum+=img[ind] * kernel[a*ks1+b];
                                    }
                                }
                      }
                    conv[x*is1+y] = sum/ksum;
                  }
            }
	}
    }

  /* Free alltype and if multiple threads were used, wait until all
     other threads finish. */
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}





void
spatialconvolveonmesh(struct meshparams *mp, float **conv)
{
  int err;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  pthread_attr_t attr;
  size_t i, nb, *chbrd;
  struct meshthreadparams *mtp;
  size_t numthreads=mp->numthreads;


  /* Allocate the arrays to keep the thread and parameters for each
     thread. */
  errno=0; mtp=malloc(numthreads*sizeof *mtp);
  if(mtp==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes in fillmesh (mesh.c) for mtp",
          numthreads*sizeof *mtp);


  /* Allocate space for the convolved array. */
  errno=0; *conv=malloc(mp->s0*mp->s1* sizeof **conv);
  if(*conv==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for convolution on mesh output.",
          mp->s0*mp->s1* sizeof **conv);


  /* Channel borders so we can see if a mesh is inside or outside a
     channel. chbrd contians for number for each channel: in order
     they are:

     1. The channel's first pixel's first C axis value.
     2. The channel's first pixel's second C axis value.
     3. The channel's last pixel's first C axis value.
     4. The channel's last pixel's second C axis value.
  */
  errno=0; chbrd=malloc(mp->nch*4*sizeof *chbrd);
  if(chbrd==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for chbrd in "
          "spatialconvolveonmesh", mp->nch*4*sizeof *chbrd);
  for(i=0;i<mp->nch;++i)
    {
      if(mp->fullconvolution)
        {
          chbrd[i*4+2]=mp->s0;
          chbrd[i*4+3]=mp->s1;
          chbrd[i*4]=chbrd[i*4+1]=0;
        }
      else
        {
          chbrd[i*4+0]=mp->start[i*mp->nmeshc]/mp->s1;
          chbrd[i*4+1]=mp->start[i*mp->nmeshc]%mp->s1;
          chbrd[i*4+2]=chbrd[i*4+0]+mp->s0/mp->nch2;
          chbrd[i*4+3]=chbrd[i*4+1]+mp->s1/mp->nch1;
        }
    }


  /* Distribute the meshes in all the threads. */
  distinthreads(mp->nmeshi, mp->numthreads, &mp->indexs, &mp->thrdcols);


  /* Spin off the threads: */
  if(numthreads==1)
    {
      mtp[0].id=0;
      mtp[0].mp=mp;
      mtp[0].conv=*conv;
      mtp[0].chbrd=chbrd;
      meshspatialconvonthreads(&mtp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
	 (that spinns off the nt threads) is also a thread, so the
	 number the barrier should be one more than the number of
	 threads spinned off. */
      if(mp->nmeshi<numthreads) nb=mp->nmeshi+1;
      else                      nb=numthreads+1;
      attrbarrierinit(&attr, &mp->b, nb);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
	if(mp->indexs[i*mp->thrdcols]!=NONTHRDINDEX)
	  {
            mtp[i].id=i;
	    mtp[i].mp=mp;
            mtp[i].conv=*conv;
            mtp[i].chbrd=chbrd;
	    err=pthread_create(&t, &attr, meshspatialconvonthreads, &mtp[i]);
	    if(err) error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
	  }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&mp->b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&mp->b);
    }

  free(mtp);
  free(chbrd);
  free(mp->indexs);
}
