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

#include "mesh.h"
#include "mode.h"
#include "forqsort.h"
#include "statistics.h"
#include "astrthreads.h"



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





void
checkcharray(struct meshparams *mp, int operationid,
             float **out1, float **out2)
{
  int needstwo=0;
  float *f, *fp, *ff, charray1, charray2;
  size_t i, row, start, *types=mp->types;
  size_t s0, s1, is1=mp->s1, *ts0=mp->ts0, *ts1=mp->ts1;

  if(operationid==MODEEQMED_AVESTD)
    needstwo=1;

  /* Allocate the array to keep the mesh indexs. Calloc is used so we
     can add all the indexes to the existing value to make sure that
     there is no overlap. */
  errno=0;
  *out1=malloc(mp->s0*mp->s1*sizeof **out1);
  if(*out1==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for out1 in checkcharray "
          "(mesh.c)", mp->s0*mp->s1*sizeof **out1);
  if(needstwo)
    {
      errno=0;
      *out2=malloc(mp->s0*mp->s1*sizeof **out2);
      if(*out2==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for out2 in checkcharray "
              "(mesh.c)", mp->s0*mp->s1*sizeof **out2);
    }

  /* Fill the indexs: */
  for(i=0;i<mp->nmeshi;++i)
    {
      row=0;
      s0=ts0[types[i]];
      s1=ts1[types[i]];
      start=mp->start[i];
      charray1=mp->charray1[i];
      if(needstwo) charray2=mp->charray2[i];
      do
        {
          fp= ( f = *out1 + start + row * is1 ) + s1;
          if(needstwo) ff= *out2 + start + row * is1;
          do
            {
              *f++ = charray1;
              if(needstwo) *ff++ = charray2;
            }
          while(f<fp);
          ++row;
        }
      while(row<s0);
    }
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

   charray:
   -------
   Or grid-array. It used for operations on the mesh grid, where one
   value is to be assigned for each mesh. It has one element for
   each mesh in the image.  Each channel has its own part of this
   larger array. The respective parts have gs0*gs1 elements. There
   are `nch' parts. in total.
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
  mp->charray1=mp->charray2=NULL;
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

  /* Distribute the meshes in all the threads. */
  distinthreads(mp->nmeshi, mp->numthreads, &mp->indexs, &mp->thrdcols);
}





void
freemesh(struct meshparams *mp)
{
  free(mp->start);
  free(mp->types);
  free(mp->chindex);
  free(mp->charray1);
  free(mp->charray2);
  free(mp->imgindex);
}




















/*********************************************************************/
/********************       Mesh operations       ********************/
/*********************************************************************/
void *
fillmeshonthreads(void *inparam)
{
  struct fillmeshparams *fcp=(struct fillmeshparams *)inparam;
  struct meshparams *mp=fcp->mp;

  float modesym;
  size_t modeindex;
  float sigcliptolerance=mp->sigcliptolerance;
  float *f, *mesh, *img, *imgend, *inimg=mp->img;
  size_t s0, s1, ind, row, num, start, is1=mp->s1;
  size_t i, *indexs=&mp->indexs[fcp->id*mp->thrdcols];
  float ave, med, std, sigclipmultip=mp->sigclipmultip;
  float mirrordist=mp->mirrordist, minmodeq=mp->minmodeq;


  /* Allocate the array that will keep this mesh's pixels
     contiguously. In fillmeshinfo the maximum s0 and s1 sizes of the
     arrays were stored so we could allocate one array for all the
     meshes irrespective of their type. */
  errno=0;
  mesh=fcp->alltypes=malloc(mp->maxs0*mp->maxs1*sizeof *fcp->alltypes);
  if(mesh==NULL)
    error(EXIT_FAILURE, errno, "Unable to allocate %lu bytes for"
          "fcp->alltypes of thread %lu in allocatetypearrays of mesh.c.",
          mp->maxs0*mp->maxs1*sizeof *fcp->alltypes, fcp->id);


  /* Start this thread's work: */
  for(i=0;indexs[i]!=NONTHRDINDEX;++i)
    {
      /* Prepare the values: */
      f=mesh;
      num=row=0;
      ind=indexs[i];
      start=mp->start[ind];
      s0=mp->ts0[mp->types[ind]];
      s1=mp->ts1[mp->types[ind]];

      /* Copy all the non-NaN pixels images pixels of this mesh into
         the mesh array. Note that currently, the spatial positioning
         of the pixels is irrelevant, so we only keep those that are
         non-NaN.*/
      do
        {
          imgend=(img = inimg + start + row++ * is1 ) + s1;
          do
            if(!isnan(*img))
            {
              ++num;
              *f++ = *img;
            }
          while(++img<imgend);
        }
      while(row<s0);

      /* Do the desired operation on the mesh: */
      qsort(mesh, num, sizeof *mesh, floatincreasing);
      modeindexinsorted(mesh, num, mirrordist, &modeindex, &modesym);
      if( modesym>MODESYMGOOD && (float)modeindex/(float)num>minmodeq )
        {
          switch(fcp->operationid)
            {

            case MODEEQMED_AVESTD: /* Average and standard devaition. */
              if(sigmaclip_converge(mesh, 1, num, sigclipmultip,
                                    sigcliptolerance, &ave, &med, &std))
                {mp->charray1[ind]=ave; mp->charray2[ind]=std;}
              else
                {mp->charray1[ind]=NAN; mp->charray2[ind]=NAN;}
              break;

            case MODEEQMED_QUANT: /* Quantile value.                 */
              mp->charray1[ind]=mesh[indexfromquantile(num, fcp->value)];
              break;

            default:
              error(EXIT_FAILURE, 0, "A bug! Please contact us at "
                    PACKAGE_BUGREPORT" so we can correct this issue. The "
                    "value to fcp->operationid is not recognized in "
                    "fillmeshonthreads (mesh.c).");
            }
        }
      else
        {
          mp->charray1[ind]=NAN;
          if(fcp->operationid==MODEEQMED_AVESTD)
            mp->charray2[ind]=NAN;
        }

    }

  /* Free alltype and if multiple threads were used, wait until all
     other threads finish. */
  free(fcp->alltypes);
  if(mp->numthreads>1)
    pthread_barrier_wait(&mp->b);
  return NULL;
}





void
fillmesh(struct meshparams *mp, int operationid, float value)
{
  int err;
  size_t i, nb;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  pthread_attr_t attr;
  struct fillmeshparams *fcp;
  size_t numthreads=mp->numthreads;

  /* Allocate the arrays to keep the thread and parameters for each
     thread. */
  errno=0; fcp=malloc(numthreads*sizeof *fcp);
  if(fcp==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes in fillmesh (mesh.c) for fcp",
          numthreads*sizeof *fcp);

  /* Allocate the array to keep the values for each  */
  errno=0; mp->charray1=malloc(mp->nmeshi*sizeof *mp->charray1);
  if(mp->charray1==NULL) error(EXIT_FAILURE, errno, "mp->charray1");
  if(operationid==MODEEQMED_AVESTD)
    {
      errno=0; mp->charray2=malloc(mp->nmeshi*sizeof *mp->charray2);
      if(mp->charray2==NULL) error(EXIT_FAILURE, errno, "mp->charray2");
    }

  /* Spin off the threads: */
  if(numthreads==1)
    {
      fcp[0].id=0;
      fcp[0].mp=mp;
      fcp[0].value=value;
      fcp[0].operationid=operationid;
      fillmeshonthreads(&fcp[0]);
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
            fcp[i].id=i;
	    fcp[i].mp=mp;
            fcp[i].value=value;
            fcp[i].operationid=operationid;
	    err=pthread_create(&t, &attr, fillmeshonthreads, &fcp[i]);
	    if(err) error(EXIT_FAILURE, 0, "Can't create thread %lu.", i);
	  }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&mp->b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&mp->b);
    }

  free(fcp);
}
