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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/fits.h>
#include <gnuastro/mesh.h>
#include <gnuastro/mode.h>
#include <gnuastro/qsort.h>
#include <gnuastro/neighbors.h>
#include <gnuastro/linkedlist.h>
#include <gnuastro/statistics.h>
#include <gnuastro/spatialconvolve.h>



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

   The programs outside of mesh.c, should only be using garray (1 or
   2). Not cgarray or fgarray. The reason for this is that the user
   might choose to ignore channels on one step and use them in
   another. All the functions in mesh.c will set garray and use
   cgarray and fgarray such that when they return, garray points to
   their output. So the caller to any of these functions doesn't have
   to worry about which one was last used, they can just ignorantly
   use garray and everything will be fine.

   Maybe in the future all these arrays can be removed and calculated
   immediately for each mesh during the processing. This might
   increase the speed and decrease the amount of necessary
   memory. However, it will add to the complexity of the code. So at
   this early stage, I am doing this job separately here for now so
   the main functions can be cleaner and more easily
   understandable. Once they have become mature enough, all these
   calculations should be done during the processing of each mesh.
*/









/*********************************************************************/
/************       Finding the proper mesh ID      ******************/
/*********************************************************************/
/* There are two interal ways to store the IDs of the meshs:

   1. By default, the mesh IDs are set based on their position within
      the channels (so the meshs in each channel are contiguous). We
      call this the channel-based mesh ID. The cgarray1 and cgarray2
      store the mesh values based on this ID. All the basic mesh
      properties, like their types and height and widths are stored
      based on this ID.

   2. If the user asks for a full image interpolation or smoothing,
      then two new arrays are created called fgarray1 and
      fgarray2. These arrays keep the mesh values as they appear on
      the image, irrespective of which channel they belong to. We call
      this the image-based ID. The image based ID is contiguous over
      the image.

  There is a third pointer called garray1 and garray2. These can point
  to either the cgarrays or the fgarrays. Ideally, the user should be
  completely ignorant to what garray is pointing to. The caller knows
  the type of IDs they are using, but they don't know what to put into
  garray. We call the ID that should be put into garray the
  garray-based ID. It can be either of the two kinds above.

  So we have made the following two functions:

    gal_mesh_ch_based_id_from_gid:

       This is useful for when you are going over the elements in
       garray (and you are completley ignorant to which one of
       cgarrays or fgarrays garray points to) and you need the channel
       based IDs to get basic mesh information like the mesh type and
       size.

    gal_mesh_gid_from_ch_based_id:

       This function is useful for the opposite case: you are going
       over the meshs through the channel-based IDs, but you need to
       know what ID to use for garray.
*/
size_t
gal_mesh_ch_based_id_from_gid(struct gal_mesh_params *mp, size_t gid)
{
  if(mp->nch==1 || mp->garray1==mp->cgarray1)
    return gid;
  else
    {
      /* The X and Y of this mesh in the full image: */
      size_t f0=gid/(mp->gs1*mp->nch1), f1=gid%(mp->gs1*mp->nch1);

      /* The channel ID: */
      size_t chid = (f0/mp->gs0) * mp->nch1 + f1/mp->gs1;

      /* The ID of this mesh in this channel: */
      size_t inchannelid = (f0%mp->gs0) * mp->gs1 + f1%mp->gs1;

      /* For a check:
      printf("%lu:\n\t(f0, f1): (%lu, %lu)\t chid: %lu\tinchannelid: %lu\n",
             gid, f0, f1, chid, inchannelid);
      */

      /* Return the channel-based ID: */
      return chid * mp->nmeshc + inchannelid;
    }

  /* This function should not reach here! If it does there is problem,
     so just return an impossible value so the root of the issue can
     be found easily: */
  return (size_t)(-1);
}





/* Get the garray-based ID from the channel-based ID. See the comments above
   gal_mesh_ch_based_id_from_gid for a complete explanation. */
size_t
gal_mesh_gid_from_ch_based_id(struct gal_mesh_params *mp, size_t chbasedid)
{
  if(mp->nch==1 || mp->garray1==mp->cgarray1)
    return chbasedid;
  else
    {
      /* The X and Y positions of this channel in the channels array: */
      size_t chx=(chbasedid/mp->nmeshc)/mp->nch1;
      size_t chy=(chbasedid/mp->nmeshc)%mp->nch1;

      /* The X and Y of this mesh in this channel: */
      size_t mx=(chbasedid%mp->nmeshc)/mp->gs1;
      size_t my=(chbasedid%mp->nmeshc)%mp->gs1;

      /* For a check:
      printf("%lu:\n\t(chx, chy): (%lu, %lu)"
             "\n\t(mx, my): (%lu, %lu)"
             "\n\t%lu\n\n",
             chbasedid, chx, chy, mx, my,
             (chx*mp->gs0+mx) * mp->nch1 + (chy*mp->gs1+my));
      */

      /* Return the : */
      return (chx*mp->gs0+mx) * mp->nch1 + (chy*mp->gs1+my);
    }

  /* This function should not reach here! If it does there is problem,
     so just return an impossible value so the root of the issue can
     be found easily: */
  return (size_t)(-1);
}





/* The user has a pixel index in the final image and wants to know the
   id it should plug into the garrays to get a value for this
   pixel. This is the job of this function. So to find the value on
   the mesh grid for a pixel at index `ind', the user should just run:

       mp->garray[imgindextomeshid(mp, ind)]
 */
size_t
gal_mesh_img_xy_to_mesh_id(struct gal_mesh_params *mp, size_t x, size_t y)
{
  /* Take the proper action. The ternary conditional is here because
     when the meshsize is not an exact multiple of the the channel
     (image) size, there might be some extra pixels in the last mesh
     in each dimension which will cause trouble in the end. So without
     these checks, a pixel lying in those extra regions will be
     thought of as belongin to another mesh (that doesn't exist). We
     have to make sure that doesn't happen. */
  if(mp->nch==1)
    return ( (x/mp->meshsize<mp->gs0 ? x/mp->meshsize : x/mp->meshsize -1)
             * mp->gs1
             + (y/mp->meshsize<mp->gs1 ? y/mp->meshsize : y/mp->meshsize -1));
  else
    {
      /* Number of pixels along each axis in all channels: */
      size_t cps0=mp->s0/mp->nch2, cps1=mp->s1/mp->nch1;

      /* The X and Y positions of this channel in the channels array: */
      size_t chx=x/cps0, chy=y/cps1;

      /* The X and Y of this mesh in this channel: */
      size_t mx = (x%cps0)/mp->meshsize;
      size_t my = (y%cps1)/mp->meshsize;

      /* If the last mesh doesn't have the same size as the rest, mx
         or my might become one larger. Note that we have already made
         sure that this pixel is in the channel specified by chx and
         chy. */
      mx = mx<mp->gs0 ? mx : mx-1;
      my = my<mp->gs1 ? my : my-1;


      /* Return the proper id to input into garray. */
      if(mp->garray1==mp->cgarray1)
        return mp->nmeshc * (chx*mp->nch1+chy) + mx * mp->gs1 + my;
      else
        return (chx*mp->gs0+mx) * mp->gs1 + (chy * mp->gs1 + my);
    }

  /* This function should not reach here! So we will just return a
     value that will always be problematic ;-). */
  return (size_t) -1;
}




















/*********************************************************************/
/********************         Full garray         ********************/
/*********************************************************************/
/* By default, the garray1 and garray2 arrays keep the meshes of each
   channel contiguously. So in practice, each channel is like a small
   independent image. This will cause problems when we want to work on
   the meshs irrespective of which channel they belong to. This
   function allocates and fills in the fgarray1 and fgarray2 arrays.

   The explanation above is for the case when reverse==0. If it is set
   equal to 1 (or any non-zero number), then
*/
void
gal_mesh_full_garray(struct gal_mesh_params *mp, int reverse)
{
  size_t nch1=mp->nch1;
  size_t ind, gs1=mp->gs1, gs0=mp->gs0;
  float *fgarray1=NULL, *fgarray2=NULL;
  float *cgarray1=mp->cgarray1, *cgarray2=mp->cgarray2;
  size_t g0, g1, f0, f1, fmeshind, chid, is1=mp->nch1*mp->gs1;

  /* First allocate the fgarrays if they were not already
     allocated. */
  if(mp->fgarray1==NULL)
    {
      /* A simple sanity check */
      if(reverse)
        error(EXIT_FAILURE, 0, "a bug!  Please contact us at %s so we can "
              "fix this problem.  For some reason, gal_mesh_full_garray "
              "has been called with the `reverse' flag set to true while "
              "fgarray is not allocated! This should not happen",
              PACKAGE_BUGREPORT);

      /* Allocate the fgarrays */
      errno=0; mp->fgarray1=malloc(mp->nmeshi*sizeof *mp->fgarray1);
      if(mp->fgarray1==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for mp->fgarray1 (mesh.c)",
              mp->nmeshi*sizeof *mp->fgarray1);
    }
  if(mp->ngarrays==2 && mp->fgarray2==NULL)
    {
      errno=0;
      mp->fgarray2=malloc(mp->nmeshi*sizeof *mp->fgarray2);
      if(mp->fgarray2==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for mp->fgarray2 (mesh.c)",
              mp->nmeshi*sizeof *mp->fgarray2);
    }
  fgarray1=mp->fgarray1;
  fgarray2=mp->fgarray2;

  /* Fill the fgarrays or cgarrays with the proper values: */
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
          if(mp->ngarrays==2)
            {
              if(reverse)
                cgarray2[ind+fmeshind]=fgarray2[(g0+f0)*is1+g1+f1];
              else
                fgarray2[(g0+f0)*is1+g1+f1]=cgarray2[ind+fmeshind];
            }
        }
    }

  /* Just for a check:
  gal_fits_array_to_file("nochannels.fits", "fgarray1", FLOAT_IMG,
                         fgarray1, mp->nch2*mp->gs0, mp->nch1*mp->gs1, 1,
                         NULL, NULL, "mesh");
  if(mp->ngarrays==2)
    gal_fits_array_to_file("nochannels.fits", "fgarray2", FLOAT_IMG,
                           fgarray2, mp->nch2*mp->gs0, mp->nch1*mp->gs1,
                           1, NULL, NULL, "mesh");
  */
}



















/*********************************************************************/
/********************     Checking functions      ********************/
/*********************************************************************/
/* Save the meshid of each pixel into an array the size of the image. */
void
gal_mesh_check_mesh_id(struct gal_mesh_params *mp, long **out)
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
    error(EXIT_FAILURE, errno,
          "the array to show mesh labels in checkmesh");

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
gal_mesh_check_garray(struct gal_mesh_params *mp, float **out1,
                      float **out2)
{
  int ngarrays=mp->ngarrays;
  size_t gid, row, start, chbasedid, *types=mp->types;
  size_t s0, s1, is1=mp->s1, *ts0=mp->ts0, *ts1=mp->ts1;
  float *f, *fp, *ff=NULL, garray1=FLT_MAX, garray2=FLT_MAX;

  /* Allocate the array to keep the mesh indexs. Calloc is used so we
     can add all the indexes to the existing value to make sure that
     there is no overlap. */
  errno=0; *out1=malloc(mp->s0*mp->s1*sizeof **out1);
  if(*out1==NULL)
    error(EXIT_FAILURE, errno,
          "%lu bytes for out1 in gal_mesh_check_garray (mesh.c)",
          mp->s0*mp->s1*sizeof **out1);
  if(ngarrays==2)
    {
      errno=0; *out2=malloc(mp->s0*mp->s1*sizeof **out2);
      if(*out2==NULL)
        error(EXIT_FAILURE, errno,
              "%lu bytes for out2 in gal_mesh_check_garray (mesh.c)",
              mp->s0*mp->s1*sizeof **out2);
    }

  /* Fill the array: */
  for(gid=0;gid<mp->nmeshi;++gid)
    {
      /* Set the proper meshid depending on what garray points to, see
         the explanation above setmeshid. */
      chbasedid = gal_mesh_ch_based_id_from_gid(mp, gid);

      /* Fill the output array with the value in this mesh. It is
         really important that `i' should be used for the garrays, not
         cmeshid. cmeshid is only for the basic mesh parameters that
         it is used for. */
      row=0;
      s0=ts0[types[chbasedid]];
      s1=ts1[types[chbasedid]];
      start=mp->start[chbasedid];
      garray1 = mp->garray1[gid];
      if(ngarrays==2)
        garray2 = mp->garray2[gid];
      do
        {
          fp= ( f = *out1 + start + row * is1 ) + s1;
          if(ngarrays==2) ff= *out2 + start + row * is1;
          do
            {
              *f++ = garray1;
              if(ngarrays==2) *ff++ = garray2;
            }
          while(f<fp);
          ++row;
        }
      while(row<s0);
    }
}





/* Save the mesh grid values into an output file. */
void
gal_mesh_value_file(struct gal_mesh_params *mp, char *filename,
                    char *extname1, char *extname2, struct wcsprm *wcs,
                    char *spack_string)
{
  float *tmp1=NULL, *tmp2=NULL;

  if(mp->meshbasedcheck)
    {
      /* We want one pixel per mesh. If the last operation was on
         cgarray, then first the fgarray has to be filled. Note that
         when more than one channel is present, only fgarray can be
         used for this job. In cgarray the meshs are ordered
         differently. */
      if(mp->garray1==mp->cgarray1) gal_mesh_full_garray(mp, 0);
      gal_fits_array_to_file(filename, extname1, FLOAT_IMG,
                             mp->fgarray1, mp->gs0*mp->nch2,
                             mp->gs1*mp->nch1, 0, wcs, NULL,
                             spack_string);
      if(mp->ngarrays==2)
        /* Note that gal_mesh_full_garray will correct both the meshs if
           there are two.*/
        gal_fits_array_to_file(filename, extname2, FLOAT_IMG,
                               mp->fgarray2, mp->gs0*mp->nch2,
                               mp->gs1*mp->nch1, 0, wcs, NULL,
                               spack_string);

    }
  else
    {
      gal_mesh_check_garray(mp, &tmp1, &tmp2);
      gal_fits_array_to_file(filename, extname1, FLOAT_IMG, tmp1,
                             mp->s0, mp->s1, 0, wcs, NULL,
                             spack_string);
      if(mp->ngarrays==2)
        gal_fits_array_to_file(filename, extname2, FLOAT_IMG, tmp2,
                               mp->s0, mp->s1, 0, wcs, NULL,
                               spack_string);
      free(tmp1);
      free(tmp2);
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
static void
fillmeshinfo(struct gal_mesh_params *mp, size_t chs0, size_t chs1,
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

  /* Main meshs (type 0) in each channel. If we have more than one row
     or column of meshes in each channel, then we will have some type1
     meshes.*/
  if(gs0>1 || gs1>1)
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

  /* Top row of meshes (type 1) in each channel. Note that when there
     is only one column of meshes, then the only one mesh remaining
     will be type 4, not type 1 (here).*/
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

  /* Left column of meshes (type 2) in each channel. When there is
     only one row of meshs, then the last single remaining mesh will
     be type 3, not type 2 (here).*/
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
    error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we can fix "
          "it. The basic information for some meshes has not been found "
          "(in fillmeshinfo of mesh.c)", PACKAGE_BUGREPORT);
}





/* In the explanation below, all the parameters named are from the
   `struct gal_mesh_params` of `mesh.h`.

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
gal_mesh_make_mesh(struct gal_mesh_params *mp)
{
  size_t meshsize=mp->meshsize, lasts0, lasts1;
  size_t i, chs0=mp->s0/mp->nch2, chs1=mp->s1/mp->nch1;

  /* Set the number of channels. */
  mp->nch=mp->nch1*mp->nch2;

  /* Incase meshsize is larger than the channel sizes, make it equal
     to the smaller of the channel axis sizes. */
  if(meshsize>chs0 || meshsize>chs1)
    mp->meshsize = meshsize = (chs0<chs1
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
  if(mp->start==NULL) error(EXIT_FAILURE, errno, "mesh starting points");
  errno=0; mp->types=malloc(mp->nmeshi*sizeof *mp->types);
  if(mp->types==NULL) error(EXIT_FAILURE, errno, "mesh types");
  errno=0; mp->chindex=malloc(mp->nmeshi*sizeof *mp->chindex);
  if(mp->chindex==NULL) error(EXIT_FAILURE, errno, "mesh in channel index");
  errno=0; mp->imgindex=malloc(mp->nmeshi*sizeof *mp->imgindex);
  if(mp->imgindex==NULL) error(EXIT_FAILURE, errno, "mesh in image index");

  /* Distribute the meshes in all the threads. */
  gal_threads_dist_in_threads(mp->nmeshi, mp->numthreads, &mp->indexs,
                              &mp->thrdcols);

  /* Fill in the information for each mesh and each type. */
  fillmeshinfo(mp, chs0, chs1, lasts0, lasts1);
}





void
gal_mesh_free_mesh(struct gal_mesh_params *mp)
{
  free(mp->start);
  free(mp->types);
  free(mp->indexs);
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

   1. A pointer to the gal_mesh_params structure that keeps all the
      information.

   2. A pointer to a function that returns and gets a `void *' as its only
      argument. This function will be directly given to
      pthread_create. Through this argument, you can choose the function to
      operate on the mesh grid.

   3. The size of each element to copy the mesh grid into, this has to
      be type size of the same type that constitutes `img' in
      gal_mesh_params. If the value to this argument is non-zero, an array
      will be allocated that can contain all the pixels in all the
      meshs and can be used by threads to manipute the pixels (for
      example sort them) in each mesh.

   4. If the value of this argument is 1, then a second garray will be
      allocated in case your operation needs it.

   5. Wether the allocated garrays should be initialized or not.
*/
void
gal_mesh_operate_on_mesh(struct gal_mesh_params *mp,
                         void *(*meshfunc)(void *), size_t oneforallsize,
                         int makegarray2, int initialize)
{
  int err;
  size_t i, nb;
  float *f, *fp;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  pthread_attr_t attr;
  struct gal_mesh_thread_params *mtp;
  size_t numthreads=mp->numthreads;

  /* Allocate the arrays to keep the thread and parameters for each
     thread. */
  errno=0; mtp=malloc(numthreads*sizeof *mtp);
  if(mtp==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes in fillmesh (mesh.c) for mtp",
          numthreads*sizeof *mtp);

  /* Set the number of garrays to operate on: */
  mp->ngarrays = makegarray2 ? 2 : 1;

  /* If cgarrays have not been allocated before, allocate them. */
  if(mp->cgarray1==NULL)
    {
      errno=0; mp->cgarray1=malloc(mp->nmeshi*sizeof *mp->cgarray1);
      if(mp->cgarray1==NULL) error(EXIT_FAILURE, errno, "mp->cgarray1");
    }
  if(mp->ngarrays==2 && mp->cgarray2==NULL)
    {
      errno=0; mp->cgarray2=malloc(mp->nmeshi*sizeof *mp->cgarray2);
      if(mp->cgarray2==NULL) error(EXIT_FAILURE, errno, "mp->cgarray2");
    }

  /* Initialize the cgarrays to NaN:*/
  if(initialize)
    {
      fp=(f=mp->cgarray1)+mp->nmeshi; do *f++=NAN; while(f<fp);
      if(mp->ngarrays==2)
        { fp=(f=mp->cgarray2)+mp->nmeshi; do *f++=NAN; while(f<fp); }
      mp->garray1=mp->cgarray1;
      mp->garray2=mp->cgarray2;
    }

  /* `oneforall' is an array with the sides equal to the maximum side
     of the meshes in the image. The purpose is to enable manipulating
     the mesh pixels (for example sorting them and so on.). One such
     array is allocated for each thread in this one full
     allocation. Each thread can then use its own portion of this
     array through the following declaration:

     float *oneforall=&mp->oneforall[mtp->id*mp->maxs0*mp->maxs1];

     In gal_mesh_params, `oneforall' is defined as a `void *', so the
     caller function, can cast it to any type it wants. The size of
     each type is given to `fillmesh' through the `oneforallsize'
     argument.
  */
  if(oneforallsize)
    {
      errno=0;
      mp->oneforall=malloc(numthreads*mp->maxs0*mp->maxs1*oneforallsize);
      if(mp->oneforall==NULL)
        error(EXIT_FAILURE, errno, "unable to allocate %lu bytes for"
              "mtp->oneforall in fillmesh of mesh.c",
              numthreads*mp->maxs0*mp->maxs1*oneforallsize);
    }

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
      gal_threads_attr_barrier_init(&attr, &mp->b, nb);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
        if(mp->indexs[i*mp->thrdcols]!=GAL_THREADS_NON_THRD_INDEX)
          {
            mtp[i].id=i;
            mtp[i].mp=mp;
            err=pthread_create(&t, &attr, meshfunc, &mtp[i]);
            if(err) error(EXIT_FAILURE, 0, "can't create thread %lu", i);
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&mp->b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&mp->b);
    }

  free(mtp);
  if(oneforallsize) free(mp->oneforall);
}




















/*********************************************************************/
/********************         Interpolate         ********************/
/*********************************************************************/
static void
preparemeshinterparrays(struct gal_mesh_params *mp)
{
  size_t bs0=mp->gs0, bs1=mp->gs1;
  size_t numthreads=mp->numthreads;


  /* If full-interpolation is to be done (ignoring channels) we will
     need two arrays to temporarily keep the actual positions of the
     meshs. Note that by default, the meshs in each channels are
     placed contiguously. In these arrays, the meshs are placed as
     they would in the image (channels are ignored).*/
  if(mp->fullinterpolation)
    {
      /* In case the previous operation was on cgarrays, then you have
         to fill in fgarray. */
      if(mp->garray1==mp->cgarray1)
        gal_mesh_full_garray(mp, 0);
      bs0=mp->nch2*mp->gs0;
      bs1=mp->nch1*mp->gs1;
      mp->garray1=mp->fgarray1;
      if(mp->ngarrays==2) mp->garray2=mp->fgarray2;
    }
  else
    {
      mp->garray1=mp->cgarray1;
      if(mp->ngarrays==2) mp->garray2=mp->cgarray2;
    }


  /* Allocate the output arrays to keep the final values. Note that we
     cannot just do the interpolaton on the same input grid, because
     the newly filled interpolated pixels will affect the later
     ones. In the end, these copies are going to replace the
     garrays. */
  errno=0; mp->outgarray1=malloc(mp->nmeshi*sizeof *mp->outgarray1);
  if(mp->outgarray1==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for outgarray1 (mesh.c)",
          mp->nmeshi*sizeof *mp->outgarray1);
  if(mp->ngarrays==2)
    {
      errno=0; mp->outgarray2=malloc(mp->nmeshi*sizeof *mp->outgarray2);
      if(mp->outgarray2==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for outgarray2 (mesh.c)",
              mp->nmeshi*sizeof *mp->outgarray2);
    }


  /* Allocate the array for the values of the nearest pixels. For each
     pixel in each thread, the nearest pixels will be put here and
     sorted to find the median value. If garray2 is not NULL, then it
     should be interpolated. Note that both arrays have blank pixels
     on the same places.*/
  errno=0;
  mp->nearest1=malloc(numthreads*mp->numnearest*sizeof *mp->nearest1);
  if(mp->nearest1==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for the array to keep the "
          "nearest1 values for interpolation (mesh.c)",
          numthreads*mp->numnearest*sizeof *mp->nearest1);
  if(mp->ngarrays==2)
    {
      errno=0;
      mp->nearest2=malloc(numthreads*mp->numnearest*sizeof *mp->nearest2);
      if(mp->nearest2==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for the array to keep the "
              "nearest2 values for interpolation (mesh.c)",
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





/* Ideally, this should be a radial distance using the square root of
   the differences. But here we are not dealing with subpixel
   distances, so manhattan distance is fine. It also avoids the square
   root function which is much slower. */
static float
manhattandistance(long ind, long xc, long yc, long s1)
{
  return labs(ind%s1 - xc) + labs(ind/s1 - yc);
}




/* Some of the variables have different names than the gal_mesh_params
   structure because they are to be fed into the
   GAL_NEIGHBORS_FILL_4_ALLIMG macro. */
static void *
meshinterponthread(void *inparams)
{
  struct gal_mesh_thread_params *mtp=(struct gal_mesh_thread_params *)inparams;
  struct gal_mesh_params *mp=mtp->mp;

  /* Basic variables used in other definitions: */
  size_t numnearest=mp->numnearest;
  size_t is0=mp->fullinterpolation ? mp->gs0*mp->nch2 : mp->gs0;
  size_t is1=mp->fullinterpolation ? mp->gs1*mp->nch1 : mp->gs1;

  /* Variables for this function: */
  struct gal_linkedlist_tosll *lQ, *sQ;
  int ngarrays=mp->ngarrays;
  size_t xc, yc, *n, *nf, currentnum, thisind;
  unsigned char *byt=&mp->byt[mtp->id*is0*is1];
  float *nearest1=&mp->nearest1[mtp->id*numnearest];
  size_t i, *indexs=&mp->indexs[mtp->id*mp->thrdcols];
  float mdist, *garray1=mp->garray1, *garray2=mp->garray2;
  size_t fmeshid, position, *ind=&position, numngb, ngb[4];
  float *outgarray1=mp->outgarray1, *outgarray2=mp->outgarray2;
  float *nearest2 = ngarrays==2 ? &mp->nearest2[mtp->id*numnearest] : NULL;


  /* Go over all the meshes for this thread. */
  for(i=0;indexs[i]!=GAL_THREADS_NON_THRD_INDEX;++i)
    {
      /* Get the index of this NaN mesh: */
      thisind=indexs[i];

      /* If this mesh is not blank and the user has only asked to
         interpolate blank pixels, then set the final values and go
         onto the next mesh. */
      if(mp->interponlyblank && !isnan(garray1[thisind]))
        {
          outgarray1[thisind]=garray1[thisind];
          if(ngarrays==2) outgarray2[thisind]=garray2[thisind];
          continue;
        }

      /* Reset the byt array: */
      memset(byt, 0, is0*is1);

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
      gal_linkedlist_add_to_tosll_end( &lQ, &sQ, *ind, 0 );

      /* Start finding the nearest filled pixels. */
      while(sQ)
        {
          /* Pop out a pixel index (p) from the queue: */
          gal_linkedlist_pop_from_tosll_start(&lQ, &sQ, ind, &mdist);

          /* If it isn't a NaN, then put it in the `nearest1' and
             `nearest2' arrays. */
          if(!isnan(garray1[*ind+fmeshid]))
            {
              nearest1[currentnum]=garray1[*ind+fmeshid];
              if(ngarrays==2) nearest2[currentnum]=garray2[*ind+fmeshid];
              if(++currentnum>=numnearest) break;
            }

          /* Check the four neighbors and if they have not already
             been checked, put them into the queue. */
          GAL_NEIGHBORS_FILL_4_ALLIMG;
          nf=(n=ngb)+numngb;
          do
            if(byt[*n]==0)
              {
                byt[*n]=1;
                gal_linkedlist_add_to_tosll_end(&lQ, &sQ, *n,
                                                manhattandistance(*n, xc,
                                                                  yc, is1));
              }
          while(++n<nf);

          /* If there are no more meshes to add to the queue, then
             this shows, there were not enough points for
             interpolation. Normally, this loop should only be exited
             through the `currentnum>=numnearest' check above. */
          if(sQ==NULL)
            error(EXIT_FAILURE, 0, "%s: only %lu mesh(s) are filled for "
                  "interpolation in channel %lu. Either set less "
                  "restrictive requirements to get more interpolation "
                  "points or decrease the number of nearest points to "
                  "use for interpolation. Problem encountered on thread "
                  "%lu, for pixel %lu. When running on multiple threads, "
                  "This message might be repeated for different threads",
                  mp->errstart, currentnum, thisind/mp->nmeshc, mtp->id,
                  thisind);
        }
      gal_linkedlist_tosll_free(lQ);  /* Rest of the queue not needed. */


      /* Find the median of the nearest neighbors and put it in: */
      qsort(nearest1, numnearest, sizeof *nearest1,
            gal_qsort_float_increasing);
      outgarray1[thisind] = ( numnearest%2 ?
                              nearest1[numnearest/2] : /* Odd.  */
                              (nearest1[numnearest/2]  /* Even. */
                               +nearest1[numnearest/2-1])/2 );
      if(ngarrays==2)
        {
          qsort(nearest2, numnearest, sizeof *nearest2,
                gal_qsort_float_increasing);
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
gal_mesh_interpolate(struct gal_mesh_params *mp, char *errstart)
{
  int err;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  size_t i, nb;
  pthread_attr_t attr;
  struct gal_mesh_thread_params *mtp;
  size_t numthreads=mp->numthreads;

  /* Prepare all the gal_mesh_params arrays: */
  mp->errstart=errstart;
  preparemeshinterparrays(mp);

  /* Allocate the arrays to keep the thread and parameters for each
     thread. */
  errno=0; mtp=malloc(numthreads*sizeof *mtp);
  if(mtp==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes in fillmesh (mesh.c) for mtp",
          numthreads*sizeof *mtp);

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
      if(mp->nmeshi<numthreads) nb=mp->nmeshi+1;
      else                      nb=numthreads+1;
      gal_threads_attr_barrier_init(&attr, &mp->b, nb);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
        if(mp->indexs[i*mp->thrdcols]!=GAL_THREADS_NON_THRD_INDEX)
          {
            mtp[i].id=i;
            mtp[i].mp=mp;
            err=pthread_create(&t, &attr, meshinterponthread, &mtp[i]);
            if(err) error(EXIT_FAILURE, 0, "can't create thread %lu", i);
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&mp->b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&mp->b);
    }

  /* Replace garray1 and garray2 with outgarray1 and outgarray2 */
  if(mp->fullinterpolation)
    { free(mp->fgarray1); mp->garray1=mp->fgarray1=mp->outgarray1; }
  else
    { free(mp->cgarray1); mp->garray1=mp->cgarray1=mp->outgarray1; }
  if(mp->ngarrays==2)
    {
      if(mp->fullinterpolation)
        { free(mp->fgarray2); mp->garray2=mp->fgarray2=mp->outgarray2; }
      else
        { free(mp->cgarray2); mp->garray2=mp->cgarray2=mp->outgarray2; }
    }
  mp->outgarray1=mp->outgarray2=NULL;

  /* For a check
  system("rm test.fits");
  gal_fits_array_to_file("test.fits", "garray1", FLOAT_IMG, mp->garray1,
                          mp->nch2*mp->gs0, mp->nch1*mp->gs1, 1, NULL, NULL,
                          "mesh");
  gal_fits_array_to_file("test.fits", "garray2", FLOAT_IMG, mp->garray2,
                         mp->nch2*mp->gs0, mp->nch1*mp->gs1, 1, NULL, NULL,
                         "mesh");
  */

  /* Clean up. */
  free(mtp);
  free(mp->byt);
  free(mp->nearest1);
  if(mp->ngarrays==2) free(mp->nearest2);
}




















/*********************************************************************/
/********************           Smooth            ********************/
/*********************************************************************/
void
gal_mesh_smooth(struct gal_mesh_params *mp)
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
      /* In case the previous operation was on the cgarrays. garray2
         should not be checked. garray1 and garray2 are always in sync
         during one operation. However, if a previous operation had
         set garray2 and this operation only uses garray1, then things
         are going to go wrong. So Based on the principle that garray1
         and garray2 are always going to be in sync with each other,
         we don't the mp->garray2==mp->cgarray2 check should not be
         done. */
      if(mp->garray1==mp->cgarray1)
        gal_mesh_full_garray(mp, 0);

      /* Do the spatial convolution */
      gal_spatialconvolve_convolve(mp->fgarray1, gs0*mp->nch2, gs1*mp->nch1,
                                   kernel, smoothwidth, smoothwidth,
                                   numthreads, 1, &sgarray1);

      free(mp->fgarray1);
      mp->garray1=mp->fgarray1=sgarray1;
      if(mp->ngarrays==2)
        {
          gal_spatialconvolve_convolve(mp->fgarray2, gs0*mp->nch2,
                                       gs1*mp->nch1, kernel, smoothwidth,
                                       smoothwidth, mp->numthreads,
                                       1, &sgarray2);
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
          gal_mesh_full_garray(mp, 1);
        mp->garray1=mp->cgarray1;
        mp->garray2=mp->cgarray2;

        charray=&mp->cgarray1[chid*nmeshc];
        gal_spatialconvolve_convolve(charray, gs0, gs1, kernel, smoothwidth,
                                     smoothwidth, numthreads, 1, &tmp);
        o=tmp; fp=(f=charray)+gs0*gs1; do *f=*o++; while(++f<fp);
        free(tmp);
        if(mp->ngarrays==2)
          {
            charray=&mp->cgarray2[chid*nmeshc];
            gal_spatialconvolve_convolve(charray, gs0, gs1, kernel,
                                         smoothwidth, smoothwidth,
                                         numthreads, 1, &tmp);
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
static void*
meshspatialconvonthreads(void *inparam)
{
  struct gal_mesh_thread_params *mtp = (struct gal_mesh_thread_params *)inparam;
  struct gal_mesh_params *mp = mtp->mp;
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
  for(i=0; (meshid=indexs[i])!=GAL_THREADS_NON_THRD_INDEX; ++i)
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
gal_mesh_spatial_convolve_on_mesh(struct gal_mesh_params *mp, float **conv)
{
  int err;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  pthread_attr_t attr;
  size_t i, nb, *chbrd;
  struct gal_mesh_thread_params *mtp;
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
    error(EXIT_FAILURE, errno, "%lu bytes for convolution on mesh output",
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
          "gal_mesh_spatial_convolve_on_mesh", mp->nch*4*sizeof *chbrd);
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
      gal_threads_attr_barrier_init(&attr, &mp->b, nb);

      /* Spin off the threads: */
      for(i=0;i<numthreads;++i)
        if(mp->indexs[i*mp->thrdcols]!=GAL_THREADS_NON_THRD_INDEX)
          {
            mtp[i].id=i;
            mtp[i].mp=mp;
            mtp[i].conv=*conv;
            mtp[i].chbrd=chbrd;
            err=pthread_create(&t, &attr, meshspatialconvonthreads, &mtp[i]);
            if(err) error(EXIT_FAILURE, 0, "can't create thread %lu", i);
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&mp->b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&mp->b);
    }

  free(mtp);
  free(chbrd);
}






/* The indexs array for correcting the convolution on inner channel edges
   has been allocated.  Note that gal_threads_dist_in_threads will
   distribute indexs from zero to numpix-1.  After it, we should fill in
   all the channels.

   The method of filling in the indexs array with the proper indexs to
   re-convolve is very similar to the method explained below in
   gal_mesh_change_to_full_convolution, where it is explained how to count
   the number of pixels that should be re-convolved. */
static void
corrconvindexs(struct gal_mesh_params *mp, size_t **indexs,
               size_t *numpix, size_t *thrdcols)
{
  size_t i, j, a, b;
  size_t numthreads=mp->numthreads;
  size_t alow, ahigh, blow, bhigh, *ind;
  size_t npch0=mp->s0/mp->nch2, npch1=mp->s1/mp->nch1;
  size_t nch1=mp->nch1, nch2=mp->nch2, ks0=mp->ks0, ks1=mp->ks1;
  size_t s0=mp->s0, s1=mp->s1, hk0=mp->ks0/2+1, hk1=mp->ks1/2+1;

  /* Find the number of pixels that must be convolved. Note that when
     we don't care about the edges, on each dimension, we have one
     less border than the number of channels in that dimention.

     On each channel's border, we want to re-convolve the ks0/2+1
     pixels before the channel edge along the first dimension. So in
     total, for each channel edge, we want ks0+1 pixels on its two
     sides. Don't forget that the kernel sides are odd.

     So on the first dimension, we have (ks0+1)*(nch2-1) pixels that
     should be re-convolved. Multiplying it by s1, we get the total
     number of pixels in 2D around the first axis internal edges.

     Along the second dimension, we have (ks1+1)*(nch1-1) pixels. But
     this time we can't just multiply by s0 to get the total number of
     2D pixels. Because of the overlap. So we only have to multily by
     the number of rows that were not accounted in the first run,
     which is s0-(ks0+1)*(nch2-1). */
  *numpix = ( (ks0+1)*(nch2-1)*s1
             + (ks1+1)*(nch1-1)*(s0-(ks0+1)*(nch2-1)) );


  /* Distribute the indexs of the desired pixels into the indexs
     array. */
  gal_threads_dist_in_threads(*numpix, numthreads, indexs, thrdcols);

  ind=*indexs;
  for(i=1;i<nch2;++i)           /* FIRST LOOP. */
    {
      alow  = i*npch0<hk0    ? 0  : i*npch0-hk0;
      ahigh = i*npch0+hk0>s0 ? s0 : i*npch0+hk0;

      /* For a check
      printf("1: alow: %-4lu ahigh: %-4lu.\t(0 -- %lu)\n",
             alow, ahigh, s1);
      */
      for(a=alow; a<ahigh; ++a)
        for(b=0;b<s1;++b)
          {
            while(*ind==GAL_THREADS_NON_THRD_INDEX) ++ind;
            *ind++=a*s1+b;
          }
    }
  for(j=1;j<nch1;++j)           /* SECOND LOOP */
    {
      blow  = j*npch1<hk1    ? 0  : j*npch1-hk1;
      bhigh = j*npch1+hk1>s1 ? s1 : j*npch1+hk1;
      for(b=blow; b<bhigh; ++b)
        {
          /* Since there might be multiple channels along the first C axis
             and we want the spaces between their borders, alow is
             initiated with 0.  */
          alow=0;

          /* Check the areas under the bordering area for each
             vertical channel row: */
          for(i=1;i<nch2;++i)
            {
              /* This `ahigh' is actually the alow in the FIRST LOOP. */
              ahigh = i*npch0<hk0    ? 0  : i*npch0-hk0;

              /* For a check:
              printf("2: alow: %-4lu ahigh: %-4lu. b: %-10lu "
                     "(%lu -- %lu)\n", alow, ahigh, b, blow, bhigh);
              */
              for(a=alow;a<ahigh;++a)
                {
                  while(*ind==GAL_THREADS_NON_THRD_INDEX) ++ind;
                  *ind++=a*s1+b;
                }
              /* This alow is actually the ahigh in the FIRST LOOP. */
              alow = i*npch0+hk0>s0 ? s0 : i*npch0+hk0;
            }

          /* All the areas under were checked, now we have to add the
             area ontop of the highest middle channel edge. Note that
             when there is no channel along the first C axis
             (nch2==1), then it will immediately come here and scan
             all the rows between the blow and bhigh columns. */
          ahigh=s0;

          /* For a check:
          printf("3: alow: %-4lu ahigh: %-4lu. b: %-10lu (%lu -- %lu)\n",
                 alow, ahigh, b, blow, bhigh);
          */
          for(a=alow;a<ahigh;++a)
            {
              while(*ind==GAL_THREADS_NON_THRD_INDEX) ++ind;
              *ind++=a*s1+b;
            }
        }
    }
}





/* Convolution is a very expensive (time consuming) operation. The
   problem is that sometimes, convolution was done on each channel
   independently, but later on, a program might need convolution over
   the full image (for example, during the program the differences
   between the channels has been calculated and removed), so it is
   very important to remove the discontinuities that the initial
   convolution on each channel can cause.

   This is the job of this function. It recieves a mesh structure and
   an already convolved image. Only the pixels lying within half of
   the PSF width of channel borders are chosen and only they are
   convolved (this time over the full image). Most of the image pixels
   (whose distance from the channel edges is more than half the PSF),
   do not need to undergo convolution again.

   Note that the pixels on the edges of the image do not need to undergo
   this correction.  Basically this function is very similar to
   gal_spatialconvolve_convolve (spatialconvolve.c), other than the fact
   that the indexs are not over the full image but only a select number of
   pixels.
*/
void
gal_mesh_change_to_full_convolution(struct gal_mesh_params *mp, float *conv)
{
  int err;
  pthread_t t;          /* All thread ids saved in this, not used. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct gal_spatialconvolve_params *scp;
  size_t i, nb, *indexs, numpix, thrdcols;

  /* If convolution was done over the full image, then there is
     nothing this function should do so just return. After this
     function, the image is fully convolved, so fullconvolution should
     be set to 1. */
  if(mp->nch==1 || mp->fullconvolution) return;
  mp->fullconvolution=1;


  /* Array keeping thread parameters for each thread.*/
  errno=0;
  scp=malloc(mp->numthreads*sizeof *scp);
  if(scp==NULL)
    error(EXIT_FAILURE, errno,
          "%lu bytes for scp in gal_mesh_change_to_full_convolution "
          "(mesh.c)", mp->numthreads*sizeof *scp);


  /* Put the indexs of the pixels to re-convolve here. */
  corrconvindexs(mp, &indexs, &numpix, &thrdcols);


  /* Start the convolution on the desired pixels. */
  if(mp->numthreads==1)
    {
      gal_spatialconvolve_pparams(mp->img, mp->s0, mp->s1, mp->kernel,
                                  mp->ks0, mp->ks1, mp->numthreads, 1,
                                  conv, indexs, &scp[0]);
      gal_spatialconvolve_thread(&scp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
         (that spinns off the nt threads) is also a thread, so the
         number the barrier should be one more than the number of
         threads spinned off. */
      if(numpix<mp->numthreads) nb=numpix+1;
      else                      nb=mp->numthreads+1;
      gal_threads_attr_barrier_init(&attr, &b, nb);

      /* Spin off the threads: */
      for(i=0;i<mp->numthreads;++i)
        if(indexs[i*thrdcols]!=GAL_THREADS_NON_THRD_INDEX)
          {
            scp[i].b=&b;
            gal_spatialconvolve_pparams(mp->img, mp->s0, mp->s1, mp->kernel,
                                        mp->ks0, mp->ks1, mp->numthreads, 1,
                                        conv, &indexs[i*thrdcols], &scp[i]);
            err=pthread_create(&t, &attr, gal_spatialconvolve_thread,
                               &scp[i]);
            if(err)
              error(EXIT_FAILURE, 0, "can't create thread %lu", i);
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }

  free(scp);
  free(indexs);
}
