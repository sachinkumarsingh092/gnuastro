/*********************************************************************
match -- Functions to match catalogs and WCS.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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

#include <gsl/gsl_sort.h>

#include <gnuastro/box.h>
#include <gnuastro/list.h>
#include <gnuastro/pointer.h>
#include <gnuastro/permutation.h>









/**********************************************************************/
/*****************   Coordinate match custom list   *******************/
/**********************************************************************/
struct match_coordinate_sfll
{
  float f;
  size_t v;
  struct match_coordinate_sfll *next;
};





static void
match_coordinate_add_to_sfll(struct match_coordinate_sfll **list,
                             size_t value, float fvalue)
{
  struct match_coordinate_sfll *newnode;

  errno=0;
  newnode=malloc(sizeof *newnode);
  if(newnode==NULL)
    error(EXIT_FAILURE, errno, "%s: new node couldn't be allocated",
          __func__);

  newnode->v=value;
  newnode->f=fvalue;
  newnode->next=*list;
  *list=newnode;
}





static void
match_coordinate_pop_from_sfll(struct match_coordinate_sfll **list,
                               size_t *value, float *fvalue)
{
  struct match_coordinate_sfll *tmp;
  tmp=*list;
  *value=tmp->v;
  *fvalue=tmp->f;
  *list=tmp->next;
  free(tmp);
}




















/********************************************************************/
/*************            Coordinate matching           *************/
/*************      Sanity checks and preparations      *************/
/********************************************************************/
/* Since these checks are repetative, its easier to have a separate
   function for both inputs. */
static void
match_coordinate_sanity_check_columns(gal_data_t *coord, char *info,
                                      int inplace, int *allf64)
{
  gal_data_t *tmp;

  for(tmp=coord; tmp!=NULL; tmp=tmp->next)
    {
      if(tmp->type!=GAL_TYPE_FLOAT64)
        {
          if(inplace)
            error(EXIT_FAILURE, 0, "%s: when 'inplace' is activated, the "
                  "input coordinates must have 'float64' type. At least "
                  "one node of the %s list has type of '%s'", __func__, info,
                  gal_type_name(tmp->type, 1));
          else
            *allf64=0;
        }

      if(tmp->ndim!=1)
        error(EXIT_FAILURE, 0, "%s: each input coordinate column must have "
              "a single dimension (be a single column). Atleast one node of "
              "the %s list has %zu dimensions", __func__, info, tmp->ndim);

      if(tmp->size!=coord->size)
        error(EXIT_FAILURE, 0, "%s: the nodes of each list of coordinates "
              "must have the same number of elements. At least one node of "
              "the %s list has %zu elements while the first has %zu "
              "elements", __func__, info, tmp->size, coord->size);
    }
}





/* To keep the main function clean, we'll do the sanity checks here. */
static void
match_coordinaes_sanity_check(gal_data_t *coord1, gal_data_t *coord2,
                              double *aperture, int inplace, int *allf64)
{
  size_t ncoord1=gal_list_data_number(coord1);

  /* Make sure both lists have the same number of datasets. NOTE: they
     don't need to have the same number of elements. */
  if( ncoord1!=gal_list_data_number(coord2) )
    error(EXIT_FAILURE, 0, "%s: the two inputs have different numbers of "
          "datasets (%zu and %zu respectively)", __func__, ncoord1,
          gal_list_data_number(coord2));


  /* This function currently only works for less than 4 dimensions. */
  if(ncoord1>3)
    error(EXIT_FAILURE, 0, "%s: %zu dimension matching requested, this "
          "function currently only matches datasets with a maximum of 3 "
          "dimensions", __func__, ncoord1);

  /* Check the column properties. */
  match_coordinate_sanity_check_columns(coord1, "first", inplace, allf64);
  match_coordinate_sanity_check_columns(coord2, "second", inplace, allf64);

  /* Check the aperture values. */
  if(aperture[0]<=0)
    error(EXIT_FAILURE, 0, "%s: the first value in the aperture (%g) cannot "
          "be zero or negative", __func__, aperture[0]);
  switch(ncoord1)
    {
    case 1:  /* We don't need any checks in a 1D match. */
      break;

    case 2:
      if( (aperture[1]<=0 || aperture[1]>1))
        error(EXIT_FAILURE, 0, "%s: the second value in the aperture (%g) "
              "is the axis ratio, so it must be larger than zero and less "
              "than 1", __func__, aperture[1]);
      break;

    case 3:
      if(aperture[1]<=0 || aperture[1]>1 || aperture[2]<=0 || aperture[2]>1)
        error(EXIT_FAILURE, 0, "%s: at least one of the second or third "
              "values in the aperture (%g and %g respectively) is smaller "
              "than zero or larger than one. In a 3D match, these are the "
              "axis ratios, so they must be larger than zero and less than "
              "1", __func__, aperture[1], aperture[2]);
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
            "the issue. The value %zu not recognized for 'ndim'", __func__,
            PACKAGE_BUGREPORT, ncoord1);
    }
}





/* To keep things clean, the sorting of each input array will be done in
   this function. */
static size_t *
match_coordinates_prepare_sort(gal_data_t *coords, size_t minmapsize)
{
  gal_data_t *tmp;
  size_t *permutation=gal_pointer_allocate(GAL_TYPE_SIZE_T, coords->size, 0,
                                           __func__, "permutation");

  /* Get the permutation necessary to sort all the columns (based on the
     first column). */
  gsl_sort_index(permutation, coords->array, 1, coords->size);

  /* Sort all the coordinates. */
  for(tmp=coords; tmp!=NULL; tmp=tmp->next)
    gal_permutation_apply(tmp, permutation);

  /* Clean up. */
  return permutation;
}





/* Do the preparations for matching of coordinates. */
static void
match_coordinates_prepare(gal_data_t *coord1, gal_data_t *coord2,
                          int sorted_by_first, int inplace, int allf64,
                          gal_data_t **A_out, gal_data_t **B_out,
                          size_t **A_perm, size_t **B_perm,
                          size_t minmapsize)
{
  gal_data_t *c, *tmp, *A=NULL, *B=NULL;

  /* Sort the datasets if they aren't sorted. If the dataset is already
     sorted, then 'inplace' is irrelevant. */
  if(sorted_by_first && allf64)
    {
      *A_out=coord1;
      *B_out=coord2;
    }
  else
    {
      /* Allocating a new list is only necessary when 'inplace==0' or all
         the columns are double. */
      if( inplace && allf64 )
        {
          *A_out=coord1;
          *B_out=coord2;
        }
      else
        {
          /* Copy the first list. */
          for(tmp=coord1; tmp!=NULL; tmp=tmp->next)
            {
              c=gal_data_copy(tmp);
              c->next=NULL;
              gal_list_data_add(&A, c);
            }

          /* Copy the second list. */
          for(tmp=coord2; tmp!=NULL; tmp=tmp->next)
            {
              c=gal_data_copy(tmp);
              c->next=NULL;
              gal_list_data_add(&B, c);
            }

          /* Reverse both lists: the copying process reversed the order. */
          gal_list_data_reverse(&A);
          gal_list_data_reverse(&B);

          /* Set the output pointers. */
          *A_out=A;
          *B_out=B;
        }

      /* Sort each dataset by the first coordinate. */
      *A_perm = match_coordinates_prepare_sort(*A_out, minmapsize);
      *B_perm = match_coordinates_prepare_sort(*B_out, minmapsize);
    }
}




















/********************************************************************/
/*************            Coordinate matching           *************/
/********************************************************************/
/* Preparations for 'match_coordinates_second_in_first'. */
static void
match_coordinates_sif_prepare(gal_data_t *A, gal_data_t *B,
                              double *aperture, size_t ndim, double **a,
                              double **b, double *dist, double *c,
                              double *s, int *iscircle)
{
  double semiaxes[3];

  /* These two are common for all dimensions. */
  a[0]=A->array;
  b[0]=B->array;

  /* See if the aperture is a circle or not. */
  switch(ndim)
    {
    case 1:
      *iscircle = 0; /* Irrelevant for 1D. */
      dist[0]=aperture[0];
      break;

    case 2:
      /* Set the main coordinate arrays. */
      a[1]=A->next->array;
      b[1]=B->next->array;

      /* See if the aperture is circular. */
      if( (*iscircle=(aperture[1]==1)?1:0)==0 )
        {
          /* Using the box that encloses the aperture, calculate the
             distance along each axis. */
          gal_box_bound_ellipse_extent(aperture[0], aperture[0]*aperture[1],
                                       aperture[2], dist);

          /* Calculate the sin and cos of the given ellipse if necessary
             for ease of processing later. */
          c[0] = cos( aperture[2] * M_PI/180.0 );
          s[0] = sin( aperture[2] * M_PI/180.0 );
        }
      else
        dist[0]=dist[1]=aperture[0];
      break;

    case 3:
      /* Set the main coordinate arrays. */
      a[1]=A->next->array;
      b[1]=B->next->array;
      a[2]=A->next->next->array;
      b[2]=B->next->next->array;

      if( (*iscircle=(aperture[1]==1 && aperture[2]==1)?1:0)==0 )
        {
          /* Using the box that encloses the aperture, calculate the
             distance along each axis. */
          semiaxes[0]=aperture[0];
          semiaxes[1]=aperture[1]*aperture[0];
          semiaxes[2]=aperture[2]*aperture[0];
          gal_box_bound_ellipsoid_extent(semiaxes, &aperture[3], dist);

          /* Calculate the sin and cos of the given ellipse if necessary
             for ease of processing later. */
          c[0] = cos( aperture[3] * M_PI/180.0 );
          s[0] = sin( aperture[3] * M_PI/180.0 );
          c[1] = cos( aperture[4] * M_PI/180.0 );
          s[1] = sin( aperture[4] * M_PI/180.0 );
          c[2] = cos( aperture[5] * M_PI/180.0 );
          s[2] = sin( aperture[5] * M_PI/180.0 );
        }
      else
        dist[0]=dist[1]=dist[2]=aperture[0];
      break;

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
            "the problem. The value %zu is not recognized for ndim",
            __func__, PACKAGE_BUGREPORT, ndim);
    }
}





static double
match_coordinates_elliptical_r_2d(double d1, double d2, double *ellipse,
                                  double c, double s)
{
  double Xr = d1 * ( c       )     +   d2 * ( s );
  double Yr = d1 * ( -1.0f*s )     +   d2 * ( c );
  return sqrt( Xr*Xr + Yr*Yr/ellipse[1]/ellipse[1] );
}





static double
match_coordinates_elliptical_r_3d(double *delta, double *ellipsoid,
                                  double *c, double *s)
{
  double Xr, Yr, Zr;
  double c1=c[0], s1=s[0];
  double c2=c[1], s2=s[1];
  double c3=c[2], s3=s[2];
  double q1=ellipsoid[1], q2=ellipsoid[2];
  double x=delta[0], y=delta[1], z=delta[2];

  Xr = x*(  c3*c1   - s3*c2*s1 ) + y*( c3*s1   + s3*c2*c1) + z*( s3*s2 );
  Yr = x*( -1*s3*c1 - c3*c2*s1 ) + y*(-1*s3*s1 + c3*c2*c1) + z*( c3*s2 );
  Zr = x*(  s1*s2              ) + y*(-1*s2*c1           ) + z*( c2    );
  return sqrt( Xr*Xr + Yr*Yr/q1/q1 + Zr*Zr/q2/q2 );
}





static double
match_coordinates_distance(double *delta, int iscircle, size_t ndim,
                           double *aperture, double *c, double *s)
{
  /* For more than one dimension, we'll need to calculate the distance from
     the deltas (differences) in each dimension. */
  switch(ndim)
    {
    case 1:
      return fabs(delta[0]);

    case 2:
      return ( iscircle
               ? sqrt( delta[0]*delta[0] + delta[1]*delta[1] )
               : match_coordinates_elliptical_r_2d(delta[0], delta[1],
                                                   aperture, c[0], s[0]) );

    case 3:
      return ( iscircle
               ? sqrt( delta[0]*delta[0]
                       + delta[1]*delta[1]
                       + delta[2]*delta[2] )
               : match_coordinates_elliptical_r_3d(delta, aperture, c, s) );

    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
            "the problem. The value %zu is not recognized for ndim",
            __func__, PACKAGE_BUGREPORT, ndim);
    }

  /* Control should not reach this point. */
  error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix the "
        "problem. Control should not reach the end of this function",
        __func__, PACKAGE_BUGREPORT);
  return NAN;
}





/* Go through both catalogs and find which records/rows in the second
   catalog (catalog b) are within the acceptable distance of each record in
   the first (a). */
static void
match_coordinates_second_in_first(gal_data_t *A, gal_data_t *B,
                                  double *aperture,
                                  struct match_coordinate_sfll **bina)
{
  /* To keep things easy to read, all variables related to catalog 1 start
     with an 'a' and things related to catalog 2 are marked with a 'b'. The
     redundant variables (those that equal a previous value) are only
     defined to make it easy to read the code.*/
  int iscircle=0;
  size_t i, ar=A->size, br=B->size;
  size_t ai, bi, blow=0, prevblow=0;
  size_t ndim=gal_list_data_number(A);
  double r, c[3]={NAN, NAN, NAN}, s[3]={NAN, NAN, NAN};
  double dist[3]={NAN, NAN, NAN}, delta[3]={NAN, NAN, NAN};
  double *a[3]={NULL, NULL, NULL}, *b[3]={NULL, NULL, NULL};


  /* Necessary preperations. */
  match_coordinates_sif_prepare(A, B, aperture,  ndim, a, b, dist, c, s,
                                &iscircle);

  /* For each row/record of catalog 'a', make a list of the nearest records
     in catalog b within the maximum distance. Note that both catalogs are
     sorted by their first axis coordinate.*/
  for(ai=0;ai<ar;++ai)
    if(blow<br)
      {
        /* Initialize 'bina'. */
        bina[ai]=NULL;

        /* Find the first (lowest first axis value) row/record in catalog
           'b' that is within the search radius for this record of catalog
           'a'. 'blow' is the index of the first element to start searching
           in the catalog 'b' for a match to 'a[][ai]' (the record in
           catalog a that is currently being searched). 'blow' is only
           based on the first coordinate, not the second.

           Both catalogs are sorted by their first coordinate, so the
           'blow' to search for the next record in catalog 'a' will be
           larger or equal to that of the previous catalog 'a' record. To
           account for possibly large distances between the records, we do
           a search here to change 'blow' if necessary before doing further
           searching.*/
        for( blow=prevblow; blow<br && b[0][blow] < a[0][ai]-dist[0]; ++blow)
          { /* This can be blank, the 'for' does all we need :-). */ }


        /* 'blow' is now found for this 'ai' and will be used unchanged to
           the end of the loop. So keep its value to help the search for
           the next entry in catalog 'a'. */
        prevblow=blow;


        /* Go through catalog 'b' (starting at 'blow') with a first axis
           value smaller than the maximum acceptable range for 'si'. */
        for( bi=blow; bi<br && b[0][bi] <= a[0][ai] + dist[0]; ++bi )
          {
            /* Only consider records with a second axis value in the
               correct range, note that unlike the first axis, the second
               axis is no longer sorted. so we have to do both lower and
               higher limit checks for each item.

               Positions can have an accuracy to a much higher order of
               magnitude than the search radius. Therefore, it is
               meaning-less to sort the second axis (after having sorted
               the first). In other words, very rarely can two first axis
               coordinates have EXACTLY the same floating point value as
               each other to easily define an independent sorting in the
               second axis. */
            if( ndim<2
                || ( b[1][bi] >= a[1][ai]-dist[1]
                     && b[1][bi] <= a[1][ai]+dist[1] ) )
              {
                /* Now, 'bi' is within the rectangular range of 'ai'. But
                   this is not enough to consider the two objects matched
                   for the following reasons:

                   1) Until now we have avoided calculations other than
                   larger or smaller on double precision floating point
                   variables for efficiency. So the 'bi' is within a square
                   of side 'dist[0]*dist[1]' around 'ai' (not within a
                   fixed radius).

                   2) Other objects in the 'b' catalog may be closer to
                   'ai' than this 'bi'.

                   3) The closest 'bi' to 'ai' might be closer to another
                   catalog 'a' record.

                   To address these problems, we will use a linked list to
                   keep the indexes of the 'b's near 'ai', along with their
                   distance. We only add the 'bi's to this list that are
                   within the acceptable distance.

                   Since we are dealing with much fewer objects at this
                   stage, it is justified to do complex mathematical
                   operations like square root and multiplication. This
                   fixes the first problem.

                   The next two problems will be solved with the list after
                   parsing of the whole catalog is complete.*/
                if( ndim<3
                    || ( b[2][bi] >= a[2][ai]-dist[2]
                         && b[2][bi] <= a[2][ai]+dist[2] ) )
                  {
                    for(i=0;i<ndim;++i) delta[i]=b[i][bi]-a[i][ai];
                    r=match_coordinates_distance(delta, iscircle, ndim,
                                                 aperture, c, s);
                    if(r<aperture[0])
                      match_coordinate_add_to_sfll(&bina[ai], bi, r);
                  }
              }
          }


        /* If there was no objects within the acceptable distance, then the
           linked list pointer will be NULL, so go on to the next 'ai'. */
        if(bina[ai]==NULL)
          continue;

        /* For checking the status of affairs uncomment this block
           {
           struct match_coordinate_sfll *tmp;
           printf("\n\nai: %lu:\n", ai);
           printf("ax: %f (%f -- %f)\n", a[0][ai], a[0][ai]-dist[0],
           a[0][ai]+dist[0]);
           printf("ay: %f (%f -- %f)\n", a[1][ai], a[1][ai]-dist[1],
           a[1][ai]+dist[1]);
           for(tmp=bina[ai];tmp!=NULL;tmp=tmp->next)
           printf("%lu: %f\n", tmp->v, tmp->f);
           }
        */
      }
}





/* In 'match_coordinates_second_in_first', we made an array of lists, here
   we want to reverse that list to fix the second two issues that were
   discussed there. */
void
match_coordinates_rearrange(gal_data_t *A, gal_data_t *B,
                            struct match_coordinate_sfll **bina)
{
  size_t bi;
  float *fp, *fpf, r, *ainb;
  size_t ai, ar=A->size, br=B->size;

  /* Allocate the space for 'ainb' and initialize it to NaN (since zero is
     meaningful) in this context (both for indexs and also for
     floats). This is a two column array that will keep the distance and
     index of the closest element in catalog 'a' for each element in
     catalog b. */
  errno=0; ainb=calloc(2*br, sizeof *ainb);
  if(ainb==NULL)
    error(EXIT_FAILURE, errno, "%s: %zu bytes for 'ainb'", __func__,
          br*sizeof *ainb);
  fpf=(fp=ainb)+2*br; do *fp++=NAN; while(fp<fpf);

  /* Go over each object in catalog 'a' and re-distribute the near objects,
     to find which ones in catalog 'a' are within the search radius of
     catalog b in a sorted manner. Note that we only need the 'ai' with the
     minimum distance to 'bi', the rest are junk.*/
  for( ai=0; ai<ar; ++ai )
    while( bina[ai] )	/* As long as its not NULL.            */
      {
	/* Pop out a 'bi' and its distance to this 'ai' from 'bina'. */
	match_coordinate_pop_from_sfll(&bina[ai], &bi, &r);

	/* If nothing has been put here (the isnan condition below is
	   true), or something exists (the isnan is false, and so it
	   will check the second OR test) with a distance that is
	   larger than this distance then just put this value in. */
	if( isnan(ainb[bi*2]) || r<ainb[bi*2+1] )
	  {
	    ainb[bi*2  ] = ai;
	    ainb[bi*2+1] = r;
	  }
      }

  /* For checking the status of affairs uncomment this block
  {
    printf("\n\nFilled ainb:\n");
    for(bi=0;bi<br;++bi)
      if( !isnan(ainb[bi*2]) )
	printf("bi: %lu: %.0f, %f\n", bi, ainb[bi*2], ainb[bi*2+1]);
  }
  */

  /* Re-fill the bina array, but this time only with the 'bi' that is
     closest to it. Note that bina was fully set to NULL after popping all
     the elements in the loop above.*/
  for( bi=0; bi<br; ++bi )
    if( !isnan(ainb[bi*2]) )
      {
	/* Just to keep the same terminology as before and easier
	   reading.*/
	r=ainb[bi*2+1];
	ai=(size_t)(ainb[bi*2]);

	/* Check if this is the first time we are associating a 'bi' to
	   this 'ai'. If so, then just allocate a single element
	   list. Otherwise, see if the distance is closer or not. If so,
	   replace the values in the single node. */
	if( bina[ai] )
	  {
	    /* If the distance of this record is smaller than the existing
	       entry, then replace the values. */
	    if( r < bina[ai]->f )
	      {
		bina[ai]->f=r;
		bina[ai]->v=bi;
	      }
	  }
	else
          match_coordinate_add_to_sfll(&bina[ai], bi, r);
      }

  /* For checking the status of affairs uncomment this block
  {
    size_t bi, counter=0;
    double *a[2]={A->array, A->next->array};
    double *b[2]={B->array, B->next->array};
    printf("Rearranged bina:\n");
    for(ai=0;ai<ar;++ai)
      if(bina[ai])
        {
          ++counter;
          bi=bina[ai]->v;
          printf("A_%lu (%.8f, %.8f) <--> B_%lu (%.8f, %.8f):\n\t%f\n",
                 ai, a[0][ai], a[1][ai], bi, b[0][bi], b[1][bi],
                 bina[ai]->f);
        }
    printf("\n-----------\nMatched: %zu\n", counter);
  }
  exit(0);
  */

  /* Clean up */
  free(ainb);
}





/* The matching has been done, write the output. */
static gal_data_t *
gal_match_coordinates_output(gal_data_t *A, gal_data_t *B, size_t *A_perm,
                             size_t *B_perm,
                             struct match_coordinate_sfll **bina,
                             size_t minmapsize, int quietmmap)
{
  float r;
  double *rval;
  gal_data_t *out;
  uint8_t *Bmatched;
  size_t ai, bi, nummatched=0;
  size_t *aind, *bind, match_i, nomatch_i;

  /* Find how many matches there were in total */
  for(ai=0;ai<A->size;++ai) if(bina[ai]) ++nummatched;


  /* If there aren't any matches, return NULL. */
  if(nummatched==0) return NULL;


  /* Allocate the output list. */
  out=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &A->size, NULL, 0,
                     minmapsize, quietmmap, "CAT1_ROW", "counter",
                     "Row index in first catalog (counting from 0).");
  out->next=gal_data_alloc(NULL, GAL_TYPE_SIZE_T, 1, &B->size, NULL, 0,
                           minmapsize, quietmmap, "CAT2_ROW", "counter",
                           "Row index in second catalog (counting "
                           "from 0).");
  out->next->next=gal_data_alloc(NULL, GAL_TYPE_FLOAT64, 1, &nummatched,
                                 NULL, 0, minmapsize, quietmmap,
                                 "MATCH_DIST", NULL,
                                 "Distance between the match.");


  /* Allocate the 'Bmatched' array which is a flag for which rows of the
     second catalog were matched. The columns that had a match will get a
     value of one while we are parsing them below. */
  Bmatched=gal_pointer_allocate(GAL_TYPE_UINT8, B->size, 1, __func__,
                                "Bmatched");


  /* Initialize the indexs. We want the first 'nummatched' indexs in both
     outputs to be the matching rows. The non-matched rows should start to
     be indexed after the matched ones. So the first non-matched index is
     at the index 'nummatched'. */
  match_i   = 0;
  nomatch_i = nummatched;


  /* Fill in the output arrays. */
  aind = out->array;
  bind = out->next->array;
  rval = out->next->next->array;
  for(ai=0;ai<A->size;++ai)
    {
      /* A match was found. */
      if(bina[ai])
        {
          /* Note that the permutation keeps the original indexs. */
          match_coordinate_pop_from_sfll(&bina[ai], &bi, &r);
          rval[ match_i   ] = r;
          aind[ match_i   ] = A_perm ? A_perm[ai] : ai;
          bind[ match_i++ ] = B_perm ? B_perm[bi] : bi;

          /* Set a '1' for this object in the second catalog. This will
             later be used to find which rows didn't match to fill in the
             output.. */
          Bmatched[ B_perm ? B_perm[bi] : bi ] = 1;
        }

      /* No match found. At this stage, we can only fill the indexs of the
         first input. The second input needs to be matched afterwards.*/
      else aind[ nomatch_i++ ] = A_perm ? A_perm[ai] : ai;
    }


  /* Complete the second input's permutation. */
  nomatch_i=nummatched;
  for(bi=0;bi<B->size;++bi)
    if( Bmatched[bi] == 0 )
      bind[ nomatch_i++ ] = bi;


  /* For a check
  printf("\nFirst input's permutation:\n");
  for(ai=0;ai<A->size;++ai)
    printf("%s%zu\n", ai<nummatched?"  ":"* ", aind[ai]+1);
  printf("\nSecond input's permutation:\n");
  for(bi=0;bi<B->size;++bi)
    printf("%s%zu\n", bi<nummatched?"  ":"* ", bind[bi]+1);
  exit(0);
  */

  /* Return the output. */
  return out;
}




















/********************************************************************/
/*************            Coordinate matching           *************/
/********************************************************************/
/* Match two positions: the two inputs ('coord1' and 'coord2') should be
   lists of coordinates (each is a list of datasets). To speed up the
   search, this function will sort the inputs by their first column. If
   both are already sorted, give a non-zero value to
   'sorted_by_first'. When sorting is necessary and 'inplace' is non-zero,
   the actual inputs will be sorted. Otherwise, an internal copy of the
   inputs will be made which will be used (sorted) and later
   freed. Therefore when 'inplace==0', the input's won't be changed.

   IMPORTANT NOTE: the output permutations will correspond to the initial
   inputs. Therefore, even when 'inplace' is non-zero (and this function
   changes the inputs' order), the output permutation will correspond to
   original inputs.

   The output is a list of 'gal_data_t' with the following columns:

       Node 1: First catalog index (counting from zero).
       Node 2: Second catalog index (counting from zero).
       Node 3: Distance between the match.                    */
gal_data_t *
gal_match_coordinates(gal_data_t *coord1, gal_data_t *coord2,
                      double *aperture, int sorted_by_first,
                      int inplace, size_t minmapsize, int quietmmap,
                      size_t *nummatched)
{
  int allf64=1;
  gal_data_t *A, *B, *out;
  size_t *A_perm=NULL, *B_perm=NULL;
  struct match_coordinate_sfll **bina;

  /* Do a small sanity check and make the preparations. After this point,
     we'll call the two arrays 'a' and 'b'.*/
  match_coordinaes_sanity_check(coord1, coord2, aperture, inplace,
                                &allf64);
  match_coordinates_prepare(coord1, coord2, sorted_by_first, inplace, allf64,
                            &A, &B, &A_perm, &B_perm, minmapsize);


  /* Allocate the 'bina' array (an array of lists). Let's call the first
     catalog 'a' and the second 'b'. This array has 'a->size' elements
     (pointers) and for each, it keeps a list of 'b' elements that are
     nearest to it. */
  errno=0;
  bina=calloc(A->size, sizeof *bina);
  if(bina==NULL)
    error(EXIT_FAILURE, errno, "%s: %zu bytes for 'bina'", __func__,
          A->size*sizeof *bina);


  /* All records in 'b' that match each 'a' (possibly duplicate). */
  match_coordinates_second_in_first(A, B, aperture, bina);


  /* Two re-arrangings will fix the issue. */
  match_coordinates_rearrange(A, B, bina);


  /* The match is done, write the output. */
  out=gal_match_coordinates_output(A, B, A_perm, B_perm, bina, minmapsize,
                                   quietmmap);


  /* Clean up. */
  free(bina);
  if(A!=coord1)
    {
      gal_list_data_free(A);
      gal_list_data_free(B);
    }
  if(A_perm) free(A_perm);
  if(B_perm) free(B_perm);


  /* Set 'nummatched' and return output. */
  *nummatched = out ?  out->next->next->size : 0;
  return out;
}
