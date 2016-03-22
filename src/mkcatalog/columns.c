/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

ABOUT THIS FILE: The information tables are fully explained in the
  comments of main.h. After the raw information is read in the first
  and second pass, the job of the functions here is to process the raw
  columns that are needed into useful knowledge and print them. for
  example those functions will only record the weighted sum of pixel
  positions and the total weight, here the weighted sum is divided by
  the total weight to yeild an average.

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
#include <string.h>
#include <stdlib.h>


#include "main.h"

#include "columns.h"
#include "mkcatalog.h"










/******************************************************************/
/*******        Information table modifications       *************/
/******************************************************************/
/* Correct the average Sky and Sky standard deviation value for
   objects and clumps. Note that during the passes, these were just
   sums of pixel values, they need to be divided by the area of the
   object/clump, which is done here. */
void
setskystd(struct mkcatalogparams *p, size_t col)
{
  size_t ac;
  double *row = p->info + p->icols;
  double *end = row + (p->icols * p->num);

  /* Only do the correction if this column is not already flagged. */
  if(p->info[col]==0.0f)
    {

      /* Set the area column: */
      ac = p->obj0clump1 ? CAREA : OAREA;

      /* Go over every row and do the correction. */
      do
        {
          row[col] = ( row[ac]>0.0f ? row[col]/row[ac] : NAN );
          row+=p->icols;
        }
      while(row<end);

      /* Set the flag so this operation is not done again. */
      p->info[col]=1.0f;
    }
}





/* Correct the average river value, after the passes, it is just the
   sum. */
void
setaveriver(struct mkcatalogparams *p)
{
  double *row = p->info + p->icols;
  double *end = row + (p->icols * p->num);

  /* Only do the correction if this column is not already flagged. */
  if(p->info[CRivAve]==0.0f)
    {

      /* Make sure the Sky values are corrected */
      setskystd(p, CSKY);

      /* Go over every row and do the correction. Note that in cases
         where the grown clumps are used instead of the clumps, we are
         not going to have any rivers (row[CRivArea]==0.0f). In such
         situations, set the per-pixel average river value to the Sky
         value under the clump. The reason is that for the clumps, Sky
         subtraction was not done on the Clump brightness, so this
         value will be used, and if there was no river, then we need
         something to replace it. */
      do
        {
          row[CRivAve] = ( row[CRivArea]>0.0f
                           ? row[CRivAve]/row[CRivArea] : CSKY );
          row+=p->icols;
        }
      while(row<end);

      /* Set the flag so this operation is not done again. */
      p->info[CRivAve]=1.0f;
    }
}





/* The clump brightness values are not Sky subtracted since the river
   values (which are also not Sky subtracted) should be subtracted
   from them. Here that job is done. */
void
setclumpbrightness(struct mkcatalogparams *p)
{
  double *row = p->info + p->icols;
  double *end = row + (p->icols * p->num);

  /* Only do the correction if this column is not already flagged. */
  if(p->info[CBrightness]==0.0f)
    {

      /* Make sure the average river value is calculated */
      setaveriver(p);

      /* On a clump, we have to subtract the average river flux
         multiplied by the the area of the clump. The value in the
         CBrightness column is simply the sum of pixels. Note that
         here we are multiplying by the area of the clump (CAREA),
         while in setaveriver(), we divided by the area of the river
         (CRivArea). */
      do
        {
          row[CBrightness] -= row[CRivAve]*row[CAREA];
          row+=p->icols;
        }
      while(row<end);

      /* Set the flag so this operation is not done again. */
      p->info[CBrightness]=1.0f;
    }
}





/* Find the geometric center of the profile (average position,
   ignoring any flux of the pixels). */
void
geoxy(struct mkcatalogparams *p, size_t col)
{
  size_t ac=-1;
  double *row = p->info + p->icols;
  double *end = row + (p->icols * p->num);

  /* Only if this column is not flagged as already done (==1.0f). */
  if(p->info[col]==0.0f)
    {

      /* First, set the columns to use for the conversion. */
      if(p->obj0clump1)                         ac=CAREA;
      else
        {
          if      (col==OGeoX || col==OGeoY)    ac=OAREA;
          else if (col==OGeoCX || col==OGeoCY)  ac=OAREAC;
          else
            error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we "
                  "can fix this. The given column in the --OBJECTS-- "
                  "information table was not recognized for calculating the "
                  "geometric X and/or Y.", PACKAGE_BUGREPORT);
        }

      /* Go over all the rows and correct this column. */
      do
        {
          /* Set the value for this row. Note that unlike the
             calculations here that start counting with zero, the FITS
             standard starts counting from 1, so add a one after
             dividing by the area. If the area is zero, then set
             NaN. */
          row[col] = row[ac]>0.0f ? row[col] / row[ac] + 1 : NAN;

          /* Go onto the next row: */
          row+=p->icols;
        }
      while(row<end);

      /* Flag this column as complete for future reference. */
      p->info[col]=1.0f;
    }
}





/* A low-level function used to find the flux weighted center, since
   it is needed by multiple columns. The geometric center for this
   axis colum (geocol) and area column (areacol) are needed for backup
   (when there might not be any positive flux pixel/data values to use
   for weight). */
void
flxwhtimg(struct mkcatalogparams *p, size_t col)
{
  size_t wc=-1, gc=-1;
  double *row = p->info + p->icols;
  double *end = row + (p->icols * p->num);;


  /* Only if this column is not flagged as already done (==1.0f). */
  if(p->info[col]==0.0f)
    {

      /* First, set the columns to use for the conversion. */
      if(p->obj0clump1)
        {
          wc=CPosBright;
          if     (col==CFlxWhtX) gc=CGeoX;
          else if(col==CFlxWhtY) gc=CGeoY;
          else
            error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we "
                  "can fix this. The given column in the --CLUMPS-- "
                  "information table was not recognized for calculating the "
                  "final flux weighted X and/or Y.", PACKAGE_BUGREPORT);
        }
      else
        {
          if (col==OFlxWhtX || col==OFlxWhtY)
            {
              wc=OPosBright;
              gc = col==OFlxWhtX ? OGeoX : OGeoY;
            }
          else if (col==OFlxWhtCX || col==OFlxWhtCY)
            {
              wc=OPosBrightC;
              gc = col==OFlxWhtCX ? OGeoCX : OGeoCY;
            }
          else
            error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we "
                  "can fix this. The given column in the --OBJECTS-- "
                  "information table was not recognized for calculating the "
                  "final flux weighted X and/or Y.", PACKAGE_BUGREPORT);
        }


      /* The geometric positions act as a backup for the flux weighted
         centers, so make sure the appropriate geometric center is
         defined. */
      geoxy(p, gc);

      /* For a check, uncomment these two lines:
      printf("\n\nInfocol: %lu (%s, %lu)\n", col,
             p->info==p->oinfo?"Objects":"Clumps", p->num);
      Then add these two lines before and after row[col] in the loop*/
      /*printf("%lu: %f --> ", (row-p->info)/p->icols, row[col]);*/
      /*printf("%f\n", row[col]);*/

      /* Go over all the rows and correct this column. */
      do
        {
          /* Set the value for this row. When a positive weight is
             present, we are adding with one (1) because of the FITS
             standard. */
          row[col] = row[wc]>0.0f ? (row[col]/row[wc])+1 : row[gc];

          /* Go onto the next row: */
          row+=p->icols;
        }
      while(row<end);

      /* Set the flag for this column to one, so this whole proces is
         not done again. */
      p->info[col]=1.0f;
    }
}





/* Fill in the RA and Dec columns, note that we will need the X and Y
   colums first for this. */
void
flxwhtwcs(struct mkcatalogparams *p, size_t col)
{
  /* Initialize all the columns to -1 (largest possible number in the
     C's unsigned char, so if there is any bugs, we get a good error. */
  size_t xc=-1, yc=-1, rc=-1, dc=-1;


  /* RA and Dec are usually needed together and must also be
     calculated together, but here, we are giving the user complete
     freedom in setting the columns in which ever order they want. So
     after calculating the RA and Dec once for either of the two,
     there is no more need to do the calculation again.  */
  if(p->info[col]==0.0f)
    {

      /* First, set the columns to use for the conversion. */
      if(p->obj0clump1)
        {
          if(col==CFlxWhtRA || col==CFlxWhtDec)
            {
              xc=CFlxWhtX;  yc=CFlxWhtY;
              rc=CFlxWhtRA; dc=CFlxWhtDec;
            }
          else
            error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we "
                  "can fix this. The given column in the --CLUMPS-- "
                  "information table was not recognized for calculating the "
                  "RA and Dec.", PACKAGE_BUGREPORT);
        }
      else
        {
          if(col==OFlxWhtCRA || col==OFlxWhtCDec)
            {
              xc=OFlxWhtCX;  yc=OFlxWhtCY;
              rc=OFlxWhtCRA; dc=OFlxWhtCDec;
            }
          else if(col==OFlxWhtRA || col==OFlxWhtDec)
            {
              xc=OFlxWhtX;  yc=OFlxWhtY;
              rc=OFlxWhtRA; dc=OFlxWhtDec;
            }
          else
            error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we "
                  "can fix this. The given column in the --OBJECT-- "
                  "information table was not recognized for calculating the "
                  "RA and Dec.", PACKAGE_BUGREPORT);
        }


      /* Finalize the relevant X and Y positions first (which are
         needed for the WCS conversion). Note that if they are ready
         to use (their flag is 1.0f), flxwhtxy will not do
         anything. But if the user hasn't already asked for X and Y,
         then these columns will be corrected here.*/
      flxwhtimg(p, xc);
      flxwhtimg(p, yc);


      /* Do the conversion. Note that the p->icols is added because
         the first row is not used by any object or colump (since
         their indexes begin from 1).*/
      xyarraytoradec(p->wcs, p->info+p->icols+xc, p->info+p->icols+rc,
                     p->num, p->icols);


      /* Set the flag of the converted columns to 1.0f, so the
         calculations are not repeated if any of the columns is needed
         again. Note that it is irrelevant which one of the RA or Dec
         were calculated, so we are not using `col' here. */
      p->info[rc]=p->info[dc]=1.0f;
    }
}




















/******************************************************************/
/***************    Add columns for printing    *******************/
/******************************************************************/
void
idcol(struct mkcatalogparams *p)
{
  size_t i;

  p->unitp=CATUNITCOUNTER;
  sprintf(p->description, "%lu: Overall %s ID",
          p->curcol, p->name);

  for(i=0;i<p->num;++i)
    p->cat[ i * p->numcols + p->curcol ] = i+1;

  p->intcols[p->intcounter++]=p->curcol;
}





/* Store IDs related to the host object:

   o1c0==1 --> hostobjid: The ID of object hosting this clump
   o1c0==0 --> idinhostobj: The ID of clump in object
 */
void
hostobj(struct mkcatalogparams *p, int o1c0)
{
  char *des;
  double counter;
  size_t i, j, n, row=0;

  p->unitp=CATUNITCOUNTER;
  des = ( o1c0 ? "ID of object hosting this clump"
          : "ID of clump in host object" );
  sprintf(p->description, "%lu: %s.",
          p->curcol, des);

  for(i=1;i<=p->numobjects;++i)
    if( (n=p->oinfo[i*OCOLUMNS+ONCLUMPS]) > 0.0f)
      {
        counter=1.0f;
        for(j=row;j<row+n;++j)
          p->cat[ j * p->numcols + p->curcol] = o1c0 ? i : counter++;
        row+=n;
      }

  p->intcols[p->intcounter++]=p->curcol;
}





void
numclumps(struct mkcatalogparams *p)
{
  size_t i;

  p->unitp=CATUNITCOUNTER;
  sprintf(p->description, "%lu: Number of clumps in this object.",
          p->curcol);

  for(i=0;i<p->numobjects;++i)
    p->cat[i * p->numcols + p->curcol ] = p->oinfo[(i+1)*OCOLUMNS+ONCLUMPS];

  p->intcols[p->intcounter++]=p->curcol;
}





void
area(struct mkcatalogparams *p, int cinobj, int isriver)
{
  char *type;
  size_t i, col;

  /* Set the proper column to use */
  if(p->obj0clump1)
    {
      if(isriver)
        {
          type="Number of river pixels around this clump";
          col=CRivArea;
        }
      else
        {
          type="Area of this clump";
          col = CAREA;
        }
    }
  else
    {
      if(cinobj)    /* It is the positions of clumps in object. */
        {
          type="Clumps in object area";
          col = OAREAC;
        }
      else          /* It is the position of the object itsself.*/
        {
          type="Full object area";
          col = OAREA;
        }
    }

  /* Set the unit and print the header. */
  p->unitp = isriver ? CATUNITCOUNTER : CATUNITPIXAREA;
  sprintf(p->description, "%lu: %s.", p->curcol, type);

  /* Fill in the output array. */
  for(i=0;i<p->num;++i)
    p->cat[i * p->numcols + p->curcol ] = p->info[(i+1)*p->icols+col];

  /* Set the precision for printing. */
  p->intcols[p->intcounter++]=p->curcol;
}





void
position(struct mkcatalogparams *p, size_t col, char *target,
         char *type, char *axis)
{
  size_t i;
  int wcsax = ( !strcmp(axis, MKCATRA) || !strcmp(axis, MKCATDEC) ) ? 1 : 0;

  /* Set the header information. */
  sprintf(p->description, "%lu: %s %s (%s).", p->curcol, target, type, axis);

  /* Prepare the respective column, set the units and also the
     printing accuracy if we are in RA/Dec mode (wcsax==1). */
  if(wcsax)
    {
      flxwhtwcs(p, col);
      p->unitp = CATUNITDEGREE;
      p->accucols[p->accucounter++]=p->curcol;
    }
  else
    {
      flxwhtimg(p, col);
      p->unitp = CATUNITPIXPOS;
    }

  /* Write respective column of the information table into the output. */
  for(i=0;i<p->num;++i)
    p->cat[i * p->numcols + p->curcol ] = p->info[(i+1)*p->icols+col];
}





void
brightnessmag(struct mkcatalogparams *p, size_t col, char *target,
              char *scale)
{
  size_t i;
  double bright, *value;

  /* Prepare other necessary columns */
  if( !strcmp(MKRIVERSSUR, target) )
    setaveriver(p);
  if( !strcmp(MKCATCLUMP, target) )
    setclumpbrightness(p);

  /* Fill the output columns: */
  for(i=0;i<p->num;++i)
    {

      /* Set the basic values: */
      bright = p->info[ (i+1) * p->icols + col ];
      value  = &p->cat[i * p->numcols + p->curcol ];


      /* Do the job: */
      if(!strcmp(MKCATMAG, scale))
        *value = bright<=0.0f ? NAN : -2.5f*log10(bright)+p->zeropoint;
      else if(!strcmp(MKCATBRIGHT, scale))
        *value = bright;
      else
        error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can "
              "fix this issue. For some reason, the value to `scale' in"
              "brightnessfluxmag (columns.c) is `%s', which is not "
              "recognized.", PACKAGE_BUGREPORT, scale);
    }

  /* Make final preparations for output. When dealing with the average
     river value, set the accuracy to high, also set the units to
     average values (per pixel). */
  if(!strcmp(MKRIVERSSUR, target))
    {
      p->unitp = CATUNITAVE;
      p->accucols[p->accucounter++]=p->curcol;
    }
  else
    p->unitp = strcmp(MKCATMAG, scale) ? CATUNITBRIGHTNESS : CATUNITMAG;
  sprintf(p->description, "%lu: %s %s.", p->curcol, target, scale);
}





void
skystd(struct mkcatalogparams *p, size_t col)
{
  size_t i;

  /* For the comments: */
  p->unitp = CATUNITAVE;
  sprintf(p->description, "%lu: Average %s under this %s.",
          p->curcol, ( (col==OSKY || col==CSKY)
                       ? "sky" : "sky standard deviation" ),
          p->name);

  /* Correct the raw values (divide them by area) if not already
     done. */
  setskystd(p, col);

  /* Fill the sky value, note that in the information array, we have
     only calculated the sum. So here, we need to divide by the area
     to find the average. */
  for(i=0;i<p->num;++i)
    p->cat[i * p->numcols + p->curcol ] = p->info[(i+1)*p->icols+col];

  /* This column should be accurate: */
  p->accucols[p->accucounter++]=p->curcol;
}





void
sncol(struct mkcatalogparams *p)
{
  size_t i;
  double I, O, Ni, errpt, *row;
  size_t stdcol        = p->obj0clump1 ? CSTD        : OSTD;
  size_t areacol       = p->obj0clump1 ? CAREA       : OAREA;
  size_t brightnesscol = p->obj0clump1 ? CBrightness : OBrightness;

  /* Do the corrections:

       1. If we are dealing with clumps, make sure the clump
          brightness is corrected first.

       2. Make sure the STD values are corrected in any case. */
  setskystd(p, stdcol);
  if(p->obj0clump1)
    setclumpbrightness(p);

  /* For the comments: */
  p->unitp = CATUNITRATIO;
  sprintf(p->description, "%lu: Signal to noise ratio.", p->curcol);

  /* Calculate the signal to noise ratio. Recall that for the objects,
     the sky value was subtracted from oinfo, but for the clumps, it
     was not subtracted. */
  for(i=0;i<p->num;++i)
    {
      /* Some convenience variables to make things readable. */
      row = p->info + ((i+1)*p->icols);     /* Pointer to this row.       */
      Ni  = row[ areacol ];                 /* Number-in                  */
      I   = row[ brightnesscol ]/Ni;        /* Inner brightness (average) */
      errpt = row[ stdcol ]*row[ stdcol ];  /* error-to-power-two         */

      /* If we are on a clump and there are actually rivers (NOTE: it
         is possible that there are no rivers, see the NoiseChisel
         dropout paper). In short, they are actually objects with no
         more than one clump. So, NoiseChisel parameters were set such
         that the objects also show up in the clumps labels. */
      if(p->obj0clump1 && row[ CRivArea ]>0.0f)
        {

          /* Another convenience variable. */
          O=row[ CRivAve ];                 /* Outer brightness (average)*/

          /* Modify the error based on the conditions. Note that the
             inner flux has already been subtracted from the average
             outer flux multiplied by the clump area in
             setclumpbrightness and was divided by the clump area
             above. It is also in per pixel units. row[CRivAve] is
             also in per pixel units. So simply by adding the two, we
             get the per pixel flux within the clump before removing
             the average river value.

             If the image was already Sky subtracted, then the Sky
             error^2 (=err) must be multiplied by 2 (we have
             implicitly used it in both estimating the inner and outer
             fluxes). Otherwise, it is multiplied by 0.0f, since we
             don't care because we are not using the Sky value
             here. */
          errpt = ( (   I+O > 0.0f ? I+O : 0.0f )
                    + ( O   > 0.0f ? O   : 0.0f )
                    + errpt * (p->skysubtracted ? 2.0f : 0.0f) );
        }
      else
        {
          /* When the flux is negative (can easily happen in matched
             photometry), then ignore the error in flux (the S/N is
             meaningless anyway) and just keep the Sky error.

             When the image was already Sky subtracted, we need two
             errpt terms, because the error in the previous Sky
             subtraction must also be included. */
          errpt = ( ( I>0 ? I : 0 )
                    + errpt * (p->skysubtracted ? 1.0f : 2.0f) );
        }

      /* Fill in the output column */
      p->cat[i * p->numcols + p->curcol ] =
        ( sqrt(Ni/p->cpscorr)*I / sqrt( errpt ) );
    }

}
