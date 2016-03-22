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
brightnessmag(struct mkcatalogparams *p, int m0b1, int cinobj, int isriver)
{
  size_t i, col;
  double bright, *value;
  char *type, *mb=m0b1 ? "brightness" :"magnitude";

  /* Set the proper column to use */
  if(p->obj0clump1) /* It is a clump.                            */
    {
      if(cinobj) { --p->curcol; return; } /* cinobj is only for objects. */
      if(isriver)
        {
          if(m0b1==0) { --p->curcol; return; }
          type="Rivers around this clump, average";
          col=CRivAve;
        }
      else
        {
          col  = CBrightness;
          type = "This clump";
        }
    }
  else
    {               /* It is an object.                         */
      if(isriver) { --p->curcol; return; }
      if(cinobj)    /* It is the positions of clumps in object. */
        {
          col  = OBrightnessC;
          type = "Clumps in object";
        }
      else          /* It is the position of the object itsself.*/
        {
          col  = OBrightness;
          type = "Full object";
        }
    }

  /* For the comments: */
  p->unitp = ( m0b1 ?
               (isriver==1 ? CATUNITAVE : CATUNITBRIGHTNESS)
               : CATUNITMAG );
  sprintf(p->description, "%lu: %s %s.",
          p->curcol, type, mb);

  /* Fill the column: */
  for(i=0;i<p->num;++i)
    {
      /* Set the basic values: */
      bright = p->info[ (i+1) * p->icols + col ];
      value  = &p->cat[i * p->numcols + p->curcol ];

      /* If we are dealing with a clump, then you have to subtract the
         average river flux multiplied by the the area of the
         clump. The value in the CBrightness column is simply the sum
         of pixels. */
      if(p->obj0clump1 && isriver==0)
        bright -= ( p->info[ (i+1) * p->icols + CRivAve ]
                    * p->info[ (i+1) * p->icols + CAREA ] );

      /* Do the job: */
      switch(m0b1)
        {
        case 0:
          *value = bright<=0.0f ? NAN : -2.5f*log10(bright)+p->zeropoint;
          break;
        case 1:
          *value = bright;
          break;
        default:
          error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can "
                "fix this issue. For some reason, the value to m0b1 in"
                "brightnessfluxmag (columns.c) is %d, which is not "
                "recognized.", PACKAGE_BUGREPORT, m0b1);
        }
    }

  /* When dealing with brightness set the column to accurate mode: */
  if(isriver)
    p->accucols[p->accucounter++]=p->curcol;
}





void
skystd(struct mkcatalogparams *p, int issky)
{
  size_t i, scol, acol;

  /* Set the proper column to use */
  acol = p->obj0clump1 ? CAREA : OAREA;
  scol = p->obj0clump1 ? (issky ? CSKY : CSTD) : (issky ? OSKY : OSTD);

  /* For the comments: */
  p->unitp = CATUNITAVE;
  sprintf(p->description, "%lu: Average %s over this %s.",
          p->curcol, issky ? "sky" : "sky standard deviation", p->name);

  /* Fill the sky value, note that in the information array, we have
     only calculated the sum. So here, we need to divide by the area
     to find the average. */
  for(i=0;i<p->num;++i)
    p->cat[i * p->numcols + p->curcol ] =
      p->info[(i+1)*p->icols+scol]/p->info[(i+1)*p->icols+acol];

  /* This column should be accurate: */
  p->accucols[p->accucounter++]=p->curcol;
}





void
sncol(struct mkcatalogparams *p)
{
  size_t i;
  double I, O, Ni, err;
  size_t stdcol        = p->obj0clump1 ? CSTD        : OSTD;
  size_t areacol       = p->obj0clump1 ? CAREA       : OAREA;
  size_t brightnesscol = p->obj0clump1 ? CBrightness : OBrightness;

  /* For the comments: */
  p->unitp = CATUNITRATIO;
  sprintf(p->description, "%lu: Signal to noise ratio.", p->curcol);

  /* Calculate the signal to noise ratio. Recall that for the objects,
     the sky value was subtracted from oinfo, but for the clumps, it
     was not subtracted. */
  for(i=0;i<p->num;++i)
    {
      Ni=p->info[(i+1)*p->icols+areacol];
      I=p->info[(i+1)*p->icols+brightnesscol]/Ni;
      err=p->info[(i+1)*p->icols+stdcol]*p->info[(i+1)*p->icols+stdcol];

      if(p->obj0clump1)
        {
          /* Set the "outer" flux (average river flux or sky
             background flux if there were no river pixels.) */
          O=p->cinfo[ (i+1) * CCOLUMNS + CRivAve ];

          /* Note that the sky error is not needed when there were
             rivers. However, if the grownclumps were used, full
             detections can be considered as a clump and there will be
             no rivers. In such cases, the average river flux is set
             to the sky value in the secondpass function
             (mkcatalog.). So here we need to use the sky error as
             it was calculated above and leave err untouched. */
          if(p->cinfo[ (i+1) * CCOLUMNS + CRivArea ] > 0.0f)
            err = (O>0?O:-1*O) + err*(p->skysubtracted ? 2.0f : 0.0f);

          p->cat[i * p->numcols + p->curcol ] =
            ( sqrt(Ni/p->cpscorr)*(I-O) / sqrt( (I>0?I:-1*I) + err ) );
        }
      else
        {
          err *= p->skysubtracted ? 1.0f : 2.0f;

          /* Note that `I' was sky subtracted for objects in firstpass
             (mkcatalog.) */
          p->cat[i * p->numcols + p->curcol ] =
            sqrt( Ni/p->cpscorr ) * I / sqrt( (I>0?I:-1*I) + err );
        }
    }

}
