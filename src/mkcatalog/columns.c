/*********************************************************************
MakeCatalog - Make a catalog from an input and labeled image.
MakeCatalog is part of GNU Astronomy Utilities (Gnuastro) package.

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


#include "main.h"

#include "columns.h"
#include "mkcatalog.h"










/* MakeCatalogs and these column functions:
   ========================================

   In order to print the information in any of the desired columns
   this is the procedure we have taken:

   There is a base array for information on objects (p->oinfo) and
   clumps (p->cinfo) which is filled in mkcatalog.c's first and second
   pass functions. Most of the necessary information can be found in
   two passes and if the desired calculation doesn't exist, adding it
   is very easy: just add to the macros of main.h and include the
   calculation in any of the passes.

   These information arrays keep the raw data that is then converted
   into the output catalog through the functions in this file. */










/******************************************************************/
/***************      Operations on columns     *******************/
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

  /* This function is only for clumps. */
  if(p->obj0clump1==0) { --p->curcol; return; }

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

  /* This function is only for objects. */
  if(p->obj0clump1) return;

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
  if(p->obj0clump1) /* It is a clump.                            */
    {
      if(cinobj) { --p->curcol; return; } /* cinobj is only for objects. */
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
    {               /* It is an object.                         */
      if(isriver) { --p->curcol; return; }
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

  p->unitp = isriver ? CATUNITCOUNTER : CATUNITPIXAREA;
  sprintf(p->description, "%lu: %s.", p->curcol, type);

  for(i=0;i<p->num;++i)
    p->cat[i * p->numcols + p->curcol ] = p->info[(i+1)*p->icols+col];

  p->intcols[p->intcounter++]=p->curcol;
}





/* Find the positions:

   x1y1==1: X axis.

   cinobj==1: (only relevant when obj0clump1==0): clump in obj info,
              not full object.

   w1i0==1: We are dealing with world coordinates not image
            coordinates. In the world coordinates x: ra and y: dec
 */
void
position(struct mkcatalogparams *p, int w1i0, int x1y0, int cinobj)
{
  char *type;
  size_t i, pcol;

  /* Set the proper column to use */
  if(p->obj0clump1) /* It is a clump.                            */
    {
      if(cinobj) { --p->curcol; return; } /* cinobj is only for objects. */
      type="This clump";
      if(w1i0)   pcol = x1y0 ? CFlxWhtRA : CFlxWhtDec;
      else       pcol = x1y0 ? CFlxWhtX : CFlxWhtY;
    }
  else
    {               /* It is an object.                         */
      if(cinobj)    /* It is the positions of clumps in object. */
        {
          type="Clumps in object";
          if(w1i0) pcol = x1y0 ? OFlxWhtCRA : OFlxWhtCDec;
          else     pcol = x1y0 ? OFlxWhtCX : OFlxWhtCY;
        }
      else          /* It is the position of the object itsself.*/
        {
          type="Full object";
          if(w1i0) pcol = x1y0 ? OFlxWhtRA : OFlxWhtDec;
          else     pcol = x1y0 ? OFlxWhtX : OFlxWhtY;
        }
    }


  p->unitp = w1i0 ? CATUNITDEGREE : CATUNITPIXPOS;
  sprintf(p->description, "%lu: %s flux weighted center (%s).",
          p->curcol, type,
          w1i0 ? (x1y0?"RA":"Dec") : (x1y0?"X":"Y") );


  for(i=0;i<p->num;++i)
    p->cat[i * p->numcols + p->curcol ] = p->info[(i+1)*p->icols+pcol];


  if(w1i0)
    p->accucols[p->accucounter++]=p->curcol;
}





void
brightnessfluxmag(struct mkcatalogparams *p, int m0b1f2, int cinobj,
                  int isriver)
{
  double bright, *value;
  size_t i, col, acol=(size_t)(-1);
  char *type, *fbm=m0b1f2 ? (m0b1f2==1 ? "brightness" : "flux") :"magnitude";

  /* Set the proper column to use */
  if(p->obj0clump1) /* It is a clump.                            */
    {
      if(cinobj) { --p->curcol; return; } /* cinobj is only for objects. */
      if(isriver)
        {
          if(m0b1f2==0 || m0b1f2==1) { --p->curcol; return; }
          type="Rivers around this clump, average";
          col=CAveRivFlux;
        }
      else
        {
          acol = CAREA;
          col  = CBrightness;
          type = "This clump";
        }
    }
  else
    {               /* It is an object.                         */
      if(isriver) { --p->curcol; return; }
      if(cinobj)    /* It is the positions of clumps in object. */
        {
          acol = OAREAC;
          col  = OBrightnessC;
          type = "Clumps in object";
        }
      else          /* It is the position of the object itsself.*/
        {
          acol = OAREA;
          col  = OBrightness;
          type = "Full object";
        }
    }

  /* For the comments: */
  p->unitp = ( m0b1f2 ?
               (m0b1f2==1 ? CATUNITBRIGHTNESS : CATUNITFLUX)
               : CATUNITMAG );
  sprintf(p->description, "%lu: %s %s.",
          p->curcol, type, fbm);

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
        bright -= ( p->info[ (i+1) * p->icols + CAveRivFlux ]
                    * p->info[ (i+1) * p->icols + CAREA ] );

      /* Do the job: */
      switch(m0b1f2)
        {
        case 0:
          *value = bright<=0.0f ? NAN : -2.5f*log10(bright)+p->zeropoint;
          break;
        case 1:
          *value = bright;
          break;
        case 2:
          *value = bright / ( isriver ? 1.0f : p->info[(i+1)*p->icols+acol] );
          break;
        default:
          error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can "
                "fix this issue. For some reason, the value to m0b1f2 in "
                "brightnessfluxmag (columns.c) is %d, which is not "
                "recognized.", PACKAGE_BUGREPORT, m0b1f2);
        }
    }

  /* When dealing with flux set the column to accurate mode: */
  if(m0b1f2==2)
    p->accucols[p->accucounter++]=p->curcol;
}





void
skystd(struct mkcatalogparams *p, int issky)
{
  size_t i, col;

  /* Set the proper column to use */
  col = p->obj0clump1 ? (issky ? CSKY : CSTD) : (issky ? OSKY : OSTD);

  /* For the comments: */
  p->unitp = CATUNITFLUX;
  sprintf(p->description, "%lu: Average %s over this %s.",
          p->curcol, issky ? "sky" : "sky standard deviation", p->name);

  /* Fill the column: */
  for(i=0;i<p->num;++i)
    p->cat[i * p->numcols + p->curcol ] = p->info[(i+1)*p->icols+col];

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
      err=p->info[(i+1)*p->icols+stdcol];
      Ni=p->info[(i+1)*p->icols+areacol];
      I=p->info[(i+1)*p->icols+brightnesscol]/Ni;

      if(p->obj0clump1)
        {
          err *= p->skysubtracted ? 2.0f*err : 0.0f;
          O=p->cinfo[ (i+1) * CCOLUMNS + CAveRivFlux ];

          p->cat[i * p->numcols + p->curcol ] =
            ( sqrt(Ni/p->cpscorr)*(I-O)
              / sqrt( (I>0?I:-1*I) + (O>0?O:-1*O) + err ) );
        }
      else
        {
          err *= p->skysubtracted ? err : 2.0f*err;
          p->cat[i * p->numcols + p->curcol ] =
            sqrt( Ni/p->cpscorr ) * I / sqrt(I+err);
        }
    }

}
