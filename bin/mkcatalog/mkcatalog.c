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
along with Gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <config.h>

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/txtarray.h>

#include <timing.h>
#include <neighbors.h>

#include "main.h"

#include "columns.h"
#include "mkcatalog.h"





/*********************************************************************/
/*****************     Fill information tables     *******************/
/*********************************************************************/

/* Macro to see if the label is indexable (belongs to an object or
   not). See the explanation in bin/noisechisel/label.h. */
#if GAL_FITS_LONG_BLANK<0
#define ISINDEXABLEOBJLABEL (objects[i]>0)
#define ISINDEXABLECLPLABEL (clumps[i]>0)
#else
#define ISINDEXABLEOBJLABEL (objects[i] && objects[i]!=GAL_FITS_LONG_BLANK)
#define ISINDEXABLECLPLABEL (clumps[i] && clumps[i]!=GAL_FITS_LONG_BLANK)
#endif





/* In the first pass, the most basic properties (mainly about the
   objects) are found. The primary reason is that we still don't know
   how many objects there are in order to be able to put the clump
   information in the proper place. This could maybe be fixed with a
   linked list in one pass, but that would drastically bring down the
   speed. */
void
firstpass(struct mkcatalogparams *p)
{
  float imgss;
  size_t i, is1=p->s1;
  double x, y, sx, sy, *thisobj;
  long *objects=p->objects, *clumps=p->clumps;

  /* Go over the pixels and fill the information array. */
  for(i=0;i<p->s0*p->s1;++i)
    if(ISINDEXABLEOBJLABEL)
      {
        /* thisobj is a pointer to the start of the row in the object
           information array (oinfo). It is mainly used to keep things
           short, simple, less bugy and most importantly: elegant. */
        imgss = p->img[i] - p->sky[i];
        thisobj = p->oinfo + objects[i]*OCOLUMNS;

        /* Set the shifts to avoid round-off errors in large numbers
           for the non-linear calculations. We are using the first
           pixel of each object as the shift parameter to keep the
           mean reasonably near to the standard deviation. Otherwise,
           when the object is far out in the image (large x and y
           positions), then roundoff errors are going to decrease the
           accuracy of the second order calculations.

           When the object and clumps columns were allocated, the
           shift columns were set to NaN so we know when to set them
           for each object/clump. For later parallelization, a mutex
           can be set up within the `if' statement below and another
           check can be done within the mutex, this way, only if
           separate threads get to the start of an object in their
           respective mesh, for an object that doesn't have a shift
           already assigned to it, at the same time (highly unlikely)
           will the process be slightly slowed down. Otherwise, there
           will be no speed penalty. */
        if(isnan(thisobj[ OPOSSHIFTX ]))
          {
            thisobj[ OPOSSHIFTX ] = i%is1+1;
            thisobj[ OPOSSHIFTY ] = i/is1+1;
          }

        /* Set the positional variables: */
        x = i%is1+1;
        y = i/is1+1;
        sx = x - thisobj[ OPOSSHIFTX ];
        sy = y - thisobj[ OPOSSHIFTY ];

        /* Properties that are independent of threshod: */
        ++thisobj[OALLAREA];
        thisobj[ OSKY ]        += p->sky[i];
        thisobj[ OSTD ]        += p->std[i];

        /* Only if the pixel is above the desired threshold.

           REASON: The reason this condition is given like this that
           we don't want to do multiple checks. The basic idea is
           this: when the user doesn't want any thresholds applied,
           then p->threshold=NAN and any conditional that involves a
           NaN will fail, so its logical negation will be positive and
           the calculations below will be done. However, if the user
           does specify a threhold and the pixel is above the
           threshold, then (imgss<p->threshold*p->std[i]) will be
           false and its logigal negation will be positive, so the
           pixel will be included.*/
        if(!(imgss<p->threshold*p->std[i]))
          {
            ++thisobj[OAREA];
            thisobj[ OGeoX ]       += x;
            thisobj[ OGeoY ]       += y;
            thisobj[ OGeoXX ]      += sx * sx;
            thisobj[ OGeoYY ]      += sy * sy;
            thisobj[ OGeoXY ]      += sx * sy;
            thisobj[ OBrightness ] += imgss;
            if(imgss>0)
              {
                thisobj[ OPosBright ]  += imgss;
                thisobj[ OFlxWhtX   ]  += imgss * x;
                thisobj[ OFlxWhtY   ]  += imgss * y;
                thisobj[ OFlxWhtXX  ]  += imgss * sx * sx;
                thisobj[ OFlxWhtYY  ]  += imgss * sy * sy;
                thisobj[ OFlxWhtXY  ]  += imgss * sx * sy;
              }
          }

        if(clumps && clumps[i]>0)
          {
            /* The largest clump ID over each object is the number of
               clumps that object has. */
            thisobj[ ONCLUMPS ] = ( clumps[i] > thisobj[ONCLUMPS]
                                    ? clumps[i] : thisobj[ONCLUMPS] );

            /* Only if the pixel is above the desired threshold, see
               explanation under similar condition above. */
            if(!(imgss<p->threshold*p->std[i]))
              {
                /* Save the information. */
                ++thisobj [ OAREAC ];
                thisobj[ OBrightnessC ]  += imgss;
                thisobj[ OGeoCX  ]       += x;
                thisobj[ OGeoCY  ]       += y;
                thisobj[ OGeoCXX ]       += sx * sx;
                thisobj[ OGeoCYY ]       += sy * sy;
                thisobj[ OGeoCXY ]       += sx * sy;
                if(imgss>0)
                  {
                    thisobj[ OPosBrightC ]  += imgss;
                    thisobj[ OFlxWhtCX   ]  += imgss * x;
                    thisobj[ OFlxWhtCY   ]  += imgss * y;
                    thisobj[ OFlxWhtCXX  ]  += imgss * sx * sx;
                    thisobj[ OFlxWhtCYY  ]  += imgss * sy * sy;
                    thisobj[ OFlxWhtCXY  ]  += imgss * sx * sy;
                  }
              }
          }
      }
}





/* In the second pass, we have the number of clumps so we can find the
   total values for the clumps. */
void
clumppass(struct mkcatalogparams *p)
{
  float imgss;
  double x, y, sx, sy;
  long wngb[2*WNGBSIZE];
  size_t ii, *n, *nf, numngb, ngb[8];
  double *thisclump, *cinfo=p->cinfo;
  long *objects=p->objects, *clumps=p->clumps;
  float *img=p->img, *sky=p->sky, *std=p->std;
  size_t i, j, *ind, *ofcrow, row=0, is0=p->s0, is1=p->s1;



  /* The job of ofcrow (object-first-clump-row) is to give the row
     number of the first clump within an object in the clumps
     information table. This value can be added with the clumpid in
     order to give the row number in the clumps information table. */
  errno=0; ofcrow=malloc((p->numobjects+1)*sizeof *ofcrow);
  if(ofcrow==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for ofcrow in seconpass "
          "(mkcatalog.c)", (p->numobjects+1)*sizeof *ofcrow);


  /* Fill ofcrow using the known number of clumps in each object. Note
     that while ofcrow counts from zero, the clump[i] values (which
     are added with this later) are added with one. So the clump
     information will start from the second row (with index 1) of the
     clumps array, not the first (with index 0). */
  for(i=1;i<=p->numobjects;++i)
    if(p->oinfo[i*OCOLUMNS+ONCLUMPS]>0.0f)
      {
        ofcrow[i]=row;
        row+=p->oinfo[i*OCOLUMNS+ONCLUMPS];
      }

  /* Go over all the pixels in the image and fill in the clump
     information: */
  for(i=0;i<p->s0*p->s1;++i)
    {
      /* We are on a clump, save its properties. */
      if(ISINDEXABLECLPLABEL)
        {
          /* This pointer really simplifies things below! */
          thisclump = ( cinfo
                        + ( ofcrow[objects[i]] + clumps[i] )
                        * CCOLUMNS );

          /* Set the shifts for this object if not already set, see
             the explanations in the firstpass function. */
          if(isnan(thisclump[ CPOSSHIFTX ]))
            {
              thisclump[ CPOSSHIFTX ] = i%is1+1;
              thisclump[ CPOSSHIFTY ] = i/is1+1;
            }

          /* Set the positional variables: */
          x = i%is1+1;
          y = i/is1+1;
          sx = x - thisclump[ CPOSSHIFTX ];
          sy = y - thisclump[ CPOSSHIFTY ];

          /* Calculations that are independent of the threshold. The
             CALLAREA is the full area of the clump irrespective of
             the threshold, this is needed to correctly implement the
             average river flux (which is defined by the. */
          ++thisclump[ CALLAREA ];
          thisclump[ CSKY ]  += sky[i];
          thisclump[ CSTD ]  += std[i];

          /* Only if the pixel is above the desired threshold, see
             explanation under similar condition above. */
          if(!(img[i]-sky[i]<p->threshold*p->std[i]))
            {
              /* Sky subtracted brightness */
              imgss=img[i]-sky[i];

              /* Fill in this clump information. IMPORTANT NOTE: The
                 Sky is not subtracted from the clump brightness or
                 river, because later, we will subtract the river flux
                 from the clump brightness and therefore we don't need
                 to know the Sky for the clump brightness. */
              ++thisclump[ CAREA ];
              thisclump[ CGeoX ]              += x;
              thisclump[ CGeoY ]              += y;
              thisclump[ CGeoXX ]             += sx * sx;
              thisclump[ CGeoYY ]             += sy * sy;
              thisclump[ CGeoXY ]             += sx * sy;
              thisclump[ CINHOSTID ]           = clumps[i];
              thisclump[ CHOSTOID ]            = objects[i];
              thisclump[ CBrightness ]        += img[i];
              thisclump[ CNoRiverBrightness ] += imgss;
              if( imgss > 0.0f )
                {
                  thisclump[ CPosBright ]   += imgss;
                  thisclump[ CFlxWhtX ]     += imgss * x;
                  thisclump[ CFlxWhtY ]     += imgss * y;
                  thisclump[ CFlxWhtXX ]    += imgss * sx * sx;
                  thisclump[ CFlxWhtYY ]    += imgss * sy * sy;
                  thisclump[ CFlxWhtXY ]    += imgss * sx * sy;
                }
            }

        }

      /* We are on a detected region but not a clump (with a negative
         label). This region can be used to find properties like the
         river fluxs in the vicinity of clumps. */
      else if (clumps[i]!=GAL_FITS_LONG_BLANK)

        /* We want to check the river pixels in each detection that
           has a clump. Recall that each detection can host more than
           one object. Since all the detection's pixels are given to
           one object, you can check if we are on an object or not
           with the clumps[i]<0 test.

           The process that keeps the labels of the clumps and their
           host object that have already been used for each river
           pixel in wngb wngb is inherited from the getclumpinfo
           function in NoiseChisel's clump.c.

           There is one big difference. When NoiseChisel was
           identifying the clumps, there were no `object's. The clumps
           were all within one detection. But here, the clumps can be
           separated by a one pixel thick river, but belong to
           different objects. So wngb has to keep two value for each
           neighboring clump: the object it belongs to and the clump
           within that object.
        */
        if(clumps[i]<0 && p->oinfo[objects[i]*OCOLUMNS+ONCLUMPS] > 0.0f)
          {
            /* Make the preparations: */
            ii=0;
            ind=&i;
            GAL_NEIGHBORS_FILL_8_ALLIMG;
            nf=(n=ngb)+numngb;
            memset(wngb, 0, sizeof(wngb));

            /* Go over the neighbors and add the flux of this river
               pixel to a clump's information if it has not been added
               already. */
            do
              if(clumps[*n]>0)
                {
                  /* Go over wngb to see if this river pixel's value
                     has been added to a segment or not. */
                  for(j=0;j<ii;++j)
                    if(wngb[j*2]==objects[*n] && wngb[j*2+1]==clumps[*n])
                      /* It is already present. break out. */
                      break;

                  /* First time we are seeing this clump for this
                     river pixel. */
                  if(j==ii)
                    {
                      /* Note that the object label has to come from
                         the object label on the neighbor, not the
                         object label of the river. This river pixel
                         might be immediately between two clumps on
                         separate objects, so it will read it
                         correctly for one clump and incorrectly for
                         the next. */
                      cinfo[ ( ofcrow[objects[*n]] + clumps[*n] )
                             * CCOLUMNS + CRivAve ] += img[i];
                      ++cinfo[ ( ofcrow[objects[*n]] + clumps[*n] )
                               * CCOLUMNS + CRivArea ];
                      wngb[ii*2]=objects[*n];
                      wngb[ii*2+1]=clumps[*n];
                      ++ii;
                    }
                }
            while(++n<nf);
          }
    }

  /* Clean up: */
  free(ofcrow);
}



















/*********************************************************************/
/*****************           Make output           *******************/
/*********************************************************************/
void
makeoutput(struct mkcatalogparams *p)
{
  double sn, pixarea;
  size_t col, *cols, tmpcol;
  char comment[COMMENTSIZE], tline[100], *target;
  int prec[2]={p->floatprecision, p->accuprecision};
  int space[3]={p->intwidth, p->floatwidth, p->accuwidth};


  /* Calculate the pixel area in arcseconds^2: */
  pixarea=gal_wcs_pixel_area_arcsec2(p->wcs);


  /* First make the objects catalog, then the clumps catalog. */
  for(p->obj0clump1=0;p->obj0clump1<2;++p->obj0clump1)
    {

      /* If no clumps image was provided, then ignore the clumps
         catalog. */
      if(p->obj0clump1==1 && p->clumps==NULL)
        continue;


      /* Do the preparations for this round: */
      p->intcounter=p->accucounter=p->curcol=0;
      p->name      = p->obj0clump1 ? "clump" : "object";
      p->icols     = p->obj0clump1 ? CCOLUMNS : OCOLUMNS;
      p->info      = p->obj0clump1 ? p->cinfo : p->oinfo;
      p->cat       = p->obj0clump1 ? p->clumpcat : p->objcat;
      target       = p->obj0clump1 ? MKCATCLUMP : MKCATOBJECT;
      p->filename  = p->obj0clump1 ? p->ccatname : p->ocatname;
      cols         = p->obj0clump1 ? p->clumpcols : p->objcols;
      p->numcols   = p->obj0clump1 ? p->clumpncols : p->objncols;
      p->num       = p->obj0clump1 ? p->numclumps : p->numobjects;



      /* Allocate the integer and accuracy arrays: */
      errno=0; p->intcols=malloc(p->numcols*sizeof *p->intcols);
      if(p->intcols==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for intcols in makeoutput "
              "(mkcatalog.c)", p->numcols*sizeof *p->intcols);
      errno=0; p->accucols=malloc(p->numcols*sizeof *p->accucols);
      if(p->accucols==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for accucols in makeoutput "
              "(mkcatalog.c)", p->numcols*sizeof *p->accucols);



      /* Write the top of the comments: */
      sprintf(comment, "# %s %s catalog.\n", SPACK_STRING, p->name);
      sprintf(p->line, "# %s started on %s", SPACK_NAME, ctime(&p->rawtime));
      strcat(comment, p->line);


      /* Write the input files: */
      strcat(comment, "#\n# Input files and information:\n"
             "# ----------------------------\n");
      sprintf(p->line, "# Input   %s (hdu: %s)\n", p->up.inputname,
              p->cp.hdu);
      strcat(comment, p->line);
      if(p->up.masknameset)
        {
          sprintf(p->line, "# Mask   %s (hdu: %s)\n", p->up.maskname,
                  p->up.mhdu);
          strcat(comment, p->line);
        }
      sprintf(p->line, "# Objects %s (hdu: %s)\n", p->up.objlabsname,
              p->up.objhdu);
      strcat(comment, p->line);
      if(p->up.clumplabsname)
        {
          sprintf(p->line, "# Clumps  %s (hdu: %s)\n", p->up.clumplabsname,
                  p->up.clumphdu);
          strcat(comment, p->line);
        }
      sprintf(p->line, "# Sky     %s (hdu: %s)\n", p->up.skyname,
              p->up.skyhdu);
      strcat(comment, p->line);
      sprintf(p->line, "# Sky STD %s (hdu: %s)\n", p->up.stdname,
              p->up.stdhdu);
      strcat(comment, p->line);
      if(p->obj0clump1==0 && p->upmask)
        {
          sprintf(p->line, "# Upper limit magnitude mask %s (hdu: %s)\n",
                  p->up.upmaskname, p->up.upmaskhdu);
          strcat(comment, p->line);
        }




      /* If a magnitude is also desired, print the zero point
         magnitude and the 5sigma magnitude: */
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n",
              "Zero point magnitude:", p->zeropoint);
      strcat(comment, p->line);
      sprintf(tline, "Pixel %g sigma surface brightness (magnitude)",
              p->nsigmag);
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n",
              tline, -2.5f*log10(p->nsigmag*p->medstd)+p->zeropoint );
      strcat(comment, p->line);

      sn = p->obj0clump1 ? p->clumpsn : p->detsn;
      if(!isnan(sn))
        {
          sprintf(tline, "%s limiting Signal to noise ratio: ",
                  p->obj0clump1 ? "Clump" : "Detection");
          sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n", tline, sn);
          strcat(comment, p->line);
          if(p->obj0clump1==0)
            strcat(comment, "# (NOTE: limits above are for detections, not "
                   "objects)\n");
        }


      /* If cpscorr was used, report it: */
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n",
              "Counts-per-second correction:", 1/p->cpscorr);
      strcat(comment, p->line);



      /* Report the area of each pixel in stradians: */
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%g\n",
              "Pixel area (arcsec^2)", pixarea);
      strcat(comment, p->line);

      /* if a threshold was used then report it: */
      if(!isnan(p->threshold))
        {
          sprintf(p->line, "# "CATDESCRIPTLENGTH"%g\n", "**IMPORTANT** "
                  "Pixel threshold (multiple of local std)",
                  p->threshold);
          strcat(comment, p->line);
        }

      /* If an upper limit was output then report its parameters: */
      if(p->obj0clump1==0 && p->up.upperlimitmagset)
        {
          sprintf(p->line, "# "CATDESCRIPTLENGTH"%lu\n", "Number of upper "
                  "limit magnitude samples", p->upnum);
          strcat(comment, p->line);
          sprintf(p->line, "# "CATDESCRIPTLENGTH"%lu\n", "Number of threads "
                  "used for upper limit magnitude",
                  p->cp.numthreads);
          strcat(comment, p->line);
          sprintf(p->line, "# "CATDESCRIPTLENGTH"%s\n", "Random number "
                  "generator type for upper limit magnitude",
                  gsl_rng_default->name);
          strcat(comment, p->line);
          sprintf(p->line, "# "CATDESCRIPTLENGTH"%lu\n", "Random number "
                  "generator seed for upper limit magnitude",
                  gsl_rng_default_seed);
          strcat(comment, p->line);
          sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n", "STD multiple for "
                  "upper limit magnitude sigma-clip", p->upsclipmultip);
          strcat(comment, p->line);
          sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n", "STD accuracy "
                  "to stop upper limit magnitude sigma-clip", p->upsclipaccu);
          strcat(comment, p->line);
          sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n", "Multiple of "
                  "sigma for final upper limit magnitude", p->upnsigma);
          strcat(comment, p->line);
        }


      /* Prepare for printing the columns: */
      strcat(comment, "#\n# Columns:\n# --------\n");


      /* Fill the catalog array, in the end set the last elements in
         intcols and accucols to -1, so gal_txtarray_array_to_txt knows
         when to stop. */
      for(p->curcol=0;p->curcol<p->numcols;++p->curcol)
        {
          col=cols[p->curcol];
          switch(col)
            {

            case CATID:
              idcol(p);
              break;

            case CATHOSTOBJID:
              hostobj(p, 1);
              break;

            case CATIDINHOSTOBJ:
              hostobj(p, 0);
              break;

            case CATNUMCLUMPS:
              numclumps(p);
              break;

            case CATAREA:
              area(p, 0, 0);
              break;

            case CATCLUMPSAREA:
              area(p, 1, 0);
              break;


            case CATX:
              tmpcol = p->obj0clump1 ? CFlxWhtX : OFlxWhtX;
              position(p, tmpcol, target, MKCATWHTC, MKCATX);
              break;

            case CATY:
              tmpcol = p->obj0clump1 ? CFlxWhtY : OFlxWhtY;
              position(p, tmpcol, target, MKCATWHTC, MKCATY);
              break;

            case CATGEOX:
              tmpcol = p->obj0clump1 ? CGeoX : OGeoX;
              position(p, tmpcol, target, MKCATGEOC, MKCATX);
              break;

            case CATGEOY:
              tmpcol = p->obj0clump1 ? CGeoY : OGeoY;
              position(p, tmpcol, target, MKCATGEOC, MKCATY);
              break;

            case CATCLUMPSX:
              position(p, OFlxWhtCX, MKCATCINO, MKCATWHTC, MKCATX);
              break;

            case CATCLUMPSY:
              position(p, OFlxWhtCY, MKCATCINO, MKCATWHTC, MKCATY);
              break;

            case CATCLUMPSGEOX:
              position(p, OGeoCX, MKCATCINO, MKCATGEOC, MKCATX);
              break;

            case CATCLUMPSGEOY:
              position(p, OGeoCY, MKCATCINO, MKCATGEOC, MKCATY);
              break;

            case CATRA:
              tmpcol = p->obj0clump1 ? CFlxWhtRA : OFlxWhtRA;
              position(p, tmpcol, target, MKCATWHTC, MKCATRA);
              break;

            case CATDEC:
              tmpcol = p->obj0clump1 ? CFlxWhtDec : OFlxWhtDec;
              position(p, tmpcol, target, MKCATWHTC, MKCATDEC);
              break;

            case CATGEORA:
              tmpcol = p->obj0clump1 ? CGeoRA : OGeoRA;
              position(p, tmpcol, target, MKCATGEOC, MKCATRA);
              break;

            case CATGEODEC:
              tmpcol = p->obj0clump1 ? CGeoDec : OGeoDec;
              position(p, tmpcol, target, MKCATGEOC, MKCATDEC);
              break;

            case CATCLUMPSRA:
              position(p, OFlxWhtCRA, MKCATCINO, MKCATWHTC, MKCATRA);
              break;

            case CATCLUMPSDEC:
              position(p, OFlxWhtCDec, MKCATCINO, MKCATWHTC, MKCATDEC);
              break;

            case CATCLUMPSGEORA:
              position(p, OGeoCRA, MKCATCINO, MKCATGEOC, MKCATRA);
              break;

            case CATCLUMPSGEODEC:
              position(p, OGeoCDec, MKCATCINO, MKCATGEOC, MKCATDEC);
              break;

            case CATBRIGHTNESS:
              tmpcol = p->obj0clump1 ? CBrightness : OBrightness;
              brightnessmag(p, tmpcol, target, MKCATBRIGHT);
              break;

            case CATCLUMPSBRIGHTNESS:
              brightnessmag(p, OBrightnessC, MKCATCINO, MKCATBRIGHT);
              break;

            case CATNORIVERBRIGHTNESS:
              brightnessmag(p, CNoRiverBrightness, target, MKCATBRIGHT);
              break;

            case CATMAGNITUDE:
              tmpcol = p->obj0clump1 ? CBrightness : OBrightness;
              brightnessmag(p, tmpcol, target, MKCATMAG);
              break;

            case CATMAGNITUDEERR:
              sncol(p, 1, target);
              break;

            case CATCLUMPSMAGNITUDE:
              brightnessmag(p, OBrightnessC, MKCATCINO, MKCATMAG);
              break;

            case CATUPPERLIMITMAG:
              upperlimitcol(p);
              break;

            case CATRIVERAVE:
              brightnessmag(p, CRivAve, MKRIVERSSUR, MKCATBRIGHT);
              break;

            case CATRIVERNUM:
              area(p, 0, 1);
              break;

            case CATSN:
              sncol(p, 0, target);
              break;

            case CATSKY:
              skystd(p, p->obj0clump1 ? CSKY : OSKY);
              break;

            case CATSTD:
              skystd(p, p->obj0clump1 ? CSTD : OSTD);
              break;

            case CATSEMIMAJOR: case CATSEMIMINOR: case CATPOSITIONANGLE:
            case CATGEOSEMIMAJOR: case CATGEOSEMIMINOR:
            case CATGEOPOSITIONANGLE:
              secondordermoment(p, col, target);
              break;

            default:
              error(EXIT_FAILURE, 0, "a bug! Please contact us at %s so we "
                    "can fix the problem. The value to cols[%lu] (%lu), is "
                    "not recognized in makeoutput (mkcatalog.c)",
                    PACKAGE_BUGREPORT, p->curcol, cols[p->curcol]);
            }

          sprintf(p->line, "# "CATDESCRIPTLENGTH"[%s]\n",
                  p->description, p->unitp);
          strcat(comment, p->line);
        }
      p->intcols[p->intcounter]=p->accucols[p->accucounter]=-1;



      /* Write the catalog to file: */
      gal_txtarray_array_to_txt(p->cat, p->num, p->numcols, comment,
                                p->intcols, p->accucols, space, prec,
                                'f', p->filename);

      /* Clean up: */
      free(p->intcols);
      free(p->accucols);
    }
}




















/*********************************************************************/
/*****************          Main function          *******************/
/*********************************************************************/
void
mkcatalog(struct mkcatalogparams *p)
{
  /* Run through the data for the first time: */
  firstpass(p);
  if(p->clumps) clumppass(p);

  /* Write the output: */
  makeoutput(p);

  /* Clean up: */
  free(p->oinfo);
  free(p->cinfo);
}
