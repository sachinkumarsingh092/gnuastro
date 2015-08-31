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
#include <float.h>
#include <string.h>
#include <stdlib.h>

#include "timing.h"
#include "neighbors.h"
#include "txtarrayvv.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "columns.h"
#include "mkcatalog.h"




/* Macro to see if the label is indexable (belongs to an object or
   not). See the explanation in src/noisechisel/label.h. */
#if FITSLONGBLANK<0
#define ISINDEXABLEOBJLABEL (objects[i]>0)
#define ISINDEXABLECLPLABEL (clumps[i]>0)
#else
#define ISINDEXABLEOBJLABEL (objects[i] && objects[i]!=FITSLONGBLANK)
#define ISINDEXABLECLPLABEL (clumps[i] && clumps[i]!=FITSLONGBLANK)
#endif




/*********************************************************************/
/*****************     Fill information tables     *******************/
/*********************************************************************/
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
  double *thisobj;
  size_t i, s1=p->s1;
  long *objects=p->objects, *clumps=p->clumps;

  /* Go over the pixels and fill the information array. */
  for(i=0;i<p->s0*p->s1;++i)
    if(ISINDEXABLEOBJLABEL)
      {
        /* thisobj is a pointer to the start of the row in the object
           information array (oinfo). It is mainly used to keep things
           short, simple, less bugy and most importantly elegant. */
        imgss = p->img[i] - p->sky[i];
        thisobj = p->oinfo + objects[i]*OCOLUMNS;

        /* Add to the flux weighted center: */
        ++thisobj[OAREA];
        thisobj[ OGeoX ]       += i%s1;
        thisobj[ OGeoY ]       += i/s1;
        thisobj[ OBrightness ] += imgss;
        thisobj[ OSKY ]        += p->sky[i];
        thisobj[ OSTD ]        += p->std[i];
        if(imgss>0)
          {
            thisobj[ OPosBright ]  += imgss;
            thisobj[ OFlxWhtX ]    += imgss * (i%s1);
            thisobj[ OFlxWhtY ]    += imgss * (i/s1);
          }

        if(clumps[i]>0)
          {
            /* The largest clump ID over each object is the number of
               clumps that object has. */
            thisobj[ ONCLUMPS ] = ( clumps[i] > thisobj[ONCLUMPS]
                                    ? clumps[i] : thisobj[ONCLUMPS] );

            /* Save the information. */
            ++thisobj [ OAREAC ];
            thisobj[ OBrightnessC ] += imgss;
            thisobj[ OGeoCX ]       += i%s1;
            thisobj[ OGeoCY ]       += i/s1;
            if(imgss>0)
              {
                thisobj[ OPosBrightC ]  += imgss;
                thisobj[ OFlxWhtCX ]    += imgss * (i%s1);
                thisobj[ OFlxWhtCY ]    += imgss * (i/s1);
              }
          }
      }

  /* Make all the correctins (for the averages): */
  for(i=1;i<=p->numobjects;++i)
    {
      /* Set the average sky and its STD: */
      thisobj = p->oinfo + i*OCOLUMNS;
      thisobj[ OSKY ] /= thisobj[ OAREA ];
      thisobj[ OSTD ] /= thisobj[ OAREA ];

      /* Set the flux weighted center of the object. The flux weighted
         center is only meaningful when there was positive flux inside
         the detection. */
      if(OPosBright>0.0f)
        {
          thisobj[ OFlxWhtX ] = thisobj[ OFlxWhtX ]/thisobj[ OPosBright ]+1;
          thisobj[ OFlxWhtY ] = thisobj[ OFlxWhtY ]/thisobj[ OPosBright ]+1;
        }
      else
        {
          thisobj[ OFlxWhtX ] = thisobj[ OGeoX ] / thisobj[ OAREA ] + 1;
          thisobj[ OFlxWhtY ] = thisobj[ OGeoY ] / thisobj[ OAREA ] + 1;
        }

      /* Set the over-all clump information: */
      if(OPosBrightC>0.0f)
        {
          thisobj[ OFlxWhtCX ] = thisobj[ OFlxWhtCX ]/thisobj[OPosBrightC]+1;
          thisobj[ OFlxWhtCY ] = thisobj[ OFlxWhtCY ]/thisobj[OPosBrightC]+1;
        }
      else
        {
          thisobj[ OFlxWhtCX ] = thisobj[ OGeoCX ] / thisobj[ OAREAC ] + 1;
          thisobj[ OFlxWhtCY ] = thisobj[ OGeoCY ] / thisobj[ OAREAC ] + 1;
        }
    }
}





/* In the second pass, we have the number of clumps so we can find
   store the total values for the clumps. In this second round we can
   also find second order moments of the objects if we want to. */
void
secondpass(struct mkcatalogparams *p)
{
  float imgss;
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

          /* Fill in this clump information:  */
          ++thisclump[ CAREA ];
          thisclump[ CGeoX ]        += i%is1;
          thisclump[ CGeoY ]        += i/is1;
          thisclump[ CBrightness ]  += img[i];
          thisclump[ CSKY ]         += sky[i];
          thisclump[ CSTD ]         += std[i];
          thisclump[ CINHOSTID ]     = clumps[i];
          thisclump[ CHOSTOID ]      = objects[i];
          if( (imgss=img[i]-sky[i]) > 0 )
            {
              thisclump[ CPosBright ]   += imgss;
              thisclump[ CFlxWhtX ]     += imgss * (i%is1);
              thisclump[ CFlxWhtY ]     += imgss * (i/is1);
            }
        }

      /* We are on a detected region but not a clump (with a negative
         label). This region can be used to find properties like the
         river fluxs in the vicinity of clumps. */
      else if (clumps[i]!=FITSLONGBLANK)

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
            FILL_NGB_8_ALLIMG;
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

  /* Make the proper corrections:

     1. Divide by total flux to get the flux weighted center.
     2. Divide the total river flux by the number of river pixels.

     Note that it might happen that there are no river pixels (when
     grown clumps were used).
  */
  for(i=1;i<=p->numclumps;++i)
    {
      /* Do the initial corrections: */
      thisclump = p->cinfo + i*CCOLUMNS;
      thisclump[ CSKY ] /= thisclump[ CAREA ];
      thisclump[ CSTD ] /= thisclump[ CAREA ];
      if(thisclump[ CRivArea ] > 0.0f)
        thisclump[ CRivAve ] /= thisclump[ CRivArea ];
      else
        thisclump[ CRivAve ] = thisclump[ CSKY ];

      if(thisclump [ CPosBright ]>0.0f)
        {
          thisclump[ CFlxWhtX ] = thisclump[CFlxWhtX]/thisclump[CPosBright]+1;
          thisclump[ CFlxWhtY ] = thisclump[CFlxWhtY]/thisclump[CPosBright]+1;
        }
      else
        {
          thisclump[ CFlxWhtX ] = thisclump[ CGeoX ] / thisclump[ CAREA ] + 1;
          thisclump[ CFlxWhtY ] = thisclump[ CGeoY ] / thisclump[ CAREA ] + 1;
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
  size_t *cols;
  double sn, pixarea;
  char comment[COMMENTSIZE], tline[100];
  int prec[2]={p->floatprecision, p->accuprecision};
  int space[3]={p->intwidth, p->floatwidth, p->accuwidth};


  /* Calculate the pixel area in arcseconds^2: */
  pixarea=pixelareaarcsec2(p->wcs);


  /* First make the objects catalog, then the clumps catalog. */
  for(p->obj0clump1=0;p->obj0clump1<2;++p->obj0clump1)
    {


      /* Do the preparations for this round: */
      p->intcounter=p->accucounter=p->curcol=0;
      p->name     = p->obj0clump1 ? "clump" : "object";
      p->icols    = p->obj0clump1 ? CCOLUMNS : OCOLUMNS;
      p->info     = p->obj0clump1 ? p->cinfo : p->oinfo;
      p->cat      = p->obj0clump1 ? p->clumpcat : p->objcat;
      p->filename = p->obj0clump1 ? p->ccatname : p->ocatname;
      cols        = p->obj0clump1 ? p->clumpcols : p->objcols;
      p->numcols  = p->obj0clump1 ? p->clumpncols : p->objncols;
      p->num      = p->obj0clump1 ? p->numclumps : p->numobjects;



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
      sprintf(p->line, "# Clumps  %s (hdu: %s)\n", p->up.clumplabsname,
              p->up.clumphdu);
      strcat(comment, p->line);
      sprintf(p->line, "# Sky     %s (hdu: %s)\n", p->up.skyname,
              p->up.skyhdu);
      strcat(comment, p->line);
      sprintf(p->line, "# Sky STD %s (hdu: %s)\n", p->up.stdname,
              p->up.stdhdu);
      strcat(comment, p->line);


      /* If a magnitude is also desired, print the zero point
         magnitude and the 5sigma magnitude: */
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n",
              "Zero point magnitude:", p->zeropoint);
      strcat(comment, p->line);
      sprintf(tline, "Pixel %g sigma surface brightness "
              "(magnitude/arcsec^2):", p->nsigmag);
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n",
              tline, -2.5f*log10(p->nsigmag*p->maxstd/pixarea)+p->zeropoint );
      strcat(comment, p->line);

      sn = p->obj0clump1 ? p->clumpsn : p->detsn;
      sprintf(tline, "%s limiting Signal to noise ratio: ",
              p->obj0clump1 ? "Clump" : "Detection");
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n", tline, sn);
      strcat(comment, p->line);
      sprintf(tline, "%s limiting magnitude: ",
              p->obj0clump1 ? "Clump" : "Detection");
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n",
              tline, -2.5f*log10(sn*p->maxstd)+p->zeropoint);
      strcat(comment, p->line);
      if(p->obj0clump1==0)
          strcat(comment, "# (NOTE: limits above are for detections, not "
                 "objects)\n");


      /* If cpscorr was used, report it: */
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%.3f\n",
              "Counts-per-second correction:", 1/p->cpscorr);
      strcat(comment, p->line);



      /* Report the area of each pixel in stradians: */
      sprintf(p->line, "# "CATDESCRIPTLENGTH"%g\n",
              "Pixel area (arcsec^2)", pixarea);
      strcat(comment, p->line);


      /* Prepare for printing the columns: */
      strcat(comment, "#\n# Columns:\n# --------\n");


      /* Fill the catalog array, in the end set the last elements in
         intcols and accucols to -1, so arraytotxt knows when to
         stop. */
      for(p->curcol=0;p->curcol<p->numcols;++p->curcol)
        {
          switch(cols[p->curcol])
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
              position(p, 0, 1, 0);
              break;

            case CATY:
              position(p, 0, 0, 0);
              break;

            case CATCLUMPSX:
              position(p, 0, 1, 1);
              break;

            case CATCLUMPSY:
              position(p, 0, 0, 1);
              break;

            case CATRA:
              position(p, 1, 1, 0);
              break;

            case CATDEC:
              position(p, 1, 0, 0);
              break;

            case CATCLUMPSRA:
              position(p, 1, 1, 1);
              break;

            case CATCLUMPSDEC:
              position(p, 1, 0, 1);
              break;

            case CATBRIGHTNESS:
              brightnessmag(p, 1, 0, 0);
              break;

            case CATCLUMPSBRIGHTNESS:
              brightnessmag(p, 1, 1, 0);
              break;

            case CATMAGNITUDE:
              brightnessmag(p, 0, 0, 0);
              break;

            case CATCLUMPSMAGNITUDE:
              brightnessmag(p, 0, 1, 0);
              break;

            case CATRIVERAVE:
              brightnessmag(p, 1, 0, 1);
              break;

            case CATRIVERNUM:
              area(p, 0, 1);
              break;

            case CATSKY:
              skystd(p, 1);
              break;

            case CATSTD:
              skystd(p, 0);
              break;

            case CATSN:
              sncol(p);
              break;

            default:
              error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we "
                    "can fix the problem. The value to cols[%lu] (%lu), is "
                    "not recognized in makeoutput (mkcatalog.c).",
                    PACKAGE_BUGREPORT, p->curcol, cols[p->curcol]);
            }

          sprintf(p->line, "# "CATDESCRIPTLENGTH"[%s]\n",
                  p->description, p->unitp);
          strcat(comment, p->line);
        }
      p->intcols[p->intcounter]=p->accucols[p->accucounter]=-1;



      /* Write the catalog to file: */
      arraytotxt(p->cat, p->num, p->numcols, comment, p->intcols,
                 p->accucols, space, prec, 'f', p->filename);

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
  /* Allocate two arrays to keep all the basic information about each
     object and clump. Note that there should be one row more than the
     total number of objects or clumps. This is because we want each
     label to be its row number and we don't have any object label of
     zero.*/
  errno=0; p->oinfo=calloc(OCOLUMNS*(p->numobjects+1), sizeof *p->oinfo);
  if(p->oinfo==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->oinfo in mkcatalog "
          "(mkcatalog.c)", OCOLUMNS*(p->numobjects+1)*sizeof *p->oinfo);
  errno=0; p->cinfo=calloc(CCOLUMNS*(p->numclumps+1), sizeof *p->cinfo);
  if(p->cinfo==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->cinfo in mkcatalog "
          "(mkcatalog.c)", CCOLUMNS*(p->numclumps+1)*sizeof *p->cinfo);



  /* Run through the data for the first time: */
  firstpass(p);
  secondpass(p);


  /* If world coordinates are needed, then do the
     transformations. Note that we are passing the pointer to the
     first X axis value so this function can immediately start calling
     wcsp2s in wcslib. */
  if(p->up.raset || p->up.decset)
    {
      xyarraytoradec(p->wcs, p->oinfo+OCOLUMNS+OFlxWhtX,
                     p->oinfo+OCOLUMNS+OFlxWhtRA, p->numobjects,
                     OCOLUMNS);
      xyarraytoradec(p->wcs, p->cinfo+CCOLUMNS+CFlxWhtX,
                     p->cinfo+CCOLUMNS+CFlxWhtRA, p->numclumps,
                     CCOLUMNS);
    }
  if(p->up.clumpsraset || p->up.clumpsdecset)
    xyarraytoradec(p->wcs, p->oinfo+OCOLUMNS+OFlxWhtCX,
                   p->oinfo+OCOLUMNS+OFlxWhtCRA, p->numobjects,
                   OCOLUMNS);


  /* Write the output: */
  makeoutput(p);

  /* Clean up: */
  free(p->oinfo);
  free(p->cinfo);
}
