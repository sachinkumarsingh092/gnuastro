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
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

#include "timing.h"	/* Includes time.h and sys/time.h   */
#include "checkset.h"
#include "txtarrayvv.h"
#include "commonargs.h"
#include "configfiles.h"
#include "fitsarrayvv.h"

#include "main.h"

#include "ui.h"		        /* Needs main.h                   */
#include "args.h"	        /* Needs main.h, includes argp.h. */


/* Set the file names of the places where the default parameters are
   put. */
#define CONFIG_FILE SPACK CONF_POSTFIX
#define SYSCONFIG_FILE SYSCONFIG_DIR "/" CONFIG_FILE
#define USERCONFIG_FILEEND USERCONFIG_DIR CONFIG_FILE
#define CURDIRCONFIG_FILE CURDIRCONFIG_DIR CONFIG_FILE










/**************************************************************/
/**************       Options and parameters    ***************/
/**************************************************************/
void
readconfig(char *filename, struct mkcatalogparams *p)
{
  int yes;
  FILE *fp;
  size_t lineno=0, len=200;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;
  char key='a';	/* Not used, just a place holder. */

  /* When the file doesn't exist or can't be opened, it is ignored. It
     might be intentional, so there is no error. If a parameter is
     missing, it will be reported after all defaults are read. */
  fp=fopen(filename, "r");
  if (fp==NULL) return;


  /* Allocate some space for `line` with `len` elements so it can
     easily be freed later on. The value of `len` is arbitarary at
     this point, during the run, getline will change it along with the
     pointer to line. */
  errno=0;
  line=malloc(len*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "ui.c: %lu bytes in readdefaults",
	  len * sizeof *line);

  /* Read the tokens in the file:  */
  while(getline(&line, &len, fp) != -1)
    {
      /* Prepare the "name" and "value" strings, also set lineno. */
      STARTREADINGLINE;




      /* Inputs: */
      if(strcmp(name, "hdu")==0)
        allocatecopyset(value, &cp->hdu, &cp->hduset);
      else if (strcmp(name, "mask")==0)
        allocatecopyset(value, &up->maskname, &up->masknameset);
      else if (strcmp(name, "mhdu")==0)
        allocatecopyset(value, &up->mhdu, &up->mhduset);
      else if (strcmp(name, "objlabs")==0)
        allocatecopyset(value, &up->objlabsname, &up->objlabsnameset);
      else if (strcmp(name, "objhdu")==0)
        allocatecopyset(value, &up->objhdu, &up->objhduset);
      else if (strcmp(name, "clumplabs")==0)
        allocatecopyset(value, &up->clumplabsname, &up->clumplabsnameset);
      else if (strcmp(name, "clumphdu")==0)
        allocatecopyset(value, &up->clumphdu, &up->clumphduset);
      else if (strcmp(name, "skyfilename")==0)
        allocatecopyset(value, &up->skyname, &up->skynameset);
      else if (strcmp(name, "skyhdu")==0)
        allocatecopyset(value, &up->skyhdu, &up->skyhduset);
      else if (strcmp(name, "stdfilename")==0)
        allocatecopyset(value, &up->stdname, &up->stdnameset);
      else if (strcmp(name, "stdhdu")==0)
        allocatecopyset(value, &up->stdhdu, &up->stdhduset);
      else if (strcmp(name, "zeropoint")==0)
        {
          if(up->zeropointset) continue;
          anyfloat(value, &p->zeropoint, name, key, SPACK, filename, lineno);
          up->zeropointset=1;
        }
      else if (strcmp(name, "skysubtracted")==0)
        {
          if(up->skysubtractedset) continue;
          intzeroorone(value, &p->skysubtracted, name, key, SPACK,
                       filename, lineno);
          up->skysubtractedset=1;
        }
      else if(strcmp(name, "threshold")==0)
        {
          if(up->thresholdset) continue;
          anydouble(value, &p->threshold, name, key, SPACK,
                   filename, lineno);
          up->thresholdset=1;
        }


      /* Outputs */
      else if(strcmp(name, "output")==0)
        allocatecopyset(value, &cp->output, &cp->outputset);
      else if(strcmp(name, "nsigmag")==0)
        {
          if(up->nsigmagset) continue;
          anydouble(value, &p->nsigmag, name, key, SPACK,
                   filename, lineno);
          up->nsigmagset=1;
        }
      else if(strcmp(name, "intwidth")==0)
        {
          if(up->intwidthset) continue;
          intlzero(value, &p->intwidth, name, key, SPACK, filename, lineno);
          up->intwidthset=1;
        }
      else if(strcmp(name, "floatwidth")==0)
        {
          if(up->floatwidthset) continue;
          intlzero(value, &p->floatwidth, name, key, SPACK, filename, lineno);
          up->floatwidthset=1;
        }
      else if(strcmp(name, "accuwidth")==0)
        {
          if(up->accuwidthset) continue;
          intlzero(value, &p->accuwidth, name, key, SPACK, filename, lineno);
          up->accuwidthset=1;
        }
      else if(strcmp(name, "floatprecision")==0)
        {
          if(up->floatprecisionset) continue;
          intlzero(value, &p->floatprecision, name, key, SPACK,
                   filename, lineno);
          up->floatprecisionset=1;
        }
      else if(strcmp(name, "accuprecision")==0)
        {
          if(up->accuprecisionset) continue;
          intlzero(value, &p->accuprecision, name, key, SPACK,
                   filename, lineno);
          up->accuprecisionset=1;
        }


      /* Catalog columns */
      else if(strcmp(name, "id")==0)
        {
          if(up->idset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATID);
          up->idset=1;
        }
      else if(strcmp(name, "hostobjid")==0)
        {
          if(up->hostobjidset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATHOSTOBJID);
          up->hostobjidset=1;
        }
      else if(strcmp(name, "idinhostobj")==0)
        {
          if(up->idinhostobjset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATIDINHOSTOBJ);
          up->idinhostobjset=1;
        }
      else if(strcmp(name, "numclumps")==0)
        {
          if(up->numclumpsset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATNUMCLUMPS);
          up->numclumpsset=1;
        }
      else if(strcmp(name, "area")==0)
        {
          if(up->areaset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATAREA);
          up->areaset=1;
        }
      else if(strcmp(name, "clumpsarea")==0)
        {
          if(up->clumpsareaset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSAREA);
          up->clumpsareaset=1;
        }
      else if(strcmp(name, "x")==0)
        {
          if(up->xset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATX);
          up->xset=1;
        }
      else if(strcmp(name, "y")==0)
        {
          if(up->yset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATY);
          up->yset=1;
        }
      else if(strcmp(name, "geox")==0)
        {
          if(up->geoxset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATGEOX);
          up->geoxset=1;
        }
      else if(strcmp(name, "geoy")==0)
        {
          if(up->geoyset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATGEOY);
          up->geoyset=1;
        }
      else if(strcmp(name, "clumpsx")==0)
        {
          if(up->clumpsxset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSX);
          up->clumpsxset=1;
        }
      else if(strcmp(name, "clumpsy")==0)
        {
          if(up->clumpsyset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSY);
          up->clumpsyset=1;
        }
      else if(strcmp(name, "clumpsgeox")==0)
        {
          if(up->clumpsgeoxset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSGEOX);
          up->clumpsgeoxset=1;
        }
      else if(strcmp(name, "clumpsgeoy")==0)
        {
          if(up->clumpsgeoyset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSGEOY);
          up->clumpsgeoyset=1;
        }
      else if(strcmp(name, "ra")==0)
        {
          if(up->raset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATRA);
          up->raset=1;
        }
      else if(strcmp(name, "dec")==0)
        {
          if(up->decset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATDEC);
          up->decset=1;
        }
      else if(strcmp(name, "geora")==0)
        {
          if(up->georaset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATGEORA);
          up->georaset=1;
        }
      else if(strcmp(name, "geodec")==0)
        {
          if(up->geodecset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATGEODEC);
          up->geodecset=1;
        }
      else if(strcmp(name, "clumpsra")==0)
        {
          if(up->clumpsraset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSRA);
          up->clumpsraset=1;
        }
      else if(strcmp(name, "clumpsdec")==0)
        {
          if(up->clumpsdecset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSDEC);
          up->clumpsdecset=1;
        }
      else if(strcmp(name, "clumpsgeora")==0)
        {
          if(up->clumpsgeoraset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSGEORA);
          up->clumpsgeoraset=1;
        }
      else if(strcmp(name, "clumpsgeodec")==0)
        {
          if(up->clumpsgeodecset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSGEODEC);
          up->clumpsgeodecset=1;
        }
      else if(strcmp(name, "brightness")==0)
        {
          if(up->brightnessset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATBRIGHTNESS);
          up->brightnessset=1;
        }
      else if(strcmp(name, "clumpsbrightness")==0)
        {
          if(up->clumpsbrightnessset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSBRIGHTNESS);
          p->up.clumpsbrightnessset=1;
        }
      else if(strcmp(name, "noriverbrightness")==0)
        {
          if(up->noriverbrightnessset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATNORIVERBRIGHTNESS);
          p->up.noriverbrightnessset=1;
        }
      else if(strcmp(name, "magnitude")==0)
        {
          if(up->magnitudeset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATMAGNITUDE);
          up->magnitudeset=1;
        }
      else if(strcmp(name, "clumpsmagnitude")==0)
        {
          if(up->clumpsmagnitudeset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATCLUMPSMAGNITUDE);
          up->clumpsmagnitudeset=1;
        }
      else if(strcmp(name, "riverave")==0)
        {
          if(up->riveraveset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATRIVERAVE);
          up->riveraveset=1;
        }
      else if(strcmp(name, "rivernum")==0)
        {
          if(up->rivernumset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATRIVERNUM);
          up->rivernumset=1;
        }
      else if(strcmp(name, "sn")==0)
        {
          if(up->snset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATSN);
          up->snset=1;
        }
      else if(strcmp(name, "sky")==0)
        {
          if(up->skyset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATSKY);
          up->skyset=1;
        }
      else if(strcmp(name, "std")==0)
        {
          if(up->stdset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATSTD);
          up->stdset=1;
        }
      else if(strcmp(name, "semimajor")==0)
        {
          if(up->semimajorset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATSEMIMAJOR);
          up->semimajorset=1;
        }
      else if(strcmp(name, "semiminor")==0)
        {
          if(up->semiminorset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATSEMIMINOR);
          up->semiminorset=1;
        }
      else if(strcmp(name, "positionangle")==0)
        {
          if(up->positionangleset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATPOSITIONANGLE);
          up->positionangleset=1;
        }
      else if(strcmp(name, "geosemimajor")==0)
        {
          if(up->geosemimajorset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATGEOSEMIMAJOR);
          up->geosemimajorset=1;
        }
      else if(strcmp(name, "geosemiminor")==0)
        {
          if(up->geosemiminorset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATGEOSEMIMINOR);
          up->geosemiminorset=1;
        }
      else if(strcmp(name, "geopositionangle")==0)
        {
          if(up->geopositionangleset) continue;
          intzeroorone(value, &yes, name, key, SPACK, filename, lineno);
          if(!yes) continue;
          add_to_sll(&p->allcolsll, CATGEOPOSITIONANGLE);
          up->geopositionangleset=1;
        }



      /* Operating modes */
      /* Read options common to all programs */
      READ_COMMONOPTIONS_FROM_CONF


      else
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      "`%s` not recognized.\n", name);
    }

  free(line);
  fclose(fp);
}





void
printvalues(FILE *fp, struct mkcatalogparams *p)
{
  int i;
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;


  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    PRINTSTINGMAYBEWITHSPACE("hdu", cp->hdu);
  if(up->masknameset)
    PRINTSTINGMAYBEWITHSPACE("mask", up->maskname);
  if(up->mhduset)
    PRINTSTINGMAYBEWITHSPACE("mhdu", up->mhdu);
  if(up->objlabsnameset)
    PRINTSTINGMAYBEWITHSPACE("objlabs", up->objlabsname);
  if(up->objhduset)
    PRINTSTINGMAYBEWITHSPACE("objhdu", up->objhdu);
  if(up->clumplabsnameset)
    PRINTSTINGMAYBEWITHSPACE("clumplabs", up->clumplabsname);
  if(up->clumphduset)
    PRINTSTINGMAYBEWITHSPACE("clumphdu", up->clumphdu);
  if(up->skynameset)
    PRINTSTINGMAYBEWITHSPACE("skyfilename", up->skyname);
  if(up->skyhduset)
    PRINTSTINGMAYBEWITHSPACE("skyhdu", up->skyhdu);
  if(up->stdnameset)
    PRINTSTINGMAYBEWITHSPACE("stdfilename", up->stdname);
  if(up->stdhduset)
    PRINTSTINGMAYBEWITHSPACE("stdhdu", up->stdhdu);
  if(up->zeropointset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "zeropoint", p->zeropoint);
  if(up->skysubtractedset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "skysubtracted", p->skysubtracted);
  if(up->thresholdset)
    fprintf(fp, CONF_SHOWFMT"%.3f\n", "threshold", p->threshold);

  /* Output: */
  fprintf(fp, "\n# Output:\n");
  if(cp->outputset)
    PRINTSTINGMAYBEWITHSPACE("output", cp->output);
  if(up->nsigmagset)
    fprintf(fp, CONF_SHOWFMT"%g\n", "nsigmag", p->nsigmag);
  if(up->intwidthset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "intwidth", p->intwidth);
  if(up->floatwidthset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "floatwidth", p->floatwidth);
  if(up->accuwidthset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "accuwidth", p->accuwidth);
  if(up->floatprecisionset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "floatprecision", p->floatprecision);
  if(up->accuprecisionset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "accuprecision", p->accuprecision);

  /* Catalog columns, since order is important. Notice that they have
     to be printed in opposite order (because of the way they are read
     through a simple linked list). */
  fprintf(fp, "\n# Catalog columns:\n");
  for(i=p->allncols-1;i>=0;--i)
    switch(p->allcols[i])
      {
      case CATID:
        fprintf(fp, CONF_SHOWFMT"%d\n", "id", 1);
        break;
      case CATHOSTOBJID:
        fprintf(fp, CONF_SHOWFMT"%d\n", "hostobjid", 1);
        break;
      case CATIDINHOSTOBJ:
        fprintf(fp, CONF_SHOWFMT"%d\n", "idinhostobj", 1);
        break;
      case CATNUMCLUMPS:
        fprintf(fp, CONF_SHOWFMT"%d\n", "numclumps", 1);
        break;
      case CATAREA:
        fprintf(fp, CONF_SHOWFMT"%d\n", "area", 1);
        break;
      case CATCLUMPSAREA:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsarea", 1);
        break;
      case CATX:
        fprintf(fp, CONF_SHOWFMT"%d\n", "x", 1);
        break;
      case CATY:
        fprintf(fp, CONF_SHOWFMT"%d\n", "y", 1);
        break;
      case CATGEOX:
        fprintf(fp, CONF_SHOWFMT"%d\n", "x", 1);
        break;
      case CATGEOY:
        fprintf(fp, CONF_SHOWFMT"%d\n", "y", 1);
        break;
      case CATCLUMPSX:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsx", 1);
        break;
      case CATCLUMPSY:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsy", 1);
        break;
      case CATCLUMPSGEOX:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsgeox", 1);
        break;
      case CATCLUMPSGEOY:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsgeoy", 1);
        break;
      case CATRA:
        fprintf(fp, CONF_SHOWFMT"%d\n", "ra", 1);
        break;
      case CATDEC:
        fprintf(fp, CONF_SHOWFMT"%d\n", "dec", 1);
        break;
      case CATGEORA:
        fprintf(fp, CONF_SHOWFMT"%d\n", "geora", 1);
        break;
      case CATGEODEC:
        fprintf(fp, CONF_SHOWFMT"%d\n", "geodec", 1);
        break;
      case CATCLUMPSRA:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsra", 1);
        break;
      case CATCLUMPSDEC:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsdec", 1);
        break;
      case CATCLUMPSGEORA:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsgeora", 1);
        break;
      case CATCLUMPSGEODEC:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsgeodec", 1);
        break;
      case CATBRIGHTNESS:
        fprintf(fp, CONF_SHOWFMT"%d\n", "brightness", 1);
        break;
      case CATCLUMPSBRIGHTNESS:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsbrightness", 1);
        break;
      case CATNORIVERBRIGHTNESS:
        fprintf(fp, CONF_SHOWFMT"%d\n", "noriverbrightness", 1);
        break;
      case CATMAGNITUDE:
        fprintf(fp, CONF_SHOWFMT"%d\n", "magnitude", 1);
        break;
      case CATCLUMPSMAGNITUDE:
        fprintf(fp, CONF_SHOWFMT"%d\n", "clumpsmagnitude", 1);
        break;
      case CATRIVERAVE:
        fprintf(fp, CONF_SHOWFMT"%d\n", "riverave", 1);
        break;
      case CATRIVERNUM:
        fprintf(fp, CONF_SHOWFMT"%d\n", "rivernum", 1);
        break;
      case CATSN:
        fprintf(fp, CONF_SHOWFMT"%d\n", "sn", 1);
        break;
      case CATSKY:
        fprintf(fp, CONF_SHOWFMT"%d\n", "sky", 1);
        break;
      case CATSTD:
        fprintf(fp, CONF_SHOWFMT"%d\n", "std", 1);
        break;
      case CATSEMIMAJOR:
        fprintf(fp, CONF_SHOWFMT"%d\n", "semimajor", 1);
        break;
      case CATSEMIMINOR:
        fprintf(fp, CONF_SHOWFMT"%d\n", "semiminor", 1);
        break;
      case CATPOSITIONANGLE:
        fprintf(fp, CONF_SHOWFMT"%d\n", "positionangle", 1);
        break;
      case CATGEOSEMIMAJOR:
        fprintf(fp, CONF_SHOWFMT"%d\n", "geosemimajor", 1);
        break;
      case CATGEOSEMIMINOR:
        fprintf(fp, CONF_SHOWFMT"%d\n", "geosemiminor", 1);
        break;
      case CATGEOPOSITIONANGLE:
        fprintf(fp, CONF_SHOWFMT"%d\n", "geopositionangle", 1);
        break;
      default:
        error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can "
              "fix the problem. For some reason p->allcols[%d] is given "
              "the value %lu which is not recognized in preparearrays "
              "(ui.c).", PACKAGE_BUGREPORT, i, p->allcols[i]);
      }

  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. */
void
checkifset(struct mkcatalogparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  int intro=0;
  if(cp->hduset==0)
    REPORT_NOTSET("hdu");
  if(up->objhduset==0)
    REPORT_NOTSET("objhdu");
  if(up->clumphduset==0)
    REPORT_NOTSET("clumphdu");
  if(up->skyhduset==0)
    REPORT_NOTSET("skyhdu");
  if(up->stdhduset==0)
    REPORT_NOTSET("stdhdu");
  if(up->zeropointset==0)
    REPORT_NOTSET("zeropoint");
  if(up->skysubtractedset==0)
    REPORT_NOTSET("skysubtracted");

  /* Output: */
  if(up->nsigmagset==0)
    REPORT_NOTSET("nsigmag");
  if(up->intwidthset==0)
    REPORT_NOTSET("intwidth");
  if(up->floatwidthset==0)
    REPORT_NOTSET("floatwidth");
  if(up->accuwidthset==0)
    REPORT_NOTSET("accuwidth");
  if(up->floatprecisionset==0)
    REPORT_NOTSET("floatprecision");
  if(up->accuprecisionset==0)
    REPORT_NOTSET("accuprecision");

  END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct mkcatalogparams *p)
{
  struct readheaderkeys keys[2];

  /* Make sure the input file exists. */
  checkfile(p->up.inputname);

  /* Set the names of the files. */
  fileorextname(p->up.inputname, p->cp.hdu, p->up.masknameset,
                &p->up.maskname, p->up.mhdu, p->up.mhduset, "mask");
  fileorextname(p->up.inputname, p->cp.hdu, p->up.objlabsnameset,
                &p->up.objlabsname, p->up.objhdu, p->up.objhduset,
                "object labels");
  fileorextname(p->up.inputname, p->cp.hdu, p->up.clumplabsnameset,
                &p->up.clumplabsname, p->up.clumphdu, p->up.clumphduset,
                "clump labels");
  fileorextname(p->up.inputname, p->cp.hdu, p->up.skynameset,
                &p->up.skyname, p->up.skyhdu, p->up.skyhduset,
                "sky value image");
  fileorextname(p->up.inputname, p->cp.hdu, p->up.stdnameset,
                &p->up.stdname, p->up.stdhdu, p->up.stdhduset,
                "sky standard deviation");

  /* Read the number of labels for the objects:  */
  keys[0].keyname="DETSN";        keys[0].datatype=TDOUBLE;
  keys[1].keyname="NOBJS";        keys[1].datatype=TLONG;
  readkeywords(p->up.objlabsname, p->up.objhdu, keys, 2);
  p->detsn=keys[0].d;
  p->numobjects=keys[1].l;

  /* Read the clumps information. Note that the datatypes don't change. */
  keys[0].keyname="CLUMPSN";
  keys[1].keyname="NCLUMPS";
  readkeywords(p->up.clumplabsname, p->up.clumphdu, keys, 2);
  p->clumpsn=keys[0].d;
  p->numclumps=keys[1].l;

  /* Read the minimum and maximum standard deviation values. */
  keys[0].keyname="MINSTD";       keys[0].datatype=TFLOAT;
  keys[1].keyname="MEDSTD";       keys[1].datatype=TFLOAT;
  readkeywords(p->up.stdname, p->up.stdhdu, keys, 2);
  p->minstd=keys[0].f;
  p->medstd=keys[1].f;
  p->cpscorr = p->minstd>1 ? 1.0f : p->minstd;

  /* When the RA and Dec are needed, make sure that the X and Y
     columns and the RA and Dec columns in the information array are
     immediately after each other and in this order. This is not set
     by the user, but by the programmer. Nevertheless, it is a very
     important issue to check. Since the enum values are constants,
     the compiler will know immediately and stop.

     NOTE: the information array is separate from the output array
  */
  if(p->up.raset || p->up.decset)
    {
      if( OFlxWhtX!=OFlxWhtY-1 || OFlxWhtRA!=OFlxWhtDec-1 )
        error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can "
              "fix the problem. X(%d) and Y(%d), or Ra(%d) and Dec(%d) "
              "columns in the information array are not immediately after "
              "each other.", PACKAGE_BUGREPORT, OFlxWhtX, OFlxWhtY,
              OFlxWhtRA, OFlxWhtDec);
    }

  /* Set the output names: */
  if(p->cp.outputset)
    {
      p->ocatname=malloccat(p->cp.output, "_o.txt");
      p->ccatname=malloccat(p->cp.output, "_c.txt");
    }
  else
    {
      automaticoutput(p->up.inputname, "_o.txt", p->cp.removedirinfo,
                      p->cp.dontdelete, &p->ocatname);
      automaticoutput(p->up.inputname, "_c.txt", p->cp.removedirinfo,
                      p->cp.dontdelete, &p->ccatname);
    }
}



















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
void
checksetlong(struct mkcatalogparams *p, char *filename, char *hdu,
             long **array)
{
  size_t s0, s1;
  int bitpix, anyblank;

  /* Read the file: */
  filetolong(filename, hdu, array, &bitpix, &anyblank, &s0, &s1);

  /* Make sure it has an integer type. */
  if(bitpix==FLOAT_IMG || bitpix==DOUBLE_IMG)
    error(EXIT_FAILURE, 0, "The labels image can be any integer type "
          "(BITPIX). However, %s (hdu: %s) is a %s precision floating "
          "point image.", filename, hdu,
          bitpix==FLOAT_IMG ? "single" : "double");

  /* Make sure it is the same size as the input image. */
  if(s0!=p->s0 || s1!=p->s1)
    error(EXIT_FAILURE, 0, "%s (hdu: %s) is %lu x %lu pixels while the "
          "%s (hdu: %s) is %lu x %lu. The images should have the same "
          "size.", filename, hdu, s1, s0, p->up.inputname,
          p->cp.hdu, p->s1, p->s0);
}





void
checksetfloat(struct mkcatalogparams *p, char *filename, char *hdu,
              float **array)
{
  size_t s0, s1;
  int bitpix, anyblank;

  /* Read the array: */
  filetofloat(filename, NULL, hdu, NULL, array, &bitpix, &anyblank, &s0, &s1);

  /* Make sure it has no blank pixels. */
  if(anyblank)
    error(EXIT_FAILURE, 0, "The Sky and Sky standard deviation images "
          "should not have any blank values. %s (hdu: %s) has blank pixels.",
          filename, hdu);

  /* Make sure it has the same size as the image. */
  if(s0!=p->s0 || s1!=p->s1)
    error(EXIT_FAILURE, 0, "%s (hdu: %s) is %lu x %lu pixels while the "
          "%s (hdu: %s) is %lu x %lu. The images should have the same "
          "size.", filename, hdu, s1, s0, p->up.inputname,
          p->cp.hdu, p->s1, p->s0);
}






void
preparearrays(struct mkcatalogparams *p)
{
  size_t i;
  int bitpix, anyblank;


  /* Prepare the columns and allocate the p->objcols and p->clumpcols
     arrays to keep the macros of what output they should keep. */
  slltoarray(p->allcolsll, &p->allcols, &p->allncols);
  if(p->allncols==0)
    error(EXIT_FAILURE, 0, "No columns specified for output.");
  errno=0; p->objcols=malloc(p->allncols*sizeof *p->objcols);
  if(p->objcols==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->objcols in "
          "preprarearrays (ui.c)", p->allncols*sizeof *p->objcols);
  errno=0; p->clumpcols=malloc(p->allncols*sizeof *p->clumpcols);
  if(p->clumpcols==NULL)
    error(EXIT_FAILURE, errno, "%lu bytes for p->clumpcols in "
          "preprarearrays (ui.c)", p->allncols*sizeof *p->clumpcols);


  /* Put all the input columns into the object or clumps columns and
     count how many columns each should have. */
  p->objncols=p->clumpncols=0;
  for(i=0;i<p->allncols;++i)
    {
      switch(p->allcols[i])
        {
        case CATID:
          p->objcols[p->objncols++]=p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATHOSTOBJID:
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATIDINHOSTOBJ:
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATNUMCLUMPS:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATAREA:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATCLUMPSAREA:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATX:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATY:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATGEOX:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATGEOY:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATCLUMPSX:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATCLUMPSY:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATCLUMPSGEOX:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATCLUMPSGEOY:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATRA:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATDEC:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATGEORA:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATGEODEC:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATCLUMPSRA:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATCLUMPSDEC:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATCLUMPSGEORA:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATCLUMPSGEODEC:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATBRIGHTNESS:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATCLUMPSBRIGHTNESS:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATNORIVERBRIGHTNESS:
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATMAGNITUDE:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATCLUMPSMAGNITUDE:
          p->objcols[p->objncols++] = p->allcols[i];
          break;
        case CATRIVERAVE:
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATRIVERNUM:
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATSN:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATSKY:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATSTD:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATSEMIMAJOR:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATSEMIMINOR:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATPOSITIONANGLE:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATGEOSEMIMAJOR:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATGEOSEMIMINOR:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        case CATGEOPOSITIONANGLE:
          p->objcols[p->objncols++] = p->allcols[i];
          p->clumpcols[p->clumpncols++] = p->allcols[i];
          break;
        default:
          error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can "
                "fix the problem. For some reason p->allcols[%lu] is given "
                "the value %lu which is not recognized in preparearrays "
                "(ui.c).", PACKAGE_BUGREPORT, i, p->allcols[i]);
        }
    }


  /* Read the input image. Note that after this step, everything
     depends on having an input filename. If the user just wants to
     check the parameters, there is no input file name. */
  if(p->up.inputname)
    {
      filetofloat(p->up.inputname, p->up.maskname, p->cp.hdu, p->up.mhdu,
                  &p->img, &bitpix, &anyblank, &p->s0, &p->s1);
      readfitswcs(p->up.inputname, p->cp.hdu, 0, 0, &p->nwcs, &p->wcs);


      /* Read and check the other arrays: */
      checksetlong(p, p->up.objlabsname, p->up.objhdu, &p->objects);
      checksetlong(p, p->up.clumplabsname, p->up.clumphdu, &p->clumps);
      checksetfloat(p, p->up.skyname, p->up.skyhdu, &p->sky);
      checksetfloat(p, p->up.stdname, p->up.stdhdu, &p->std);


      /* Allocate the catalog arrays: */
      if(p->objncols>0 && p->numobjects>0)
        {
          errno=0;
          p->objcat=malloc(p->objncols*p->numobjects*sizeof *p->objcat);
          if(p->objcat==NULL)
            error(EXIT_FAILURE, errno, "%lu bytes for p->objcat in "
                  "preprarearrays (ui.c)",
                  p->objncols*p->numobjects*sizeof *p->objcat);
        }
      else p->objcat=NULL;
      if(p->clumpncols>0 && p->numclumps>0)
        {
          errno=0;
          p->clumpcat=malloc(p->clumpncols*p->numclumps*sizeof *p->clumpcat);
          if(p->clumpcat==NULL)
            error(EXIT_FAILURE, errno, "%lu bytes for p->clumpcat in "
                  "preprarearrays (ui.c)",
                  p->clumpncols*p->numclumps*sizeof *p->clumpcat);
        }
      else p->clumpcat=NULL;
    }


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

  /* Set the shift values to NAN for all the objects and clumps: */
  for(i=1;i<=p->numobjects;++i)
    p->oinfo[i*OCOLUMNS+OPOSSHIFTX]=p->oinfo[i*OCOLUMNS+OPOSSHIFTY]=NAN;
  for(i=1;i<=p->numclumps;++i)
    p->cinfo[i*CCOLUMNS+CPOSSHIFTX]=p->cinfo[i*CCOLUMNS+CPOSSHIFTY]=NAN;

  /* Clean up: */
  freesll(p->allcolsll);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct mkcatalogparams *p)
{
  struct commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = DP_NUMTHREADS;
  cp->removedirinfo = 1;

  p->threshold      = NAN;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "Parsing arguments");


  /* Add the user default values and save them if asked. */
  CHECKSETCONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Do a sanity check. Note that if the user just wants to see the
     parameters and hasn't given any file name, sanity check is
     useless, because in MakeProfiles, sanitycheck just checks the
     file names. So we first have to check if an input */
  if(p->up.inputname)
    sanitycheck(p);


  /* Make the array of input images. */
  preparearrays(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    REPORT_PARAMETERS_SET;

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));
      printf("  - Input   %s (hdu: %s)\n", p->up.inputname, p->cp.hdu);
      if(p->up.masknameset)
        printf("  - Mask   %s (hdu: %s)\n", p->up.maskname, p->up.mhdu);
      printf("  - Objects %s (hdu: %s)\n", p->up.objlabsname,
             p->up.objhdu);
      printf("  - Clumps  %s (hdu: %s)\n", p->up.clumplabsname,
             p->up.clumphdu);
      printf("  - Sky     %s (hdu: %s)\n", p->up.skyname, p->up.skyhdu);
      printf("  - Sky STD %s (hdu: %s)\n", p->up.stdname, p->up.stdhdu);
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct mkcatalogparams *p, struct timeval *t1)
{
  /* Free all the allocated spaces: */
  free(p->sky);
  free(p->std);
  free(p->objcat);
  free(p->cp.hdu);
  free(p->clumps);
  free(p->allcols);
  free(p->objects);
  free(p->objcols);
  free(p->clumpcat);
  free(p->ocatname);
  free(p->ccatname);
  free(p->clumpcols);
  free(p->up.objhdu);
  free(p->cp.output);
  free(p->up.skyhdu);
  free(p->up.stdhdu);
  free(p->up.clumphdu);
  if(p->up.mhduset) free(p->up.mhdu);
  if(p->wcs) wcsvfree(&p->nwcs, &p->wcs);
  if(p->up.skynameset) free(p->up.skyname);
  if(p->up.stdnameset) free(p->up.stdname);
  if(p->up.masknameset) free(p->up.maskname);
  if(p->up.objlabsnameset) free(p->up.objlabsname);
  if(p->up.clumplabsnameset) free(p->up.clumplabsname);

  /* Print the final message. */
  if(p->cp.verb)
    reporttiming(t1, SPACK_NAME" finished in", 0);
}
