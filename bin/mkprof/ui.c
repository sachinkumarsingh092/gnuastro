/*********************************************************************
MakeProfiles - Create mock astronomical profiles.
MakeProfiles is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <time.h>
#include <errno.h>
#include <error.h>
#include <stdio.h>
#include <stdlib.h>

#include <nproc.h>               /* From Gnulib.                   */

#include <gnuastro/box.h>
#include <gnuastro/wcs.h>
#include <gnuastro/fits.h>
#include <gnuastro/array.h>
#include <gnuastro/timing.h>     /* Includes time.h and sys/time.h */
#include <gnuastro/checkset.h>
#include <gnuastro/txtarray.h>
#include <gnuastro/statistics.h>
#include <gnuastro/commonargs.h>
#include <gnuastro/configfiles.h>

#include "main.h"
#include "mkprof.h"

#include "ui.h"                  /* Needs main.h.                  */
#include "args.h"                /* Needs main.h, includes argp.h. */
#include "oneprofile.h"          /* Needs main.h and mkprof.h.     */


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
readconfig(char *filename, struct mkprofparams *p)
{
  FILE *fp;
  char *line, *name, *value;
  struct uiparams *up=&p->up;
  size_t lineno=0, len=200, tmp;
  struct gal_commonparams *cp=&p->cp;
  char key='a';        /* Not used, just a place holder. */

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
      GAL_CONFIGFILES_START_READING_LINE;

      /* Inputs: */
      if(strcmp(name, "hdu")==0)
        gal_checkset_allocate_copy_set(value, &cp->hdu, &cp->hduset);



      /* Outputs: */
      else if(strcmp(name, "output")==0)
        gal_checkset_allocate_copy_set(value, &cp->output, &cp->outputset);

      else if(strcmp(name, "naxis1")==0)
        {
          if(up->naxis1set) continue;
          gal_checkset_sizet_l_zero(value, &tmp, name, key, SPACK,
                                    filename, lineno);
          p->naxes[0]=tmp;
          up->naxis1set=1;
        }
      else if(strcmp(name, "inputascanvas")==0)
        {
          if(up->inputascanvasset) continue;
          gal_checkset_int_zero_or_one(value, &up->inputascanvas, name,
                                       key, SPACK, filename, lineno);
          up->inputascanvasset=1;
        }
      else if(strcmp(name, "naxis2")==0)
        {
          if(up->naxis2set) continue;
          gal_checkset_sizet_l_zero(value, &tmp, name, key, SPACK,
                                    filename, lineno);
          p->naxes[1]=tmp;
          up->naxis2set=1;
        }
      else if(strcmp(name, "oversample")==0)
        {
          if(up->oversampleset) continue;
          gal_checkset_sizet_l_zero(value, &p->oversample, name, key,
                                    SPACK, filename, lineno);
          up->oversampleset=1;
        }
      else if(strcmp(name, "replace")==0)
        {
          if(up->replaceset) continue;
          gal_checkset_int_zero_or_one(value, &p->replace, name,
                                       key, SPACK, filename, lineno);
          up->replaceset=1;
        }
      else if(strcmp(name, "type")==0)
        {
          if(p->up.typeset) continue;
          gal_checkset_known_types(value, &p->up.type, filename, lineno);
          p->up.typeset=1;
        }




      /* Profiles: */
      else if(strcmp(name, "tunitinp")==0)
        {
          if(up->tunitinpset) continue;
          gal_checkset_int_zero_or_one(value, &p->tunitinp, name, key,
                                       SPACK, filename, lineno);
          up->tunitinpset=1;
        }
      else if(strcmp(name, "numrandom")==0)
        {
          if(up->numrandomset) continue;
          gal_checkset_sizet_l_zero(value, &p->numrandom, name, key, SPACK,
                                    filename, lineno);
          up->numrandomset=1;
        }
      else if(strcmp(name, "tolerance")==0)
        {
          if(up->toleranceset) continue;
          gal_checkset_float_l_0(value, &p->tolerance, name, key, SPACK,
                                 filename, lineno);
          up->toleranceset=1;
        }
      else if(strcmp(name, "zeropoint")==0)
        {
          if(up->zeropointset) continue;
          gal_checkset_any_float(value, &p->zeropoint, name, key, SPACK,
                                 filename, lineno);
          up->zeropointset=1;
        }
      else if(strcmp(name, "prepforconv")==0)
        {
          if(up->prepforconvset) continue;
          gal_checkset_int_zero_or_one(value, &p->up.prepforconv, name, key,
                                       SPACK, filename, lineno);
          up->prepforconvset=1;
        }
      else if(strcmp(name, "xshift")==0)
        {
          if(up->xshiftset) continue;
          gal_checkset_sizet_el_zero(value, &tmp, name, key, SPACK,
                                     filename, lineno);
          p->shift[0]=tmp;
          up->xshiftset=1;
        }
      else if(strcmp(name, "yshift")==0)
        {
          if(up->yshiftset) continue;
          gal_checkset_sizet_el_zero(value, &tmp, name, key, SPACK,
                                     filename, lineno);
          p->shift[1]=tmp;
          up->yshiftset=1;
        }
      else if(strcmp(name, "circumwidth")==0)
        {
          if(up->circumwidthset) continue;
          gal_checkset_double_l_value(value, &p->circumwidth, name, key,
                                      SPACK, MINCIRCUMWIDTH, filename,
                                      lineno);
          up->circumwidthset=1;
        }




      /* Catalog: */
      else if(strcmp(name, "xcol")==0)
        {
          if(up->xcolset) continue;
          gal_checkset_sizet_el_zero(value, &p->xcol, name, key, SPACK,
                                     filename, lineno);
          up->xcolset=1;
        }
      else if(strcmp(name, "ycol")==0)
        {
          if(up->ycolset) continue;
          gal_checkset_sizet_el_zero(value, &p->ycol, name, key, SPACK,
                                     filename, lineno);
          up->ycolset=1;
        }
      else if(strcmp(name, "racol")==0)
        {
          if(up->racolset) continue;
          gal_checkset_sizet_el_zero(value, &p->racol, name, key, SPACK,
                                     filename, lineno);
          up->racolset=1;
        }
      else if(strcmp(name, "deccol")==0)
        {
          if(up->deccolset) continue;
          gal_checkset_sizet_el_zero(value, &p->deccol, name, key, SPACK,
                                     filename, lineno);
          up->deccolset=1;
        }
      else if(strcmp(name, "fcol")==0)
        {
          if(up->fcolset) continue;
          gal_checkset_sizet_el_zero(value, &p->fcol, name, key, SPACK,
                                     filename, lineno);
          up->fcolset=1;
        }
      else if(strcmp(name, "rcol")==0)
        {
          if(up->rcolset) continue;
          gal_checkset_sizet_el_zero(value, &p->rcol, name, key, SPACK,
                                     filename, lineno);
          up->rcolset=1;
        }
      else if(strcmp(name, "ncol")==0)
        {
          if(up->ncolset) continue;
          gal_checkset_sizet_el_zero(value, &p->ncol, name, key, SPACK,
                                     filename, lineno);
          up->ncolset=1;
        }
      else if(strcmp(name, "pcol")==0)
        {
          if(up->pcolset) continue;
          gal_checkset_sizet_el_zero(value, &p->pcol, name, key, SPACK,
                                     filename, lineno);
          up->pcolset=1;
        }
      else if(strcmp(name, "qcol")==0)
        {
          if(up->qcolset) continue;
          gal_checkset_sizet_el_zero(value, &p->qcol, name, key, SPACK,
                                     filename, lineno);
          up->qcolset=1;
        }
      else if(strcmp(name, "mcol")==0)
        {
          if(up->mcolset) continue;
          gal_checkset_sizet_el_zero(value, &p->mcol, name, key, SPACK,
                                     filename, lineno);
          up->mcolset=1;
        }
      else if(strcmp(name, "tcol")==0)
        {
          if(up->tcolset) continue;
          gal_checkset_sizet_el_zero(value, &p->tcol, name, key, SPACK,
                                     filename, lineno);
          up->tcolset=1;
        }
      else if(strcmp(name, "mforflatpix")==0)
        {
          if(up->mforflatpixset) continue;
          gal_checkset_int_zero_or_one(value, &p->mforflatpix, name, key,
                                       SPACK, filename, lineno);
          up->mforflatpixset=1;
        }





      /* WCS: */
      else if(strcmp(name, "crpix1")==0)
        {
          if(up->crpix1set) continue;
          gal_checkset_any_double(value, &p->crpix[0], name, key, SPACK,
                                  filename, lineno);
          up->crpix1set=1;
        }
      else if(strcmp(name, "crpix2")==0)
        {
          if(up->crpix2set) continue;
          gal_checkset_any_double(value, &p->crpix[1], name, key, SPACK,
                                  filename, lineno);
          up->crpix2set=1;
        }
      else if(strcmp(name, "crval1")==0)
        {
          if(up->crval1set) continue;
          gal_checkset_any_double(value, &p->crval[0], name, key, SPACK,
                                  filename, lineno);
          up->crval1set=1;
        }
      else if(strcmp(name, "crval2")==0)
        {
          if(up->crval2set) continue;
          gal_checkset_any_double(value, &p->crval[1], name, key, SPACK,
                                  filename, lineno);
          up->crval2set=1;
        }
      else if(strcmp(name, "resolution")==0)
        {
          if(up->resolutionset) continue;
          gal_checkset_any_float(value, &p->resolution, name, key, SPACK,
                                 filename, lineno);
          up->resolutionset=1;
        }



      /* Operating modes: */
      /* Read options common to all programs */
      GAL_CONFIGFILES_READ_COMMONOPTIONS_FROM_CONF


      else
        error_at_line(EXIT_FAILURE, 0, filename, lineno,
                      "`%s` not recognized.\n", name);
    }

  free(line);
  fclose(fp);
}






void
printvalues(FILE *fp, struct mkprofparams *p)
{
  struct uiparams *up=&p->up;
  struct gal_commonparams *cp=&p->cp;

  fprintf(fp, "\n# Input:\n");
  if(cp->hduset)
    GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("hdu", cp->hdu);

  fprintf(fp, "\n# Output:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);
  if(up->naxis1set)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "naxis1", p->naxes[0]);
  if(up->naxis2set)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "naxis2", p->naxes[1]);
  if(up->inputascanvasset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "inputascanvas", up->inputascanvas);
  if(up->oversampleset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "oversample", p->oversample);
  if(up->replaceset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "replace", p->replace);
  if(up->typeset)
    gal_configfiles_print_type(fp, p->up.type);

  fprintf(fp, "\n# Profiles:\n");
  if(up->tunitinpset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "tunitinp", p->tunitinp);
  if(up->numrandomset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "numrandom", p->numrandom);
  if(up->toleranceset)
    fprintf(fp, CONF_SHOWFMT"%.2f\n", "tolerance", p->tolerance);
  if(up->zeropointset)
    fprintf(fp, CONF_SHOWFMT"%.2f\n", "zeropoint", p->zeropoint);
  if(up->circumwidthset)
    fprintf(fp, CONF_SHOWFMT"%.2f\n", "circumwidth", p->circumwidth);

  fprintf(fp, "\n# Catalog:\n");
  if(up->xcolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "xcol", p->xcol);
  if(up->ycolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "ycol", p->ycol);
  if(up->racolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "racol", p->racol);
  if(up->deccolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "deccol", p->deccol);
  if(up->fcolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "fcol", p->fcol);
  if(up->rcolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "rcol", p->rcol);
  if(up->ncolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "ncol", p->ncol);
  if(up->pcolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "pcol", p->pcol);
  if(up->qcolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "qcol", p->qcol);
  if(up->mcolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "mcol", p->mcol);
  if(up->tcolset)
    fprintf(fp, CONF_SHOWFMT"%lu\n", "tcol", p->tcol);
  if(up->mforflatpixset)
    fprintf(fp, CONF_SHOWFMT"%d\n", "mforflatpix", p->mforflatpix);

  fprintf(fp, "\n# WCS:\n");
  if(up->crpix1set)
    fprintf(fp, CONF_SHOWFMT"%g\n", "crpix1", p->crpix[0]);
  if(up->crpix2set)
    fprintf(fp, CONF_SHOWFMT"%g\n", "crpix2", p->crpix[1]);
  if(up->crval1set)
    fprintf(fp, CONF_SHOWFMT"%g\n", "crval1", p->crval[0]);
  if(up->crval2set)
    fprintf(fp, CONF_SHOWFMT"%g\n", "crval2", p->crval[1]);
  if(up->resolutionset)
    fprintf(fp, CONF_SHOWFMT"%g\n", "resolution", p->resolution);

  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating modes:\n");
  GAL_CONFIGFILES_PRINT_COMMONOPTIONS;
}





void
checkifset(struct mkprofparams *p)
{
  struct uiparams *up=&p->up;

  int intro=0;
  if(p->cp.hduset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("hdu");
  if(up->tunitinpset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("tunitinp");
  if(up->numrandomset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("numrandom");
  if(up->toleranceset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("tolerance");
  if(up->zeropointset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("zeropoint");
  if(up->fcolset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("fcol");
  if(up->rcolset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("rcol");
  if(up->ncolset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("ncol");
  if(up->pcolset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("pcol");
  if(up->qcolset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("qcol");
  if(up->mcolset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("mcol");
  if(up->oversampleset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("oversample");
  if(up->circumwidthset==0)
    GAL_CONFIGFILES_REPORT_NOTSET("circumwidth");


  /* The output image size, type, and WCS are only necessary if the user
     has not provided an input image. */
  if(p->up.backname==NULL)
    {
      if(up->typeset==0)
        GAL_CONFIGFILES_REPORT_NOTSET("type");
      if(up->naxis1set==0)
        GAL_CONFIGFILES_REPORT_NOTSET("naxis1");
      if(up->naxis2set==0)
        GAL_CONFIGFILES_REPORT_NOTSET("naxis2");
      if(up->crpix1set==0)
        GAL_CONFIGFILES_REPORT_NOTSET("crpix1");
      if(up->crpix2set==0)
        GAL_CONFIGFILES_REPORT_NOTSET("crpix2");
      if(up->crval1set==0)
        GAL_CONFIGFILES_REPORT_NOTSET("crval1");
      if(up->crval2set==0)
        GAL_CONFIGFILES_REPORT_NOTSET("crval2");
      if(up->resolutionset==0)
        GAL_CONFIGFILES_REPORT_NOTSET("resolution");
    }


  /* The X and Y columns are only needed when the RA and Dec columns have
     not been given. */
  if(up->racolset==0 && up->deccolset==0)
    {
      if(up->xcolset==0)
        GAL_CONFIGFILES_REPORT_NOTSET("xcol");
      if(up->ycolset==0)
        GAL_CONFIGFILES_REPORT_NOTSET("ycol");
    }
  /* The `if' statment above made sure that at least one of the RA and Dec
     columns have been specified. So make sure that both are specified. */
  else if(up->racolset!=up->deccolset)
    {
      if(up->racolset==0)
        GAL_CONFIGFILES_REPORT_NOTSET("racol");
      if(up->deccolset==0)
        GAL_CONFIGFILES_REPORT_NOTSET("deccol");
    }

  GAL_CONFIGFILES_END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/
void
sanitycheck(struct mkprofparams *p)
{
  int d0f1;
  long width[2]={1,1};
  struct uiparams *up=&p->up;
  double truncr, *cat=p->cat, *row;
  size_t i, j, columns[9], cs1=p->cs1;
  char *tmpname=NULL, *xcolstr, *ycolstr;



  /* Check if over-sampling is an odd number, then set/modify the
     respective values.*/
  if(p->oversample%2==0)
    error(EXIT_FAILURE, 0, "the value to the `--oversample' (`-s') option "
          "must be an odd number. Please run the following command for a "
          "complete explanation:\n\n  info gnuastro \"Oversampling\"\n\nOr "
          "See the \"Oversampling\" section of the Gnuastro book.");
  p->halfpixel = 0.5f/p->oversample;
  p->naxes[0] *= p->oversample;
  p->naxes[1] *= p->oversample;



  /* When the RA and Dec columns have been given use them for the profile
     positions instead of the X and Y columns. In the next step we are
     going to convert the RAs and Decs to Xs and Ys and until then, we are
     just dealing with the columns, not the actual values, so it is safe
     (and greatfly simplifies the sanity checks below) to set xcol to ra
     column and ycol to deccol. Also use this check to set the string that
     should be printed if the column is not within the catalog's number of
     columsn*/
  if(up->racolset)
    {
      p->xcol=p->racol;
      p->ycol=p->deccol;
      xcolstr="racol";
      ycolstr="deccol";
    }
  else
    {
      xcolstr="xcol";
      ycolstr="ycol";
    }


  /* If the column numbers are not equal. */
  columns[0]=p->xcol; columns[1]=p->ycol; columns[2]=p->fcol;
  columns[3]=p->rcol; columns[4]=p->ncol; columns[5]=p->pcol;
  columns[6]=p->qcol; columns[7]=p->mcol; columns[8]=p->tcol;
  for(i=0;i<9;++i)
    for(j=0;j<9;++j)
      if(i!=j && columns[i]==columns[j])
        error(EXIT_FAILURE, 0, "at least two of the specified columns "
              "are set to %lu! By adding the `-P` or `--printparams` "
              "option you can check the final column numbers. They "
              "all have to be different", columns[i]);


  /* If all the columns are within the catalog and have proper values. */
  GAL_CHECKSET_CHECK_COL_IN_CAT(p->xcol, xcolstr);
  GAL_CHECKSET_CHECK_COL_IN_CAT(p->ycol, ycolstr);
  GAL_CHECKSET_CHECK_COL_IN_CAT(p->fcol, "fcol");
  GAL_CHECKSET_CHECK_COL_IN_CAT(p->rcol, "rcol");
  GAL_CHECKSET_CHECK_COL_IN_CAT(p->ncol, "ncol");
  GAL_CHECKSET_CHECK_COL_IN_CAT(p->pcol, "pcol");
  GAL_CHECKSET_CHECK_COL_IN_CAT(p->qcol, "qcol");
  GAL_CHECKSET_CHECK_COL_NUM_IN_CAT(p->mcol, "mcol");
  GAL_CHECKSET_CHECK_COL_IN_CAT(p->tcol, "tcol");


  /* If there were terms that gal_txtarray_txt_to_array could not read,
     delete the log file. Note that don't care about the whole input
     catalog, we just want the columns that are important to
     MakeProfiles. The GAL_CHECKSET_CHECK_COL_IN_CAT tests above checked
     those columns and they are fine. For example the user might have some
     alphabetic information in the input file, but as long as the columns
     we want are correct, we have no problem. We can't remove this line
     before those tests, because in their errors, those tests guide the
     reader to check the `txtarray.log' file. */
  gal_checkset_check_remove_file(GAL_TXTARRAY_LOG, 0);


  /* Check if all the profile codes are within the desired range: */
  for(i=0;i<p->cs0;++i)
    if(cat[i*cs1+p->fcol]<0 || cat[i*cs1+p->fcol]>MAXIMUMCODE)
      error(EXIT_FAILURE, 0, "%s: In row %lu, the function code should"
            "be positive and smaller or equal to %d",
            p->up.catname, i+1, MAXIMUMCODE);


  /* If any of xshift or yshift is non-zero, the other should be too!
     Note that conditional operators return 1 if true and 0 if false,
     so if one is non-zero while the other is zero, then sum will be
     1. Otherwise the sum will either be 0 or 2.*/
  switch ( (p->shift[0]!=0) + (p->shift[1]!=0) )
    {
      /* If prepforconv is called, then xshift and yshift should be
         zero. Also a Moffat or Gaussian profile should exist in the
         image. */
    case 0:
      if(p->up.prepforconv)
        {
          /* Check if there is at least one Moffat or Gaussian profile. */
          j=0;
          for(i=0;i<p->cs0;++i)
            if(ispsf(cat[i*cs1+p->fcol]))
              {
                j=i;
                break;
              }

          /* If there is no PSF in the catalog, then you can ignore
             prepforconv. */
          if(i<p->cs0)
            {
              /* Set the row, to simplify: */
              row=&cat[j*cs1];

              /* Find the correct xshift and yshift using the first Moffat
                 or Gaussian profile (in row 'j'). Note that the output of
                 encloseellipse will be the total width, we only want half
                 of it for the shift.*/
              truncr = ( p->tunitinp ? row[p->tcol] :
                         row[p->tcol] * row[p->rcol]/2);
              gal_box_ellipse_in_box(truncr, row[p->qcol]*truncr,
                                     row[p->pcol]*DEGREESTORADIANS, width);
              p->shift[0]  = (width[0]/2)*p->oversample;
              p->shift[1]  = (width[1]/2)*p->oversample;
            }
        }
      break;

    case 1:
      error(EXIT_FAILURE, 0, "at least one of `--xshift` (`-X`) or "
            "`--yshift` (`-Y`) are zero");
      break;

    case 2:
      p->shift[0] *= p->oversample;
      p->shift[1] *= p->oversample;
      break;

    default:
      error(EXIT_FAILURE, 0, "a bug in sanitycheck (ui.c)! In checks "
            "for shifts. Please contact us so we can fix it");
    }
  p->naxes[0] += 2*p->shift[0];
  p->naxes[1] += 2*p->shift[1];


  /* Check the output name: */
  d0f1=gal_checkset_dir_0_file_1(p->cp.output, p->cp.dontdelete);
  if(d0f1)                        /* --output is a file name. */
    {
      p->mergedimgname=p->cp.output;
      p->outdir=gal_checkset_dir_part(p->mergedimgname);
    }
  else                                /* --output is a directory name. */
    {
      errno=0;
      p->outdir=malloc((strlen(p->cp.output)+1)*sizeof *p->outdir);
      if(p->outdir==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for p->outdir in ui.c",
              (strlen(p->cp.output)+1)*sizeof *p->outdir);
      strcpy(p->outdir, p->cp.output);
      gal_checkset_check_dir_write_add_slash(&p->outdir);
      gal_checkset_automatic_output(p->up.catname, ".fits",
                                    p->cp.removedirinfo, p->cp.dontdelete,
                                    &tmpname);
      p->mergedimgname=gal_checkset_malloc_cat(p->outdir, tmpname);
      free(tmpname);
    }
  p->basename=gal_checkset_not_dir_part(p->mergedimgname);
}




















/**************************************************************/
/***************       Preparations         *******************/
/**************************************************************/
void
preparewcs(struct mkprofparams *p)
{
  int status;
  struct wcsprm *wcs;
  long os=p->oversample;

  /* Allocate the memory necessary for the wcsprm structure. */
  errno=0;
  wcs=p->wcs=malloc(sizeof *wcs);
  if(wcs==NULL)
    error(EXIT_FAILURE, errno, "%lu for wcs in preparewcs", sizeof *wcs);

  /* Initialize the structure (allocate all the arrays). */
  wcs->flag=-1;
  if( (status=wcsini(1, 2, wcs)) )
    error(EXIT_FAILURE, 0, "wcsinit error %d: %s",
          status, wcs_errmsg[status]);

  /* Correct the CRPIX values based on oversampling and shifting. */
  p->crpix[0] = p->crpix[0]*os + p->shift[0] - os/2;
  p->crpix[1] = p->crpix[1]*os + p->shift[1] - os/2;

  /* Fill in all the important WCS structure parameters. */
  wcs->equinox=2000.0f;
  wcs->crpix[0]=p->crpix[0];
  wcs->crpix[1]=p->crpix[1];
  wcs->crval[0]=p->crval[0];
  wcs->crval[1]=p->crval[1];
  wcs->pc[0]=-1.0f*p->resolution/3600/p->oversample;
  wcs->pc[3]=p->resolution/3600/p->oversample;
  wcs->pc[1]=wcs->pc[2]=0.0f;
  wcs->cdelt[0]=wcs->cdelt[1]=1.0f;
  strcpy(wcs->cunit[0], "deg");
  strcpy(wcs->cunit[1], "deg");
  strcpy(wcs->ctype[0], "RA---TAN");
  strcpy(wcs->ctype[1], "DEC--TAN");

  /* Set up the wcs structure with the constants defined above. */
  status=wcsset(wcs);
  if(status)
    error(EXIT_FAILURE, 0, "wcsset error %d: %s", status,
          wcs_errmsg[status]);

  /* When individual mode is requested, write the WCS structure to a header
     string to speed up the process: if we don't do it here, this process
     will be necessary on every individual profile's output. */
  if(p->individual)
    {
      status=wcshdo(WCSHDO_safe, wcs, &p->wcsnkeyrec, &p->wcsheader);
      if(status)
        error(EXIT_FAILURE, 0, "wcshdo error %d: %s", status,
              wcs_errmsg[status]);
    }
}





void
preparearrays(struct mkprofparams *p)
{
  void *array=NULL;
  size_t i, naxes[2];
  double *wcstoimg=NULL;
  struct uiparams *up=&p->up;

  /* Allocate space for the log file: */
  p->log=malloc(p->cs0*LOGNUMCOLS*sizeof *p->log);
  if(p->log==NULL)
    error(EXIT_FAILURE, 0, "Allocating %lu bytes for log file",
          p->cs0*LOGNUMCOLS*sizeof *p->log);


  /* If a background image is specified, then use that as the output
     image to build the profiles over. */
  if(up->backname)
    {
      /* Read the input WCS. */
      gal_fits_read_wcs(up->backname, p->cp.hdu, 0, 0, &p->nwcs, &p->wcs);

      /* Read in the background image and its coordinates: */
      p->anyblank=gal_fits_hdu_to_array(up->backname, p->cp.hdu,
                                        &p->bitpix, &array,
                                        &naxes[1], &naxes[0]);
      p->naxes[1]=naxes[1];
      p->naxes[0]=naxes[0];

      /* The the type of the input image is not float, then convert it
         to float to add the mock profile. */
      if(p->bitpix==FLOAT_IMG)
        p->out=array;
      else
        {
          gal_fits_change_type(array, p->bitpix, p->naxes[1]*p->naxes[0],
                               p->anyblank, (void **)(&p->out),
                               FLOAT_IMG);
          free(array);
        }

      /* If the user just wanted the headers, then change all non-NaN
         pixels to 0.0f. */
      if(up->inputascanvas)
        gal_array_freplace_nonnans(p->out, naxes[0]*naxes[1], 0.0f);
    }


  /* Make the WCS structure if it has not been set so far. */
  if(p->wcs==NULL)
    preparewcs(p);

  /* Set the output image type when a background image is not specified,
     or when inputascanvas is called (with a background image). */
  if(!up->backname || up->inputascanvas)
    p->bitpix=up->type;

  /* Convert the RA and Dec to X and Y. We will make a temporary array with
     four rows and the same number of columns, then use that array and the
     WCS structure to get X and Y values. Those X and Y values will then be
     replaced with the RAs and Decs of the original catalog. Recall that
     this is only done in the RAM, not the actual input file, so there is
     no problem. One of the reasons we are doing this is that WCSLIB needs
     the X and Ys and the RA and Decs to be in touching pieces of memory
     and we can't guarantee that (the user might have these columns in
     any order). */
  if(up->racolset)
    {
      /* Allocate the space for the temporary array. */
      errno=0;
      wcstoimg=malloc(4*p->cs0*sizeof *wcstoimg);
      if(wcstoimg==NULL)
        error(EXIT_FAILURE, errno, "%lu bytes for wcstoimg",
              4*p->cs0*sizeof *wcstoimg);


      /* Fill in the first two columns with the RA and Dec. */
      for(i=0;i<p->cs0;++i)
        {
          wcstoimg[ i*4   ] = p->cat[ i*p->cs1+p->racol  ];
          wcstoimg[ i*4+1 ] = p->cat[ i*p->cs1+p->deccol ];
          wcstoimg[ i*4+2 ] = wcstoimg[ i*4+3 ] = NAN;
        }


      /* Convert the X and Y to RA and Dec. */
      gal_wcs_radec_array_to_xy(p->wcs, wcstoimg, wcstoimg+2, p->cs0, 4);


      /* Write the produced X and Y into the input catalog, note that in
         `sanitycheck', xcol became identical to racol and ycol to
         deccol. Note that oversampling has been applied to the WCS
         structure. However, when X and Y are given, oversampling is not
         applied at this point, so we have to correct for the WCS's
         oversampling.*/
      for(i=0;i<p->cs0;++i)
        {
          p->cat[ i*p->cs1+p->xcol  ] = wcstoimg[ i*4+2 ] / p->oversample;
          p->cat[ i*p->cs1+p->ycol  ] = wcstoimg[ i*4+3 ] / p->oversample;
        }


      /* Free the temporary array space. */
      free(wcstoimg);
    }

  /* Allocate the random number generator: */
  gsl_rng_env_setup();
  p->rng=gsl_rng_alloc(gsl_rng_default);
}



















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct mkprofparams *p)
{
  char *jobname;
  struct timeval t1;
  struct gal_commonparams *cp=&p->cp;
  char message[GAL_TIMING_VERB_MSG_LENGTH_V];

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack            = SPACK;
  cp->verb             = 1;
  cp->numthreads       = num_processors(NPROC_CURRENT);
  cp->removedirinfo    = 1;

  p->out               = NULL;
  p->wcs               = NULL;
  p->mforflatpix       = 0;
  p->up.inputascanvas  = 0;

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "parsing arguments");

  /* Add the user default values and save them if asked. */
  GAL_CONFIGFILES_CHECK_SET_CONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    GAL_CONFIGFILES_REPORT_PARAMETERS_SET;

  /* Read catalog if given. */
  gal_txtarray_txt_to_array(p->up.catname, &p->cat, &p->cs0, &p->cs1);

  /* If cp->output was not specified on the command line or in any of
     the configuration files, then automatic output should be used, in
     which case, cp->output should be the current directory. */
  if(p->cp.outputset==0)
    {
      gal_checkset_allocate_copy("./", &p->cp.output);
      p->cp.outputset=1;
    }

  /* Do a sanity check, then remove the possibly existing log file
     created by gal_txtarray_txt_to_array. */
  gettimeofday(&t1, NULL);
  sanitycheck(p);

  /* Prepare the necessary arrays: */
  preparearrays(p);

  /* Everything is ready, notify the user of the program starting. */
  if(cp->verb)
    {
      printf(SPACK_NAME" started on %s", ctime(&p->rawtime));

      errno=0;
      jobname=malloc(strlen(p->up.catname)+100*sizeof *jobname);
      if(jobname==NULL)
        error(EXIT_FAILURE, errno, "jobname in ui.c");
      sprintf(jobname, "%lu profile%sread from %s", p->cs0,
              p->cs0>1?"s ":" ", p->up.catname);
      gal_timing_report(&t1, jobname, 1);
      free(jobname);

      sprintf(message, "Random number generator (RNG) type: %s",
              gsl_rng_name(p->rng));
      gal_timing_report(NULL, message, 1);
      if(p->envseed)
        {
          sprintf(message, "RNG seed for all profiles: %lu",
                  gsl_rng_default_seed);
          gal_timing_report(NULL, message, 1);
        }

      sprintf(message, "Using %lu threads.", cp->numthreads);
      gal_timing_report(NULL, message, 1);
    }
}




















/**************************************************************/
/************      Free allocated, report         *************/
/**************************************************************/
void
freeandreport(struct mkprofparams *p, struct timeval *t1)
{
  int status;

  /* Free all the allocated arrays. */
  free(p->cat);
  free(p->cp.hdu);
  free(p->outdir);
  free(p->basename);
  if(p->individual==0) free(p->log);

  /* p->cp.output might be equal to p->mergedimgname. In this case, if
     we simply free them after each other, there will be a double free
     error. So after freeing output, we set it to NULL since
     free(NULL) is ok.*/
  if(p->cp.output==p->mergedimgname)
    free(p->cp.output);
  else
    {
      free(p->cp.output);
      free(p->mergedimgname);
    }

  /* Free the WCS headers string that was defined for individual mode. */
  if(p->individual)
    free(p->wcsheader);

  /* Free the WCS structure. */
  if( (status=wcsvfree(&p->nwcs, &p->wcs)) )
    error(EXIT_FAILURE, 0, "wcsfree error %d: %s", status,
          wcs_errmsg[status]);

  /* Free the random number generator: */
  gsl_rng_free(p->rng);

  /* Print the final message. */
  if(p->cp.verb)
    gal_timing_report(t1, SPACK_NAME" finished in: ", 0);
}
