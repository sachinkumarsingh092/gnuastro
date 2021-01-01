/*********************************************************************
Gaia Query: retrieve tables from Gaia catalog.
Query is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2020-2021, Free Software Foundation, Inc.

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
#include <string.h>

#include <gnuastro/wcs.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "ui.h"





static void
gaia_sanitycheck(struct queryparams *p)
{
  /* Make sure that atleast one type of constraint is specified. */
  if(p->query==NULL && p->center==NULL && p->overlapwith==NULL)
    error(EXIT_FAILURE, 0, "no '--query', '--center' or '--overlapwith' "
          "specified. At least one of these options are necessary in the "
          "Gaia dataset");

  /* If '--center' is given, '--radius' is also necessary. */
  if(p->center || p->overlapwith)
    {
      /* Make sure the radius is given, and that it isn't zero. */
      if(p->overlapwith==NULL && p->radius==NULL && p->width==NULL)
        error(EXIT_FAILURE, 0, "the '--radius' ('-r') or '--width' ('-w') "
              "options are necessary with the '--center' ('-C') option");

      /* Make sure a dataset is also given. */
      if( p->datasetstr==NULL)
        error(EXIT_FAILURE, 0, "the '--dataset' ('-s') option is necessary "
              "with the '--center' ('-C') option");

      /* Use simpler names for the commonly used datasets. */
      if( !strcmp(p->datasetstr, "edr3") )
        {
          free(p->datasetstr);
          gal_checkset_allocate_copy("gaiaedr3.gaia_source", &p->datasetstr);
        }
      else if( !strcmp(p->datasetstr, "dr2") )
        {
          free(p->datasetstr);
          gal_checkset_allocate_copy("gaiadr2.gaia_source", &p->datasetstr);
        }
      else if( !strcmp(p->datasetstr, "dr1") )
        {
          free(p->datasetstr);
          gal_checkset_allocate_copy("gaiadr1.gaia_source", &p->datasetstr);
        }
      else if( !strcmp(p->datasetstr, "hipparcos") )
        {
          free(p->datasetstr);
          gal_checkset_allocate_copy("public.hipparcos", &p->datasetstr);
        }
      else if( !strcmp(p->datasetstr, "tyco2") )
        {
          free(p->datasetstr);
          gal_checkset_allocate_copy("public.tyco2", &p->datasetstr);
        }
    }
}





void
gaia_query(struct queryparams *p)
{
  size_t ndim;
  gal_data_t *tmp;
  double width2, *center, *darray;
  char *tmpstr, *regionstr, *rangestr=NULL;
  char *command, *columns, allcols[]="*", *querystr;
  double *ocenter=NULL, *owidth=NULL, *omin=NULL, *omax=NULL;

  /* Make sure everything is fine. */
  gaia_sanitycheck(p);


  /* If the raw query has been given, use it. */
  if(p->query)
    querystr=p->query;
  else
    {
      /* If certain columns have been requested use them, otherwise
         download all existing columns.*/
      columns = p->columns ? ui_strlist_to_str(p->columns) : allcols;

      /* If the user wanted an overlap with an image, then calculate it. */
      if(p->overlapwith)
        {
          /* Calculate the Sky coverage of the overlap dataset. */
          gal_wcs_coverage(p->overlapwith, p->cp.hdu, &ndim, &ocenter,
                           &owidth, &omin, &omax);

          /* Make sure a WCS existed in the file. */
          if(owidth==NULL)
            error(EXIT_FAILURE, 0, "%s (hdu %s): contains no WCS to "
                  "derive the sky coverage", p->overlapwith, p->cp.hdu);
        }

      /* For easy reading. */
      center = p->overlapwith ? ocenter : p->center->array;

      /* Write the region. */
      if(p->radius)
        {
          darray=p->radius->array;
          if( asprintf(&regionstr, "CIRCLE('ICRS', %.8f, %.8f, %g)",
                       center[0], center[1], darray[0])<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation ('regionstr')",
                  __func__);
        }
      else if(p->width || p->overlapwith)
        {
          darray = p->overlapwith ? owidth : p->width->array;
          width2 = ( (p->overlapwith || p->width->size==2)
                     ? darray[1] : darray[0] );
          if( asprintf( &regionstr, "BOX('ICRS', %.8f, %.8f, %.8f, %.8f)",
                        center[0], center[1], darray[0], width2 )<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation ('regionstr')",
                  __func__);
        }

      /* Set the range criteria on the requested columns. */
      if(p->range)
        for(tmp=p->range; tmp!=NULL; tmp=tmp->next)
          {
            darray=tmp->array;
            if( asprintf(&tmpstr, "%s%sAND %s>=%g AND %s<=%g",
                         rangestr==NULL ? "" : rangestr,
                         rangestr==NULL ? "" : " ",
                         tmp->name, darray[0], tmp->name, darray[1]) < 0 )
              error(EXIT_FAILURE, 0, "%s: asprintf allocation ('tmpstr')",
                    __func__);
            free(rangestr);
            rangestr=tmpstr;
          }

      /* Write the automatically generated query string. */
      if( asprintf(&querystr,  "SELECT %s "
                   "FROM %s "
                   "WHERE 1=CONTAINS( POINT('ICRS', ra, dec), %s ) %s",
                   columns, p->datasetstr, regionstr,
                   rangestr ? rangestr : "")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation ('querystr')",
              __func__);

      /* Clean up. */
      free(regionstr);
      if(columns!=allcols) free(columns);
      if(p->overlapwith)
        {free(ocenter); free(owidth); free(omin); free(omax);}
    }


  /* Build the calling command. */
  if( asprintf(&command, "curl -o%s --form LANG=ADQL --form FORMAT=fits "
               "--form REQUEST=doQuery --form QUERY=\"%s\" "
               "https://gea.esac.esa.int/tap-server/tap/sync", p->downloadname,
               querystr)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation ('command')", __func__);


  /* Print the calling command for the user to know. */
  if(p->cp.quiet==0)
    printf("Running: %s\n", command);

  /* Run the command. */
  if(system(command))
    error(EXIT_FAILURE, 0, "the query download command %sfailed%s\n",
          p->cp.quiet==0 ? "printed above " : "",
          p->cp.quiet==0 ? "" : " (the command can be printed "
          "if you don't use the option '--quiet', or '-q')");

  /* Clean up. */
  free(command);
}
