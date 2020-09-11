/*********************************************************************
Query - Retreive data from a remote data server.
Query is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2020, Free Software Foundation, Inc.

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
#include <string.h>

#include <gnuastro/pointer.h>

#include <gnuastro-internal/checkset.h>

#include "main.h"
#include "query.h"





char *
query_strlist_to_str(gal_list_str_t *input)
{
  char *out=NULL;
  gal_list_str_t *node;
  size_t n, nn, nnodes=0, alllen=0;

  /* First calculate the full length of all nodes. */
  for(node=input; node!=NULL; node=node->next)
    {
      /* We'll add two extra for each. One for the ',' that must come in
         between it and the next one. One just for a buffer, incase we
         haven't accounted for something. */
      alllen += strlen(node->v) + 2;
      ++nnodes;
    }

  /* Allocate the output string. */
  out=gal_pointer_allocate(GAL_TYPE_STRING, alllen, 1, "out", __func__);

  /* Write all the strings into the allocated space. */
  n=nn=0;
  for(node=input; node!=NULL; node=node->next)
    {
      if(nn++==nnodes-1)
        sprintf(out+n, "%s", node->v);
      else
        n += sprintf(out+n, "%s,", node->v);
    }

  /* Return the merged string. */
  return out;
}





/* Gaia database. */
void
query_gaia_sanitycheck(struct queryparams *p)
{
  /* Make sure that atleast one type of constraint is specified. */
  if(p->center==NULL && p->query==NULL)
    error(EXIT_FAILURE, 0, "no '--center' or '--query' specified. At least "
          "one of these options are necessary in the Gaia dataset");

  /* If '--center' is given, '--radius' is also necessary. */
  if(p->center)
    {
      /* Make sure the radius is given, and that it isn't zero. */
      if( isnan(p->radius) )
        error(EXIT_FAILURE, 0, "the '--radius' ('-r') option is necessary "
              "with the '--center' ('-c') option");

      /* Make sure a dataset is also given. */
      if( p->datasetstr==NULL)
        error(EXIT_FAILURE, 0, "the '--dataset' ('-s') option is necessary "
              "with the '--center' ('-c') option");


      /* Use simpler names for the commonly used datasets. */
      if( !strcmp(p->datasetstr, "dr2") )
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
query_gaia(struct queryparams *p)
{
  double *center;
  char *command, *columns, allcols[]="*", *querystr;

  /* Make sure everything is fine. */
  query_gaia_sanitycheck(p);


  /* If the raw query has been given, use it. */
  if(p->query)
    querystr=p->query;
  else
    {
      /* For easy reading. */
      center=p->center->array;

      /* If certain columns have been requested use them, otherwise
         download all existing columns.

         columns="source_id,ra,dec,phot_g_mean_mag";
      */
      columns = p->columns ? query_strlist_to_str(p->columns) : allcols;

      /* Write the automatically generated query string. */
      if( asprintf(&querystr,  "SELECT %s "
                   "FROM %s "
                   "WHERE 1=CONTAINS( "
                   "POINT('ICRS', ra, dec), "
                   "CIRCLE('ICRS', %.8f, %.8f, %g) )", columns,
                   p->datasetstr, center[0], center[1], p->radius)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation ('querystr')", __func__);

      /* Clean up. */
      if(columns!=allcols) free(columns);
    }


  /* Build the calling command. */
  if( asprintf(&command, "curl -o%s --form LANG=ADQL --form FORMAT=fits "
               "--form REQUEST=doQuery --form QUERY=\"%s\" "
               "https://gea.esac.esa.int/tap-server/tap/sync", p->cp.output,
               querystr)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation ('command')", __func__);


  /* Print the calling command for the user to know. */
  if(p->cp.quiet==0)
    printf("Running: %s\n", command);
  //exit(0);

  /* Run the command. */
  if(system(command))
    error(EXIT_FAILURE, 0, "the query download command %sfailed%s\n",
          p->cp.quiet==0 ? "printed above " : "",
          p->cp.quiet==0 ? "" : " (the command can be printed "
          "if you don't use the option '--quiet', or '-q')");
}


















void
query_check_download(struct queryparams *p)
{
  size_t len;
  int status=0;
  char *logname;
  fitsfile *fptr;

  /* Open the FITS file and if the status value is still zero, it means
     everything worked properly. */
  fits_open_file(&fptr, p->cp.output, READONLY, &status);
  if(status==0) fits_close_file(fptr, &status);
  else
    {
      /* Add a '.log' suffix to the output filename. */
      len=strlen(p->cp.output);
      logname=gal_pointer_allocate(GAL_TYPE_STRING, len+10, 1,
                                   __func__, "logname");
      sprintf(logname, "%s.log", p->cp.output);

      /* Rename the output file to the logname file and let the user
         know. */
      rename(p->cp.output, logname);
      error(EXIT_FAILURE, 0, "the requested dataset could not be retreived! "
            "For more, please see '%s'", logname);
    }
}





void
query(struct queryparams *p)
{
  /* Download the dataset. */
  switch(p->database)
    {
    case QUERY_DATABASE_GAIA: query_gaia(p); break;
    default:
      error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to address "
            "the problem. '%d' is not a recognized database code", __func__,
            PACKAGE_BUGREPORT, p->database);
    }

  /* Make sure that the result is a readable FITS file, otherwise, abort
     with an error. */
  query_check_download(p);

  /* Let the user know that things went well. */
  if(p->cp.quiet==0)
    printf("Query output written to: %s\n", p->cp.output);
}
