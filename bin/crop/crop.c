/*********************************************************************
Crop - Crop a given size from one or multiple images.
Crop is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2015-2020, Free Software Foundation, Inc.

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
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/threads.h>

#include <gnuastro-internal/timing.h>
#include <gnuastro-internal/checkset.h>

#include "main.h"

#include "onecrop.h"
#include "wcsmode.h"




/* Write the log entry for each crop.

   A maximum length of FILENAME_BUFFER_IN_VERB characters is set for the
   filename to be displayed in stdout in verbose mode. This length is set
   to make the output on the user's terminal reasonable (in one line). So
   when the filename is longer than this, its first set of characters are
   truncated. In the log-file there is no truncation, therefore the log
   file should be used for checking the outputs, not the outputs printed on
   the screen. */
static void
crop_verbose_info(struct onecropparams *crp)
{
  char *filestatus, *msg;
  size_t outnamelen=strlen(crp->name);;

  /* Human readable values. */
  filestatus = ( crp->centerfilled==0
                 ? ( crp->numimg == 0
                     ? "no overlap"
                     : "removed (blank center)" )
                 : "created");

  /* Define the output string based on the length of the output file. */
  if ( outnamelen > FILENAME_BUFFER_IN_VERB )
    {
      if( asprintf(&msg, "...%s %s: %zu input%s.",
                   &crp->name[ outnamelen - FILENAME_BUFFER_IN_VERB + 3 ],
                   filestatus, crp->numimg, crp->numimg==1 ?  "" :"s")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
    }
  else
    {
      if( asprintf(&msg, "%-*s %s: %zu input%s.", FILENAME_BUFFER_IN_VERB,
                   crp->name, filestatus, crp->numimg,
                   crp->numimg==1 ? "" : "s")<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
    }

  /* Print the results. */
  gal_timing_report(NULL, msg, 2);
  free(msg);
}





/* Print final statistics in verbose mode. */
static void
crop_verbose_final(struct cropparams *p)
{
  char *msg;
  gal_data_t *tmp;
  size_t i, counter=0, numcrops=0, numstitched=0, numcfilled=0;

  /* This function is only useful in verbose (non-quiet) mode. */
  if(p->cp.quiet) return;

  /* The information is only available if the user asks for a log file. */
  if(p->catname && p->log)
    {
      /* Get the basic counts. */
      for(tmp=p->log; tmp!=NULL; tmp=tmp->next)
        switch(++counter)
          {
          case 2:
            for(i=0;i<p->numout;++i)
              if( ((unsigned short *)(tmp->array))[i] > 1) ++numstitched;
            break;
          case 3:
            /* When the center wasn't checked it has a value of -1, and
               when it was checked and the center was filled, it has a
               value of 1. So if 'array[i]==0', we know that the file was
               removed. */
            for(i=0;i<p->numout;++i)
              {
                if( ((unsigned char *)(tmp->array))[i] )      ++numcrops;
                if( ((unsigned char *)(tmp->array))[i] == 1 ) ++numcfilled;
              }

            break;
          }

      /* Print the basic information. */
      if( asprintf(&msg, "%zu crops created.", numcrops)<0 )
        error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
      gal_timing_report(NULL, msg, 1);
      free(msg);

      /* Only if the user wanted to check the center. */
      if(p->checkcenter)
        {
          if( asprintf(&msg,"%zu pixels filled in the center.",numcfilled)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_timing_report(NULL, msg, 1);
          free(msg);
        }

      /* Only if there were stitched images. */
      if(numstitched)
        {
          if( asprintf(&msg, "%zu crops used more than one input.",
                       numstitched)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_timing_report(NULL, msg, 1);
          free(msg);
        }
    }

}





static void
crop_write_to_log(struct onecropparams *crp)
{
  char **strarr;
  gal_data_t *tmp;
  size_t counter=0;

  for(tmp=crp->p->log; tmp!=NULL; tmp=tmp->next)
    {
      switch(++counter)
        {
        case 1:
          strarr=tmp->array;
          gal_checkset_allocate_copy(crp->name, &strarr[crp->out_ind]);
          break;

        case 2:
          ((unsigned short *)(tmp->array))[crp->out_ind]=crp->numimg;
          break;

        case 3:
          ((unsigned char *)(tmp->array))[crp->out_ind]=crp->centerfilled;
          break;

        default:
          error(EXIT_FAILURE, 0, "%s: a bug! Please contact us at %s to fix "
                "the problem. The value of %zu is not valid for 'counter'",
                __func__, PACKAGE_BUGREPORT, counter);
        }
    }
}





static void *
crop_mode_img(void *inparam)
{
  struct onecropparams *crp=(struct onecropparams *)inparam;
  struct cropparams *p=crp->p;

  size_t i;
  int status;
  struct inputimgs *img;

  /* In image mode, we always only have one image. */
  crp->in_ind=0;

  /* The whole catalog is from one image, so you can get the
     information here:*/
  img=&p->imgs[crp->in_ind];
  crp->infits=gal_fits_hdu_open_format(img->name, p->cp.hdu, 0);

  /* Go over all the outputs that are assigned to this thread: */
  for(i=0; crp->indexs[i]!=GAL_BLANK_SIZE_T; ++i)
    {
      /* Set all the output parameters: */
      crp->out_ind=crp->indexs[i];
      crp->outfits=NULL;
      crp->numimg=1;   /* In Image mode there is only one input image. */
      onecrop_name(crp);

      /* Crop the image. */
      onecrop(crp);

      /* If there was no overlap, then no FITS pointer is created, so
         'numimg' should be set to zero. */
      if(crp->outfits==NULL) crp->numimg=0;

      /* Check the final output: */
      if(crp->numimg)
        {
          /* Check if the center of the crop is filled or not. */
          crp->centerfilled=onecrop_center_filled(crp);

          /* Add the final headers and close output FITS image: */
          gal_fits_key_write_version_in_ptr(NULL, NULL, crp->outfits);
          status=0;
          if( fits_close_file(crp->outfits, &status) )
            gal_fits_io_error(status, "CFITSIO could not close "
                              "the opened file");

          /* Remove the output image if its center was not filled. */
          if(crp->centerfilled==0)
            {
              errno=0;
              if(unlink(crp->name))
                error(EXIT_FAILURE, errno, "can't delete %s (center"
                      "was blank)", crp->name);
            }

        }
      else crp->centerfilled=0;

      /* Report the status on stdout if verbose mode is requested. */
      if(!p->cp.quiet) crop_verbose_info(crp);
      if(p->cp.log)    crop_write_to_log(crp);
    }

  /* Close the input image. */
  status=0;
  if( fits_close_file(crp->infits, &status) )
    gal_fits_io_error(status, "could not close FITS file");

  /* Wait until all other threads finish. */
  if(p->cp.numthreads>1)
    pthread_barrier_wait(crp->b);

  return NULL;
}





static void *
crop_mode_wcs(void *inparam)
{
  struct onecropparams *crp=(struct onecropparams *)inparam;
  struct cropparams *p=crp->p;

  size_t i;
  int status;


  /* Go over all the output objects for this thread. */
  for(i=0; crp->indexs[i]!=GAL_BLANK_SIZE_T; ++i)
    {
      /* Set all the output parameters: */
      crp->out_ind=crp->indexs[i];
      crp->outfits=NULL;
      crp->name=NULL;
      crp->numimg=0;


      /* Set the sides of the crop in RA and Dec */
      wcsmode_crop_corners(crp);


      /* Go over all the images to see if this target is within their
         range or not. */
      crp->in_ind=0;
      do
        if(wcsmode_overlap(crp))
          {
            /* Open the input FITS file. */
            crp->infits=gal_fits_hdu_open_format(p->imgs[crp->in_ind].name,
                                                 p->cp.hdu, 0);

            /* If a name isn't set yet, set it. */
            if(crp->name==NULL) onecrop_name(crp);

            /* Increment the number of images used (necessary for the
               header keywords that are written in 'onecrop'). Then do the
               crop. */
            ++crp->numimg;
            onecrop(crp);

            /* Close the file. */
            status=0;
            if( fits_close_file(crp->infits, &status) )
              gal_fits_io_error(status, "could not close FITS file");
          }
      while ( ++(crp->in_ind) < p->numin );


      /* Correct in_ind. The loop above went until 'in_ind' is one more
         than the index for the last input image (that is how it exited the
         loop). But 'crp->in_ind' is needed later, so correct it here. */
      --crp->in_ind;


      /* Check the final output: */
      if(crp->numimg)
        {
          /* See if the center is filled. */
          crp->centerfilled=onecrop_center_filled(crp);

          /* Write all the dependency versions and close the file. */
          gal_fits_key_write_version_in_ptr(NULL, NULL, crp->outfits);
          status=0;
          if( fits_close_file(crp->outfits, &status) )
            gal_fits_io_error(status, "CFITSIO could not close the "
                                     "opened file");

          if(crp->centerfilled==0)
            {
              errno=0;
              if(unlink(crp->name))
                error(EXIT_FAILURE, errno, "%s", crp->name);
            }
        }
      else
        {
          onecrop_name(crp);
          crp->centerfilled=0;
        }



      /* Report the status on stdout if verbose mode is requested. */
      if(!p->cp.quiet) crop_verbose_info(crp);
      if(p->cp.log)    crop_write_to_log(crp);
    }

  /* Wait until all other threads finish, then return. */
  if(p->cp.numthreads>1)
    pthread_barrier_wait(crp->b);
  return NULL;
}




















/*******************************************************************/
/**************           Output function           ****************/
/*******************************************************************/
/* Main function for the Image Mode. It is assumed that if only one
   crop box from each input image is desired, the first and last
   pixels are already set, irrespective of how the user specified that
   box.  */
void
crop(struct cropparams *p)
{
  int err=0;
  char *tmp;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct onecropparams *crp;
  size_t i, *indexs, thrdcols;
  gal_list_str_t *comments=NULL;
  size_t nt=p->cp.numthreads, nb;
  void *(*modefunction)(void *)=NULL;


  /* Set the function to run: */
  modefunction = p->mode==IMGCROP_MODE_IMG ? &crop_mode_img : &crop_mode_wcs;


  /* Allocate the array of structures to keep the thread and parameters for
     each thread. */
  errno=0;
  crp=malloc(nt*sizeof *crp);
  if(crp==NULL)
    error(EXIT_FAILURE, errno, "%s: allocating %zu bytes for 'crp'",
          __func__, nt*sizeof *crp);


  /* Distribute the indexs into the threads (for clarity, this is needed
     even if we only have one object). */
  gal_threads_dist_in_threads(p->catname ? p->numout : 1, nt,
                              &indexs, &thrdcols);


  /* Run the job, if there is only one thread, don't go through the
     trouble of spinning off a thread! */
  if(nt==1)
    {
      crp[0].p=p;
      crp[0].indexs=indexs;
      modefunction(&crp[0]);
    }
  else
    {
      /* Initialize the attributes. Note that this running thread
         (that spinns off the nt threads) is also a thread, so the
         number the barrier should be one more than the number of
         threads spinned off. */
      if(p->numout<nt) nb=p->numout+1;
      else             nb=nt+1;
      gal_threads_attr_barrier_init(&attr, &b, nb);

      /* Spin off the threads: */
      for(i=0;i<nt;++i)
        if(indexs[i*thrdcols]!=GAL_BLANK_SIZE_T)
          {
            crp[i].p=p;
            crp[i].b=&b;
            crp[i].indexs=&indexs[i*thrdcols];
            err=pthread_create(&t, &attr, modefunction, &crp[i]);
            if(err)
              error(EXIT_FAILURE, 0, "%s: can't create thread %zu",
                    __func__, i);
          }

      /* Wait for all threads to finish and free the spaces. */
      pthread_barrier_wait(&b);
      pthread_attr_destroy(&attr);
      pthread_barrier_destroy(&b);
    }


  /* Print the log file. */
  if(p->cp.log)
    {
      if(p->checkcenter)
        {
          if( asprintf(&tmp, "Width of central check box (in pixels): %zu",
                       p->checkcenter)<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          gal_list_str_add(&comments, tmp, 0);
        }
      gal_checkset_writable_remove(LOGFILENAME, 0, p->cp.dontdelete);
      gal_table_write_log(p->log, PROGRAM_STRING, &p->rawtime, comments,
                          LOGFILENAME, p->cp.quiet);
      gal_list_str_free(comments, 1);
    }

  /* Print the final verbose info, save log, and clean up: */
  crop_verbose_final(p);
  free(indexs);
  free(crp);
}
