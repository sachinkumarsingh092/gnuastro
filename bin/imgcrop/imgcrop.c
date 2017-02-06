/*********************************************************************
ImageCrop - Crop a given size from one or multiple images.
ImageCrop is part of GNU Astronomy Utilities (Gnuastro) package.

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

#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/threads.h>

#include <timing.h>
#include <checkset.h>

#include "main.h"

#include "crop.h"
#include "wcsmode.h"




/* Write the log entry for each crop.

   A maximum length of FILENAME_BUFFER_IN_VERB characters is set for the
   filename to be displayed in stdout in verbose mode. This length is set
   to make the output on the user's terminal reasonable (in one line). So
   when the filename is longer than this, its first set of characters are
   truncated. In the log-file there is no truncation, therefore the log
   file should be used for checking the outputs, not the outputs printed on
   the screen. */
void
imgcrop_verbose_info(struct cropparams *crp)
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
    asprintf(&msg, "...%s %s: %zu input%s.",
             &crp->name[ outnamelen - FILENAME_BUFFER_IN_VERB + 3 ],
             filestatus, crp->numimg, crp->numimg==1 ?  "" :"s");
  else
    asprintf(&msg, "%-*s %s: %zu input%s.", FILENAME_BUFFER_IN_VERB,
             crp->name, filestatus, crp->numimg, crp->numimg==1 ? "" : "s");

  /* Print the results. */
  gal_timing_report(NULL, msg, 2);
  free(msg);
}





/* Print final statistics in verbose mode. */
void
imgcrop_verbose_final(struct imgcropparams *p)
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
               value of 1. So if `array[i]==0', we know that the file was
               removed. */
            for(i=0;i<p->numout;++i)
              {
                if( ((unsigned char *)(tmp->array))[i] )      ++numcrops;
                if( ((unsigned char *)(tmp->array))[i] == 1 ) ++numcfilled;
              }

            break;
          }

      /* Print the basic information. */
      asprintf(&msg, "%zu crops created.", numcrops);
      gal_timing_report(NULL, msg, 1);
      free(msg);

      /* Only if the user wanted to check the center. */
      if(p->checkcenter)
        {
          asprintf(&msg, "%zu filled in the center.",
                   numcfilled);
          gal_timing_report(NULL, msg, 1);
          free(msg);
        }

      /* Only if there were stitched images. */
      if(numstitched)
        {
          asprintf(&msg, "%zu crops used more than one input.",
                  numstitched);
          gal_timing_report(NULL, msg, 1);
          free(msg);
        }
    }

}





void
imgcrop_write_to_log(struct cropparams *crp)
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
          error(EXIT_FAILURE, 0, "a bug! Please contact us at %s to fix the "
                "problem. For some reason `counter' has become %zu in "
                "`imgcrop_write_to_log'", PACKAGE_BUGREPORT, counter);
        }
    }
}





void *
imgmodecrop(void *inparam)
{
  struct cropparams *crp=(struct cropparams *)inparam;
  struct imgcropparams *p=crp->p;

  size_t i;
  int status;
  struct inputimgs *img;

  /* In image mode, we always only have one image. */
  crp->in_ind=0;

  /* The whole catalog is from one image, so you can get the
     information here:*/
  img=&p->imgs[crp->in_ind];
  crp->infits=gal_fits_read_hdu(img->name, p->cp.hdu, 0);

  /* Go over all the outputs that are assigned to this thread: */
  for(i=0;crp->indexs[i]!=GAL_THREADS_NON_THRD_INDEX;++i)
    {
      /* Set all the output parameters: */
      crp->out_ind=crp->indexs[i];
      crp->outfits=NULL;
      crp->numimg=0;
      cropname(crp);

      /* Crop the image. */
      onecrop(crp);

      /* Check the final output: */
      if(crp->numimg)
        {
          /* Check if the center of the crop is filled or not. */
          crp->centerfilled=iscenterfilled(crp);

          /* Add the final headers and close output FITS image: */
          gal_fits_write_keys_version(crp->outfits, NULL, PROGRAM_STRING);
          status=0;
          if( fits_close_file(crp->outfits, &status) )
            gal_fits_io_error(status, "CFITSIO could not close "
                                   "the opened file");

          /* Remove the output image if its center was not filled. */
          if(crp->centerfilled==0)
            {
              errno=0;
              if(unlink(crp->name))
                error(EXIT_FAILURE, errno, "can't delet %s (center"
                      "was blank)", crp->name);
            }

        }
      else crp->centerfilled=0;

      /* Report the status on stdout if verbose mode is requested. */
      if(!p->cp.quiet) imgcrop_verbose_info(crp);
      if(p->cp.log)    imgcrop_write_to_log(crp);
    }

  /* Close the input image. */
  status=0;
  if( fits_close_file(crp->infits, &status) )
    gal_fits_io_error(status, "imgmode.c: imgcroponthreads could "
                      "not close FITS file");

  /* Wait until all other threads finish. */
  if(p->cp.numthreads>1)
    pthread_barrier_wait(crp->b);

  return NULL;
}





void *
wcsmodecrop(void *inparam)
{
  struct cropparams *crp=(struct cropparams *)inparam;
  struct imgcropparams *p=crp->p;

  size_t i;
  int status;


  /* Go over all the output objects for this thread. */
  for(i=0;crp->indexs[i]!=GAL_THREADS_NON_THRD_INDEX;++i)
    {
      /* Set all the output parameters: */
      crp->out_ind=crp->indexs[i];
      crp->outfits=NULL;
      crp->name=NULL;
      crp->numimg=0;


      /* Set the sides of the crop in RA and Dec */
      setcsides(crp);


      /* Go over all the images to see if this target is within their
         range or not. */
      crp->in_ind=0;
      do
        if(radecoverlap(crp))
          {
            /* Open the input FITS file. */
            crp->infits=gal_fits_read_hdu(p->imgs[crp->in_ind].name,
                                          p->cp.hdu, 0);

            /* If a name isn't set yet, set it. */
            if(crp->name==NULL) cropname(crp);

            /* Do the crop. */
            onecrop(crp);

            /* Close the file. */
            status=0;
            if( fits_close_file(crp->infits, &status) )
              gal_fits_io_error(status, "imgmode.c: imgcroponthreads "
                                     "could not close FITS file");
          }
      while ( ++(crp->in_ind) < p->numin );


      /* Check the final output: */
      if(crp->numimg)
        {
          crp->centerfilled=iscenterfilled(crp);

          gal_fits_write_keys_version(crp->outfits, NULL, PROGRAM_STRING);
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
          cropname(crp);
          crp->centerfilled=0;
        }



      /* Report the status on stdout if verbose mode is requested. */
      if(!p->cp.quiet) imgcrop_verbose_info(crp);
      if(p->cp.log)    imgcrop_write_to_log(crp);
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
imgcrop(struct imgcropparams *p)
{
  int err=0;
  char *comments;
  pthread_t t; /* We don't use the thread id, so all are saved here. */
  pthread_attr_t attr;
  pthread_barrier_t b;
  struct cropparams *crp;
  size_t i, *indexs, thrdcols;
  size_t nt=p->cp.numthreads, nb;
  void *(*modefunction)(void *)=NULL;


  /* Set the function to run: */
  modefunction = p->mode==IMGCROP_MODE_IMG ? &imgmodecrop : &wcsmodecrop;


  /* Allocate the array of structures to keep the thread and parameters for
     each thread. */
  errno=0;
  crp=malloc(nt*sizeof *crp);
  if(crp==NULL)
    error(EXIT_FAILURE, errno,
          "%zu bytes in imgcrop (imgcrop.c) for crp", nt*sizeof *crp);


  /* Distribute the indexs into the threads (this is needed even if we
     only have one object where p->cs0 is not defined): */
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
        if(indexs[i*thrdcols]!=GAL_THREADS_NON_THRD_INDEX)
          {
            crp[i].p=p;
            crp[i].b=&b;
            crp[i].indexs=&indexs[i*thrdcols];
            err=pthread_create(&t, &attr, modefunction, &crp[i]);
            if(err)
              error(EXIT_FAILURE, 0, "can't create thread %zu", i);
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
        asprintf(&comments, "# Width of central check box: %zu\n#",
                 p->checkcenter);
      else
        comments=NULL;
      gal_options_print_log(p->log, PROGRAM_STRING, &p->rawtime, comments,
                            LOGFILENAME, &p->cp);
      if(comments) free(comments);
    }

  /* Print the final verbose info, save log, and clean up: */
  imgcrop_verbose_final(p);
  free(indexs);
  free(crp);
}
