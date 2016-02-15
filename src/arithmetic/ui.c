/*********************************************************************
ImageArithmetic - Do arithmetic operations on images.
ImageArithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include "commonargs.h"
#include "configfiles.h"
#include "fitsarrayvv.h"
#include "fixedstringmacros.h"

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
readconfig(char *filename, struct imgarithparams *p)
{
  FILE *fp;
  int junkset;
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
      junkset=0;



      /* Inputs: */
      if(strcmp(name, "hdu")==0)
        allocatecopyset(value, &cp->hdu, &cp->hduset);

      else if (strcmp(name, "mask")==0)
        allocatecopyset(value, &up->maskname, &up->masknameset);

      else if (strcmp(name, "mhdu")==0)
        allocatecopyset(value, &up->mhdu, &up->mhduset);
      else if (strcmp(name, "hdu1")==0)
        {
          if(up->hdus[1]) junkset=1;
          allocatecopyset(value, &up->hdus[1], &junkset);
        }
      else if (strcmp(name, "hdu2")==0)
        {
          if(up->hdus[2]) junkset=1;
          allocatecopyset(value, &up->hdus[2], &junkset);
        }
      else if (strcmp(name, "hdu3")==0)
        {
          if(up->hdus[3]) junkset=1;
          allocatecopyset(value, &up->hdus[3], &junkset);
        }
      else if (strcmp(name, "hdu4")==0)
        {
          if(up->hdus[4]) junkset=1;
          allocatecopyset(value, &up->hdus[4], &junkset);
        }
      else if (strcmp(name, "hdu5")==0)
        {
          if(up->hdus[5]) junkset=1;
          allocatecopyset(value, &up->hdus[5], &junkset);
        }
      else if (strcmp(name, "hdu6")==0)
        {
          if(up->hdus[6]) junkset=1;
          allocatecopyset(value, &up->hdus[6], &junkset);
        }
      else if (strcmp(name, "hdu7")==0)
        {
          if(up->hdus[7]) junkset=1;
          allocatecopyset(value, &up->hdus[7], &junkset);
        }
      else if (strcmp(name, "hdu8")==0)
        {
          if(up->hdus[8]) junkset=1;
          allocatecopyset(value, &up->hdus[8], &junkset);
        }
      else if (strcmp(name, "hdu9")==0)
        {
          if(up->hdus[9]) junkset=1;
          allocatecopyset(value, &up->hdus[9], &junkset);
        }
      else if (strcmp(name, "hdu10")==0)
        {
          if(up->hdus[10]) junkset=1;
          allocatecopyset(value, &up->hdus[10], &junkset);
        }
      else if (strcmp(name, "hdu11")==0)
        {
          if(up->hdus[11]) junkset=1;
          allocatecopyset(value, &up->hdus[11], &junkset);
        }
      else if (strcmp(name, "hdu12")==0)
        {
          if(up->hdus[12]) junkset=1;
          allocatecopyset(value, &up->hdus[12], &junkset);
        }
      else if (strcmp(name, "hdu13")==0)
        {
          if(up->hdus[13]) junkset=1;
          allocatecopyset(value, &up->hdus[13], &junkset);
        }
      else if (strcmp(name, "hdu14")==0)
        {
          if(up->hdus[14]) junkset=1;
          allocatecopyset(value, &up->hdus[14], &junkset);
        }
      else if (strcmp(name, "hdu15")==0)
        {
          if(up->hdus[15]) junkset=1;
          allocatecopyset(value, &up->hdus[15], &junkset);
        }
      else if (strcmp(name, "hdu16")==0)
        {
          if(up->hdus[16]) junkset=1;
          allocatecopyset(value, &up->hdus[16], &junkset);
        }
      else if (strcmp(name, "hdu17")==0)
        {
          if(up->hdus[17]) junkset=1;
          allocatecopyset(value, &up->hdus[17], &junkset);
        }
      else if (strcmp(name, "hdu18")==0)
        {
          if(up->hdus[18]) junkset=1;
          allocatecopyset(value, &up->hdus[18], &junkset);
        }
      else if (strcmp(name, "hdu19")==0)
        {
          if(up->hdus[19]) junkset=1;
          allocatecopyset(value, &up->hdus[19], &junkset);
        }



      /* Outputs */
      else if(strcmp(name, "output")==0)
        allocatecopyset(value, &cp->output, &cp->outputset);


      /* Operating modes: */
      /* Read options common to all programs */
      READ_COMMONOPTIONS_FROM_CONF


      else
	error_at_line(EXIT_FAILURE, 0, filename, lineno,
		      "`%s` not recognized.\n", name);
    }

  /* Set the --hdu value in up->hdus[0]  */
  up->hdus[0]=cp->hdu;

  free(line);
  fclose(fp);
}





void
printvalues(FILE *fp, struct imgarithparams *p)
{
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  /* Print all the options that are set. Separate each group with a
     commented line explaining the options in that group. */
  fprintf(fp, "\n# Input image:\n");
  if(cp->hduset)
    PRINTSTINGMAYBEWITHSPACE("hdu", cp->hdu);
  if(up->masknameset)
    PRINTSTINGMAYBEWITHSPACE("mask", up->maskname);
  if(up->mhdu)
    PRINTSTINGMAYBEWITHSPACE("mhdu", up->mhdu);

  if(up->hdus[1])
    PRINTSTINGMAYBEWITHSPACE("hdu1", up->hdus[1]);
  if(up->hdus[2])
    PRINTSTINGMAYBEWITHSPACE("hdu2", up->hdus[2]);
  if(up->hdus[3])
    PRINTSTINGMAYBEWITHSPACE("hdu3", up->hdus[3]);
  if(up->hdus[4])
    PRINTSTINGMAYBEWITHSPACE("hdu4", up->hdus[4]);
  if(up->hdus[5])
    PRINTSTINGMAYBEWITHSPACE("hdu5", up->hdus[5]);
  if(up->hdus[6])
    PRINTSTINGMAYBEWITHSPACE("hdu6", up->hdus[6]);
  if(up->hdus[7])
    PRINTSTINGMAYBEWITHSPACE("hdu7", up->hdus[7]);
  if(up->hdus[8])
    PRINTSTINGMAYBEWITHSPACE("hdu8", up->hdus[8]);
  if(up->hdus[9])
    PRINTSTINGMAYBEWITHSPACE("hdu9", up->hdus[9]);
  if(up->hdus[10])
    PRINTSTINGMAYBEWITHSPACE("hdu10", up->hdus[10]);
  if(up->hdus[1])
    PRINTSTINGMAYBEWITHSPACE("hdu11", up->hdus[11]);
  if(up->hdus[2])
    PRINTSTINGMAYBEWITHSPACE("hdu12", up->hdus[12]);
  if(up->hdus[3])
    PRINTSTINGMAYBEWITHSPACE("hdu13", up->hdus[13]);
  if(up->hdus[4])
    PRINTSTINGMAYBEWITHSPACE("hdu14", up->hdus[14]);
  if(up->hdus[5])
    PRINTSTINGMAYBEWITHSPACE("hdu15", up->hdus[15]);
  if(up->hdus[6])
    PRINTSTINGMAYBEWITHSPACE("hdu16", up->hdus[16]);
  if(up->hdus[7])
    PRINTSTINGMAYBEWITHSPACE("hdu17", up->hdus[17]);
  if(up->hdus[8])
    PRINTSTINGMAYBEWITHSPACE("hdu18", up->hdus[18]);
  if(up->hdus[9])
    PRINTSTINGMAYBEWITHSPACE("hdu19", up->hdus[19]);

  /* Output: */
  fprintf(fp, "\n# Output:\n");
  if(cp->outputset)
    fprintf(fp, CONF_SHOWFMT"%s\n", "output", cp->output);


  /* For the operating mode, first put the macro to print the common
     options, then the (possible options particular to this
     program). */
  fprintf(fp, "\n# Operating mode:\n");
  PRINT_COMMONOPTIONS;
}






/* Note that numthreads will be used automatically based on the
   configure time. Note that those options which are not mandatory
   must not be listed here. */
void
checkifset(struct imgarithparams *p)
{
  struct stll *t;
  int intro=0, counter=0;
  struct uiparams *up=&p->up;
  struct commonparams *cp=&p->cp;

  for(t=p->tokens; t!=NULL; t=t->next)
    {
      if(nameisfits(t->v))
        {
          p->firstname=t->v;
          switch(counter)
            {
            case 0:
              if(cp->hduset==0) REPORT_NOTSET("hdu"); break;
            case 1:
              if(up->hdus[1]==NULL) REPORT_NOTSET("hdu1"); break;
            case 2:
              if(up->hdus[2]==NULL) REPORT_NOTSET("hdu2"); break;
            case 3:
              if(up->hdus[3]==NULL) REPORT_NOTSET("hdu3"); break;
            case 4:
              if(up->hdus[4]==NULL) REPORT_NOTSET("hdu4"); break;
            case 5:
              if(up->hdus[5]==NULL) REPORT_NOTSET("hdu5"); break;
            case 6:
              if(up->hdus[6]==NULL) REPORT_NOTSET("hdu6"); break;
            case 7:
              if(up->hdus[7]==NULL) REPORT_NOTSET("hdu7"); break;
            case 8:
              if(up->hdus[8]==NULL) REPORT_NOTSET("hdu8"); break;
            case 9:
              if(up->hdus[9]==NULL) REPORT_NOTSET("hdu9"); break;
            case 10:
              if(up->hdus[9]==NULL) REPORT_NOTSET("hdu10"); break;
            case 11:
              if(up->hdus[1]==NULL) REPORT_NOTSET("hdu11"); break;
            case 12:
              if(up->hdus[2]==NULL) REPORT_NOTSET("hdu12"); break;
            case 13:
              if(up->hdus[3]==NULL) REPORT_NOTSET("hdu13"); break;
            case 14:
              if(up->hdus[4]==NULL) REPORT_NOTSET("hdu14"); break;
            case 15:
              if(up->hdus[5]==NULL) REPORT_NOTSET("hdu15"); break;
            case 16:
              if(up->hdus[6]==NULL) REPORT_NOTSET("hdu16"); break;
            case 17:
              if(up->hdus[7]==NULL) REPORT_NOTSET("hdu17"); break;
            case 18:
              if(up->hdus[8]==NULL) REPORT_NOTSET("hdu18"); break;
            case 19:
              if(up->hdus[9]==NULL) REPORT_NOTSET("hdu19"); break;

            default:
              error(EXIT_FAILURE, 0, "Only %d FITS HDUs are given as  "
                    "options, but there are more input FITS images. Please "
                    "specify the HDU values for those images with the "
                    "--hduN options (where N stands for the image number).",
                    counter);
            }
          ++counter;
        }
      else
        /* The token is not a file name, also no operators begin with
           `n', so if the token starts with `n', it must be a negative
           number. */
        if(t->v[0]==NEGDASHREPLACE && isdigit(t->v[1])) t->v[0]='-';
    }


  END_OF_NOTSET_REPORT;
}




















/**************************************************************/
/***************       Sanity Check         *******************/
/**************************************************************/

/* The dash of a negative number will cause problems for the users,
   so to work properly we will go over all the options/arguments and
   if any one starts with a dash and is followed by a number, then
   the dash is replaced by NEGDASHREPLACE. */
void
correctnegativedash(int argc, char *argv[])
{
  size_t i;
  for(i=0;i<argc;++i)
    if(argv[i][0]=='-' && isdigit(argv[i][1]))
      argv[i][0]=NEGDASHREPLACE;
}





/* Standard sanity checks. */
void
sanitycheck(struct imgarithparams *p)
{
  char *token;
  struct stll *correctorder=NULL;

  /* An output file name or mask are only necessary when there is a
     FITS image in the arguments. */
  if(p->firstname)
    {
      /* Set the p->up.maskname accordingly: */
      fileorextname(p->firstname, p->cp.hdu, p->up.masknameset,
                    &p->up.maskname, p->up.mhdu, p->up.mhduset, "mask");

      /* Set the names of the output files: */
      if(p->cp.outputset)
        checkremovefile(p->cp.output, p->cp.dontdelete);
      else
        automaticoutput(p->firstname, "_arith.fits", p->cp.removedirinfo,
                        p->cp.dontdelete, &p->cp.output);
    }
  /* Reorder the tokens so the first one that pops out in imgarith.c
     is the first one that was written by the user. */
  while(p->tokens!=NULL)
    {
      pop_from_stll(&p->tokens, &token);
      add_to_stll(&correctorder, token);
    }
  p->tokens=correctorder;
}




















/**************************************************************/
/************         Set the parameters          *************/
/**************************************************************/
void
setparams(int argc, char *argv[], struct imgarithparams *p)
{
  size_t i;
  struct commonparams *cp=&p->cp;

  /* Set the non-zero initial values, the structure was initialized to
     have a zero value for all elements. */
  cp->spack         = SPACK;
  cp->verb          = 1;
  cp->numthreads    = DP_NUMTHREADS;
  cp->removedirinfo = 1;

  p->tokens         = NULL;
  p->firstname      = NULL;
  p->up.maskname    = NULL;
  for(i=0;i<MAXNUMIMAGES;++i) p->up.hdus[i]=NULL;

  /* The hyphen of a negative number can be confused with a dash, so
     we will temporarily replace such hyphens with other
     characters. */
  correctnegativedash(argc, argv);

  /* Read the arguments. */
  errno=0;
  if(argp_parse(&thisargp, argc, argv, 0, 0, p))
    error(EXIT_FAILURE, errno, "Parsing arguments");

  /* Add the user default values and save them if asked. */
  CHECKSETCONFIG;

  /* Check if all the required parameters are set. */
  checkifset(p);

  /* Print the values for each parameter. */
  if(cp->printparams)
    REPORT_PARAMETERS_SET;

  /* Do a sanity check. */
  sanitycheck(p);
}
