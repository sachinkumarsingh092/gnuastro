/*********************************************************************
defaults -- Similar functions for defaults in all utilities.
This is part of GNU Astronomy Utilities (AstrUtils) package.

Copyright (C) 2013-2014 Mohammad Akhlaghi
Tohoku University Astronomical Institute, Sendai, Japan.
http://astr.tohoku.ac.jp/~akhlaghi/

AstrUtils is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

AstrUtils is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with AstrUtils. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#include <time.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>		/* For mkdir permission flags. */

#include "config.h"
#include "defaults.h"
#include "checkset.h"




/* Add the HOME environment variable to the given directory. */
char *
addhomedir(char *dir)
{
  char *home;

  /* Get the home environment variable. */
  home=getenv("HOME");
  if(home==NULL)
    error(EXIT_FAILURE, 0, "The HOME environment variable "
	  "is not defined.");

  /* Concatenate the two strings together: */
  return malloccat(home, dir);
}




FILE *
writelocaldefaultstop(char *indir, char *filename, char *spack,
		      char *spack_name, char **outfilename)
{
  DIR *dp;
  FILE *fp;
  time_t rawtime;

  errno=0;
  time(&rawtime);

  /* Make sure the directory exists, if it doesn't, try to make it.*/
  dp=opendir(indir);
  if(dp==NULL)			/* The directory could not be opened. */
    {
      if(errno==ENOENT)
	{
	  errno=0;
	  if(mkdir(indir, S_IRWXU)==-1)
	    error(EXIT_FAILURE, errno, "%s: Could not be created. Try "
		  "running `mkdir -p %s` to built it and run your "
		  "previous command again.", indir, indir);
	}
      else
	error(EXIT_FAILURE, errno, indir);
    }
  else
    {
      errno=0;
      if (closedir(dp)==-1)
	error(EXIT_FAILURE, errno, indir);
    }


  /* Make the local defaults file and put the top information in
     it. */
  *outfilename=malloccat(indir, filename);

  /* Check if the file opening was successful: */
  errno=0;
  fp=fopen(*outfilename, "w");
  if (fp==NULL)
    error(EXIT_FAILURE, errno, *outfilename);

  /* write the comments: */
  fprintf(fp, "# Default parameters for %s (%s).\n"
	  "# %s is part of GNU Astronomy Utitlies.\n"
	  "# This file was created on %s#\n"
	  "# Use the long option name of each paramter followed by\n"
	  "# a value. The name and value should be separated by\n"
	  "# atleast one of the following charaacters:\n"
	  "# space, `,`, `=` or `:`.\n#\n"
	  "# Run `%s --help` or `info %s`\n"
	  "# for more information.\n#\n"
	  "# NOTE I:  All counting is from zero, not one.\n"
	  "# NOTE II: Lines starting with `#` are ignored.\n",
	  spack_name, spack, spack_name, ctime(&rawtime), spack, spack);

  return fp;
}
