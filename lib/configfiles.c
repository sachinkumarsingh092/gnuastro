/*********************************************************************
configfiles -- Read configuration files for each program.
This is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include "checkset.h"
#include "configfiles.h"




/* Add the HOME environment variable to the given directory. */
char *
gal_configfiles_add_home_dir(char *dir)
{
  char *home;

  /* Get the home environment variable. */
  home=getenv("HOME");
  if(home==NULL)
    error(EXIT_FAILURE, 0, "The HOME environment variable "
	  "is not defined.");

  /* Concatenate the two strings together: */
  return gal_checkset_malloc_cat(home, dir);
}





void
gal_configfiles_read_name_value(char *line, char *filename, size_t lineno,
                                char **name, char **value)
{
  int notyetfinished=1, inword=0, inquote=0;

  /* Initialize name and value: */
  *name=NULL;
  *value=NULL;

  /* Go through the characters and set the values: */
  do
    switch(*line)
      {
      case ' ': case '\t': case '\v': case '\n': case '\r':
	if(inword) /* Only considered in a word, not in a quote*/
	  {
	    inword=0;
	    *line='\0';
	    if(*value && inquote==0)
	      notyetfinished=0;
	  }
	break;
      case '#':
	notyetfinished=0;
	break;
      case '"':
	if(inword)
	  error_at_line(EXIT_FAILURE, 0, filename, lineno,
			"Quotes have to be surrounded by whitespace "
			"characters (space, tab, new line, etc).");
	if(inquote)
	  {
	    *line='\0';
	    inquote=0;
	    notyetfinished=0;
	  }
	else
	  {
	    if(*name==NULL)
	      error_at_line(EXIT_FAILURE, 0, filename, lineno,
			    "Parameter name should not start with "
			    "double quotes (\").");
	    inquote=1;
	    *value=line+1;
	  }
	break;
      default:
	if(inword==0 && inquote==0)
	  {
	    if(*name==NULL)
	      *name=line;
	    else  /* name is set, now assign *value. */
	      *value=line;
	    inword=1;
	  }
	break;
      }
  while(*(++line)!='\0' && notyetfinished);

  /* In the last line of the file, there is no new line to be
     converted to a '\0' character! So if value has been assigned, we
     are not in a quote and the line has finished, it means the given
     value has also finished. */
  if(*line=='\0' && *value && inquote==0)
    notyetfinished=0;

  /* This was a blank line: */
  if(*name==NULL && *value==NULL)
    return;

  /* Name or value were set but not yet finished. */
  if(notyetfinished)
    error_at_line(EXIT_FAILURE, 0, filename, lineno,
		  "line finished before parameter name and "
		  "value could be read.");
}





FILE *
gal_configfiles_write_local_config_stop(char *indir, char *filename,
                                        char *spack, char *spack_name,
                                        char **outfilename)
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
		  "running:\n\n    mkdir -p %s\n\nto built it and run "
                  "your previous command again.", indir, indir);
	}
      else
	error(EXIT_FAILURE, errno, "%s", indir);
    }
  else
    {
      errno=0;
      if (closedir(dp)==-1)
	error(EXIT_FAILURE, errno, "%s", indir);
    }


  /* Make the local defaults file and put the top information in
     it. */
  *outfilename=gal_checkset_malloc_cat(indir, filename);

  /* Check if the file opening was successful: */
  errno=0;
  fp=fopen(*outfilename, "w");
  if (fp==NULL)
    error(EXIT_FAILURE, errno, "%s", *outfilename);

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
