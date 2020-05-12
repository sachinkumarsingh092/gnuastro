/*********************************************************************
BuildProgram: Compile and run programs using Gnuastro's library
BuildProgram is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2020, Free Software Foundation, Inc.

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

#include <gnuastro/list.h>
#include <gnuastro/pointer.h>

#include <main.h>


/* Write the given list into  */
char *
buildprog_as_one_string(char *opt, gal_list_str_t *list)
{
  char *out;
  size_t len=0;
  gal_list_str_t *tmp;

  /* Only if we have a list. */
  if(list)
    {
      /* For every node in the list, we want the 'opt' and a space along with
         the actual string. */
      for(tmp=list; tmp!=NULL; tmp=tmp->next)
        len += 1 + (opt ? strlen(opt) : 0) + strlen(tmp->v);

      /* Allocate space for the string. */
      out=gal_pointer_allocate(GAL_TYPE_UINT8, len+1, 0, __func__, "out");

      /* Write all the strings into the allocated space. */
      len=0;
      for(tmp=list; tmp!=NULL; tmp=tmp->next)
        len += sprintf(&out[len], "%s%s ", opt ? opt : "", tmp->v);
    }
  else out=NULL;

  /* Return the final string. */
  return out;
}





int
buildprog(struct buildprogparams *p)
{
  /* Note that the first node of 'sourceargs' is the acutal source and the
     rest are arguments to be run later. */
  int retval;
  char *fullla;
  char *ldflags=NULL, *cppflags=NULL;
  char *command, *optimize=NULL, *warning=NULL;
  char *include   = buildprog_as_one_string("-I", p->include);
  char *linkdir   = buildprog_as_one_string("-L", p->linkdir);
  char *linklib   = buildprog_as_one_string("-l", p->linklib);
  char *arguments = buildprog_as_one_string(NULL, p->sourceargs->next);

  /* If not in quiet mode, let the user know. */
  if(!p->cp.quiet)
    {
      printf("\nCompiling and linking the program\n");
      printf("---------------------------------\n");
    }

  /* If environment should be read, read it. */
  if(p->noenv==0)
    {
      ldflags=getenv("LDFLAGS");
      cppflags=getenv("CPPFLAGS");
    }

  /* Compiler options with values: */
  if(p->warning)
    if( asprintf(&warning,  "-W%s", p->warning)<0 )
      error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
  if(p->optimize)
    if( asprintf(&optimize, "-O%s", p->optimize)<0 )
      error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);

  /* Libtool '.la' file: */
  if(p->la) fullla=p->la;
  else
    if( asprintf(&fullla, "%s/libgnuastro.la", LIBDIR)<0 )
      error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);

  /* Write the full Libtool command into a string (to run afterwards). */
  if( asprintf(&command, "%s -c \"%s %s %s%s --mode=link %s %s %s "
               "%s %s %s %s %s %s %s -I%s %s -o %s\"",
               GAL_CONFIG_GNULIBTOOL_SHELL,
               GAL_CONFIG_GNULIBTOOL_EXEC,
               p->cp.quiet ? "--quiet" : "",
               p->tag      ? "--tag="  : "",
               p->tag      ? p->tag    : "",
               p->cc,
               warning     ? warning   : "",
               p->debug    ? "-g"      : "",
               optimize    ? optimize  : "",
               include     ? include   : "",
               cppflags    ? cppflags  : "",
               linkdir     ? linkdir   : "",
               ldflags     ? ldflags   : "",
               p->sourceargs->v,
               linklib     ?linklib    : "",
               INCLUDEDIR,
               fullla,
               p->cp.output)<0 )
    error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);

  /* Compile (and link): */
  retval=system(command);
  if(retval!=EXIT_SUCCESS)
    error(EXIT_FAILURE, 0, "failed to build, see libtool error above");
  else if(p->onlybuild==0)
    {
      /* Free the initial command. */
      free(command);

      /* Wright the command to run the program. Note that if the output
         value doesn't start with a directory, we'll have to put one for
         it. */
      switch(p->cp.output[0])
        {
        case '.':
        case '/':
          if( asprintf(&command, "%s %s", p->cp.output,
                       arguments?arguments:"")<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
          break;

        default:
          if( asprintf(&command, "./%s %s", p->cp.output,
                       arguments?arguments:"")<0 )
            error(EXIT_FAILURE, 0, "%s: asprintf allocation", __func__);
        }

      /* Print the executed command if necessary, then run it. */
      if(!p->cp.quiet)
        {
          printf("\nRun the compiled program\n");
          printf("------------------------\n");
          printf("%s\n", command);
        }
      retval=system(command);

      /* Delete the compiled program after running it. */
      if(p->deletecompiled)
        {
          errno=0;
          if( remove(p->cp.output) == -1 )
            error(EXIT_FAILURE, 0, "unable to delete %s", p->cp.output);
        }
    }

  /* Clean up and return. */
  free(include);
  free(linkdir);
  free(linklib);
  free(command);
  if(warning) free(warning);
  if(optimize) free(optimize);
  if(!p->la && fullla) free(fullla);
  return retval;
}
