/*********************************************************************
BuildProgram: Compile and run programs using Gnuastro's library
BuildProgram is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2017-2018, Free Software Foundation, Inc.

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
      /* For every node in the list, we want the `opt' and a space along with
         the actual string. */
      for(tmp=list; tmp!=NULL; tmp=tmp->next)
        len += 1 + (opt ? strlen(opt) : 0) + strlen(tmp->v);

      /* Allocate space for the string. */
      out=gal_data_malloc_array(GAL_TYPE_UINT8, len+1, __func__, "out");

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
  /* Note that the first node of `sourceargs' is the acutal source and the
     rest are arguments to be run later. */
  int retval;
  char *fullla;
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

  /* Compiler options with values: */
  if(p->warning)   asprintf(&warning,  "-W%s", p->warning);
  if(p->optimize)  asprintf(&optimize, "-O%s", p->optimize);

  /* Libtool `.la' file: */
  if(p->la) fullla=p->la;
  else      asprintf(&fullla, "%s/libgnuastro.la", LIBDIR);

  /* Put the command to run into a string. */
  asprintf(&command, "%s %s %s%s --mode=link gcc %s %s %s %s %s %s %s "
           "-I%s %s -o %s",
           GAL_CONFIG_GNULIBTOOL_EXEC,
           p->cp.quiet ? "--quiet" : "",
           p->tag      ? "--tag="   : "",
           p->tag      ? p->tag    : "",
           warning     ? warning   : "",
           p->debug    ? "-g"      : "",
           optimize    ? optimize  : "",
           include     ? include   : "",
           linkdir     ? linkdir   : "",
           p->sourceargs->v,
           linklib     ?linklib    : "",
           INCLUDEDIR,
           fullla,
           p->cp.output);

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
          asprintf(&command, "%s %s", p->cp.output, arguments?arguments:"");
          break;

        default:
          asprintf(&command, "./%s %s", p->cp.output, arguments?arguments:"");
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
