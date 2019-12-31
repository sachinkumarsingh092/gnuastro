/*********************************************************************
Functions to deal with Git version controlled directories.
This is part of GNU Astronomy Utilities (Gnuastro) package.

Original author:
     Mohammad Akhlaghi <mohammad@akhlaghi.org>
Contributing author(s):
Copyright (C) 2016-2020, Free Software Foundation, Inc.

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

#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/git.h>

/* Return the result of Git describe and  */
char *
gal_git_describe(void)
{
  char *describe_return=NULL;

/* Only actually work on Git if the library is present. */
#if GAL_CONFIG_HAVE_LIBGIT2 == 1
  git_buf buf={0};
  git_repository *repo;
  git_describe_result *describe_result;
  git_describe_format_options format_options;
  git_describe_options describe_options = GIT_DESCRIBE_OPTIONS_INIT;

  /* Initialize Git's global setup. */
  git_libgit2_init();

  /* Initialize the descriptions. */
  git_describe_init_format_options(&format_options,
				   GIT_DESCRIBE_FORMAT_OPTIONS_VERSION);

  /* Set some values manually: */
  describe_options.show_commit_oid_as_fallback=1;
  format_options.dirty_suffix="-dirty";

  /* Open the Git repository. */
  if( !git_repository_open_ext(&repo, ".", 0, NULL)
      && !git_describe_workdir(&describe_result, repo, &describe_options)
      && !git_describe_format(&buf, describe_result, &format_options) )
    {
      /* Allocate space for and copy the description. */
      errno=0;
      describe_return=malloc(strlen(buf.ptr)+1);
      if(describe_return==NULL)
        error(EXIT_FAILURE, errno, "%s: allocating %zu bytes to copy "
              "Git's describe", __func__, strlen(buf.ptr)+1);
      strcpy(describe_return, buf.ptr);
    }

  /* Clean up. */
  git_repository_free(repo);
  git_libgit2_shutdown();
#endif

  /* Return */
  return describe_return;
}
