/*********************************************************************
Common parameters for all the utilities.
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
#ifndef __GAL_COMMONPARAMS_H__
#define __GAL_COMMONPARAMS_H__



/* The structure keeping all the common parameters in gnuastro. This
   could be part of commonargs.h, but since main.h in all those
   programs will need this structure, there will be problems about the
   unused `commonargp` structure.*/
struct gal_commonparams
{
  char        *spack;  /* Subpackage name.                          */

  /* Input/Output: */
  char       *output;  /* Directory containg output.                */
  char          *hdu;  /* Image extension.                          */
  int     dontdelete;  /* ==1: Don't delete existing.               */
  int  removedirinfo;  /* ==1: Remove directory information.        */

  /* Operating modes: */
  int           verb;  /* ==1: report steps. ==0 don't.             */
  int    printparams;  /* Only print the used values.               */
  int     setdirconf;  /* ==1: Set the current directory config.    */
  int     setusrconf;  /* ==1: Set the user default values.         */
  size_t  numthreads;  /* Number of threads to use.                 */
  int    onlydirconf;  /* Only check current directory conf. file.  */
  char  *onlyversion;  /* The string of the requested version.      */
  size_t  minmapsize;  /* The minimum bytes necessary to use mmap.  */
  int          nolog;  /* ==1: do not make a log file.              */

  /* Check: */
  int       quietset;  /* If the verbose flag is called.            */
  int  numthreadsset;  /* If the number of threads are set.         */
  int onlyversionset;  /* If the only version option is set.        */
  int onlydirconfset;  /* If --onlydirconf was set before this.     */
  int         hduset;  /* If the input image extension is set.      */
  int      outputset;  /* If the output is set.                     */
  int       nologset;  /* If nolog is set.                          */
  int  minmapsizeset;  /* If minmapsize is set.                     */
  int  dontdeleteset;  /* If the --dontdelete option was called.    */
  int removedirinfoset;  /* If --keepinputdir was called.           */
};

#endif
