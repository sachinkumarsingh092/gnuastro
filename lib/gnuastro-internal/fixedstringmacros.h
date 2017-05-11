/*********************************************************************
Fixed strings for use in all utilities.
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
#ifndef __GAL_FIXEDSTRINGMACROS_H__
#define __GAL_FIXEDSTRINGMACROS_H__



#define GAL_STRINGS_SHORT_COPYRIGHT                                     \
  "Copyright (C) 2015-2017, Free Software Foundation, Inc."


#define GAL_STRINGS_SHORT_LICENSE                                       \
  "License GPLv3+: GNU General public license version 3 or later."


#define GAL_STRINGS_COPYRIGHT                                           \
  GAL_STRINGS_SHORT_COPYRIGHT"\n"GAL_STRINGS_SHORT_LICENSE"\n"          \
  "This is free software: you are free to change and redistribute "     \
  "it.\nThere is NO WARRANTY, to the extent permitted by law."          \


#define GAL_STRINGS_TOP_HELP_INFO                                       \
  "\n"PROGRAM_NAME" is part of "PACKAGE_STRING".\n"



/* This is fixed for all the packages. */
#define GAL_STRINGS_MORE_HELP_INFO                                      \
  "\nFor more information, please run any of the "                      \
  "following commands. They will respectively show you the `Invoking "  \
  PROGRAM_NAME"' subsection, the complete `"PROGRAM_NAME"' section, or "\
  "the full "PACKAGE_NAME" manual. In particular the first contains a " \
  "very comprehensive explanation of "PROGRAM_NAME"'s invocation and "  \
  "all the options.\n\n"                                                \
  "     $ info "PROGRAM_EXEC"\n\n"                                      \
  "     $ info "PROGRAM_NAME"\n\n"                                      \
  "     $ info "PACKAGE_TARNAME"\n\n"                                   \
  "If you couldn't find your answer in the manual, you can get "        \
  "direct help from experienced Gnuastro users and developers. "        \
  "For more information, please run:\n\n"                               \
  "     $ info help-gnuastro\n\n"                                       \
  PROGRAM_NAME" options:"                                               \


/* This can be used in the end of error messages related to option
   values. */
#define GAL_STRINGS_HOW_TO_CHECK_VALUES                                 \
  " You can check all the input values with the `--printparams' "       \
  "(-P) option."



#endif           /* __GAL_FIXEDSTRINGMACROS_H__ */
