/*********************************************************************
forall_defaults -- Similar functions for defaults in all utilities.
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
#ifndef FORALL_DEFAULTS_H
#define FORALL_DEFAULTS_H



/**************************************************************/
/************               Macros                *************/
/**************************************************************/
/* Simple macros: */
#define DEFAULT_DELIMITERS " ,=:\t\n"




/* Functional macros: These are not actual functions, because they
   depend on functions that are different for different programs. So
   they have to be written into the functions with a macro. */
#define SAVE_LOCAL_DEFAULTS(INDIR) {					\
    FILE *fp;								\
    int freeindir=0;							\
    char *indir=(INDIR), *outfilename, *command;			\
    if(strcmp(indir, CURDIRDATA_DIR))					\
      {freeindir=1; indir=addhomedir(INDIR);}				\
    fp=writelocaldefaultstop(indir, DEFAULT_FILE, SPACK,		\
			     SPACK_NAME, &outfilename);			\
    if(freeindir) free(indir);						\
    printvalues(fp, p);							\
    errno=0;								\
    if(fclose(fp)==-1)							\
      error(EXIT_FAILURE, errno, outfilename);				\
    command=malloccat("cat ", outfilename);				\
    printf("Default values saved in %s:\n\n", outfilename);		\
    system(command);							\
    free(outfilename);							\
    free(command);							\
    exit(EXIT_SUCCESS);							\
  }





#define CHECKSETDEFAULTS {				                \
    char *userdata_dir, *userdefault_file;				\
									\
    readdefaults(CURDIRDEFAULT_FILE, p);				\
    if(cp->dirdefaults)							\
      SAVE_LOCAL_DEFAULTS(CURDIRDATA_DIR);				\
									\
    userdata_dir=addhomedir(USERDATA_DIR);				\
    userdefault_file=addhomedir(USERDEFAULT_FILEEND);			\
    readdefaults(userdefault_file, p);					\
    if(cp->userdefaults)						\
      SAVE_LOCAL_DEFAULTS(userdata_dir);				\
    free(userdefault_file);						\
    free(userdata_dir);							\
    									\
    readdefaults(SYSDEFAULT_FILE, p);					\
  }







#define REPORT_NOTSET(var_name) {					\
    if(intro==0)							\
      {									\
	fprintf(stderr, SPACK": Parameter(s) not set: %s", (var_name));	\
	intro=1;							\
      }									\
    else								\
      fprintf(stderr, ", %s", (var_name));				\
  }





#define END_OF_NOTSET_REPORT {			                        \
    if(intro)								\
      {									\
	char *userdefault_file;						\
	fprintf(stderr, ".\n\n");					\
	fprintf(stderr, "You can assign values in the local, user or "	\
		"system wide default files. Otherwise you have to "	\
		"explicitly call them each time. See `"SPACK" --help` "	\
		"or `info "SPACK"` for more information.\n\n");		\
	userdefault_file=addhomedir(USERDEFAULT_FILEEND);		\
	fprintf(stderr, "Default files checked (existing or not):\n"	\
		"   %s\n   %s\n   %s\n", CURDIRDEFAULT_FILE,		\
		userdefault_file, SYSDEFAULT_FILE);			\
	free(userdefault_file);						\
	exit(EXIT_FAILURE);						\
      }									\
  }

#define REPORT_PARAMETERS_SET {			                        \
    fprintf(stdout, SPACK_STRING"\n");				\
    fprintf(stdout, "Configured on "CONFIGDATE" at "CONFIGTIME"\n");	\
    printvalues(stdout, p);						\
    exit(EXIT_SUCCESS);							\
  }









/**************************************************************/
/************       Function declarations         *************/
/**************************************************************/
char *
addhomedir(char *dir);

FILE *
writelocaldefaultstop(char *indir, char *filename, char *spack,
		      char *spack_name, char **outfilename);

#endif
