/*********************************************************************
configfiles -- Read configuration files for each program.
This is part of GNU Astronomy Utilities (gnuastro) package.

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
along with gnuastro. If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/
#ifndef CONFIGFILES_H
#define CONFIGFILES_H



/**************************************************************/
/************               Macros                *************/
/**************************************************************/
/* Simple macros: */
#define CONFIG_DELIMITERS " ,=:\t\n"



#define STARTREADINGLINE {						\
  ++lineno;								\
  if(*line=='#') continue;						\
  else readnamevalue(line, filename, lineno, &name, &value);		\
  if(name==NULL && value==NULL) continue;				\
}



/* Functional macros: These are not actual functions, because they
   depend on functions that are different for different programs. So
   they have to be written into the functions with a macro. */
#define SAVE_LOCAL_CONFIG(INDIR) {					\
    FILE *fp;								\
    char *outfilename, *command;		 			\
    fp=writelocalconfigstop(INDIR, CONFIG_FILE, SPACK,			\
			     SPACK_NAME, &outfilename);			\
    printvalues(fp, p);							\
    errno=0;								\
    if(fclose(fp)==-1)							\
      error(EXIT_FAILURE, errno, "%s", outfilename);                    \
    command=malloccat("cat ", outfilename);				\
    printf("Values saved in %s:\n\n", outfilename);			\
    if(system(command))                                                 \
      error(EXIT_FAILURE, 0, "The `%s` command could not be run or "    \
            "failed.", command);                                        \
    free(outfilename);							\
    free(command);							\
    exit(EXIT_SUCCESS);							\
  }





#define CHECKSETCONFIG {				                \
    char *userconfig_dir, *userconfig_file;				\
    userconfig_dir=addhomedir(USERCONFIG_DIR);				\
    userconfig_file=addhomedir(USERCONFIG_FILEEND);			\
    if(cp->setdirconf || cp->setusrconf)				\
      {									\
	if(cp->setdirconf)						\
	  {								\
	    readconfig(CURDIRCONFIG_FILE, p);				\
	    SAVE_LOCAL_CONFIG(CURDIRCONFIG_DIR);			\
	  }								\
	if(cp->setusrconf)						\
	  {								\
	    readconfig(userconfig_file, p);				\
	    SAVE_LOCAL_CONFIG(userconfig_dir);				\
	  }								\
      }									\
    else								\
      {									\
	readconfig(CURDIRCONFIG_FILE, p);			        \
	readconfig(userconfig_file, p);					\
	readconfig(SYSCONFIG_FILE, p);					\
      }									\
    free(userconfig_file);						\
    free(userconfig_dir);						\
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
	char *userconfig_file;						\
	fprintf(stderr, ".\n\n");					\
	fprintf(stderr, "You can assign values in the local, user or "	\
		"system wide default files. Otherwise you have to "	\
		"explicitly call them each time. See `"SPACK" --help` "	\
		"or `info "SPACK"` for more information.\n\n");		\
	userconfig_file=addhomedir(USERCONFIG_FILEEND);			\
	fprintf(stderr, "Default files checked (existing or not):\n"	\
		"   %s\n   %s\n   %s\n", CURDIRCONFIG_FILE,		\
		userconfig_file, SYSCONFIG_FILE);			\
	free(userconfig_file);						\
	exit(EXIT_FAILURE);						\
      }									\
  }

#define REPORT_PARAMETERS_SET {			                        \
    fprintf(stdout, "# "SPACK_STRING"\n");				\
    fprintf(stdout, "# Configured on "CONFIGDATE" at "CONFIGTIME"\n");	\
    fprintf(stdout, "# Written on %s", ctime(&p->rawtime));	        \
    printvalues(stdout, p);						\
    exit(EXIT_SUCCESS);							\
  }









/**************************************************************/
/************       Function declarations         *************/
/**************************************************************/
char *
addhomedir(char *dir);

void
readnamevalue(char *line, char *filename, size_t lineno,
	      char **name, char **value);

FILE *
writelocalconfigstop(char *indir, char *filename, char *spack,
		     char *spack_name, char **outfilename);

#endif
