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
#ifndef __GAL_CONFIGFILES_H__
#define __GAL_CONFIGFILES_H__

/* Include other headers if necessary here. Note that other header files
   must be included before the C++ preparations below */



/* C++ Preparations */
#undef __BEGIN_C_DECLS
#undef __END_C_DECLS
#ifdef __cplusplus
# define __BEGIN_C_DECLS extern "C" {
# define __END_C_DECLS }
#else
# define __BEGIN_C_DECLS                /* empty */
# define __END_C_DECLS                  /* empty */
#endif
/* End of C++ preparations */



/* Actual header contants (the above were for the Pre-processor). */
__BEGIN_C_DECLS  /* From C++ preparations */



/**************************************************************/
/************               Macros                *************/
/**************************************************************/
/* Simple macros: */
#define GAL_CONFIGFILES_DELIMITERS " ,=:\t\n"



#define GAL_CONFIGFILES_START_READING_LINE {                            \
    ++lineno;                                                           \
  if(*line=='#') continue;                                              \
  else gal_configfiles_read_name_value(line, filename, lineno,          \
                                       &name, &value);                  \
  if(name==NULL && value==NULL) continue;                               \
}



/* Functional macros: These are not actual functions, because they
   depend on functions that are different for different programs. So
   they have to be written into the functions with a macro. */
#define GAL_CONFIGFILES_SAVE_LOCAL_CONFIG(INDIR) {                      \
    FILE *fp;                                                           \
    char *outfilename, *command;                                        \
    fp=gal_configfiles_write_local_config_stop(INDIR, CONFIG_FILE,      \
                                               SPACK, SPACK_NAME,       \
                                               &outfilename);           \
    printvalues(fp, p);                                                 \
    errno=0;                                                            \
    if(fclose(fp)==-1)                                                  \
      error(EXIT_FAILURE, errno, "%s", outfilename);                    \
    command=gal_checkset_malloc_cat("cat ", outfilename);               \
    printf("Values saved in %s:\n\n", outfilename);                     \
    if(system(command))                                                 \
      error(EXIT_FAILURE, 0, "the `%s` command could not be run or "    \
            "failed", command);                                         \
    free(outfilename);                                                  \
    free(command);                                                      \
    exit(EXIT_SUCCESS);                                                 \
  }





#define GAL_CONFIGFILES_CHECK_SET_CONFIG {                              \
    char *userconfig_dir, *userconfig_file;                             \
                                                                        \
    readconfig(CURDIRCONFIG_FILE, p);                                   \
    if(cp->setdirconf)                                                  \
      GAL_CONFIGFILES_SAVE_LOCAL_CONFIG(CURDIRCONFIG_DIR);              \
    if(cp->onlyversionset && strcmp(cp->onlyversion, PACKAGE_VERSION))  \
      error(EXIT_FAILURE, 0, "you are currently running Gnuastro %s. "  \
            "However, this run should be with version `%s'.\n\n"        \
            "To resolve the situation, use the the '--onlyversion' "    \
            "option, either on the command-line or in a configuration " \
            "file. For example, set it to `%s' by repeating the "       \
            "previous command with:\n\n"                                \
            "    --onlyversion=%s\n\n"                                  \
            "Alternatively, you can install Gnuastro %s.\n\n"           \
            "NOTE: If you didn't set this option on the command-line, " \
            "it was probably intended for reproducability. If so, it "  \
            "is advised to install Gnuastro %s",                        \
            PACKAGE_VERSION, cp->onlyversion, PACKAGE_VERSION,          \
            PACKAGE_VERSION, cp->onlyversion, cp->onlyversion);         \
                                                                        \
    if(cp->onlydirconf==0)                                              \
      {                                                                 \
        userconfig_dir=gal_configfiles_add_home_dir(USERCONFIG_DIR);    \
        userconfig_file=                                                \
          gal_configfiles_add_home_dir(USERCONFIG_FILEEND);             \
        readconfig(userconfig_file, p);                                 \
        if(cp->setusrconf)                                              \
          GAL_CONFIGFILES_SAVE_LOCAL_CONFIG(userconfig_dir);            \
        readconfig(SYSCONFIG_FILE, p);                                  \
        free(userconfig_file);                                          \
        free(userconfig_dir);                                           \
      }                                                                 \
  }







#define GAL_CONFIGFILES_REPORT_NOTSET(var_name) {                       \
    if(intro==0)                                                        \
      {                                                                 \
        fprintf(stderr, SPACK": Parameter(s) not set: `%s'",            \
                (var_name));                                            \
        intro=1;                                                        \
      }                                                                 \
    else                                                                \
      fprintf(stderr, ", `%s'", (var_name));                            \
  }





#define GAL_CONFIGFILES_END_OF_NOTSET_REPORT {                          \
    if(intro)                                                           \
      {                                                                 \
        char *userconfig_file;                                          \
        fprintf(stderr, ".\n\n");                                       \
        fprintf(stderr, "You can assign values in the local, user or "  \
                "system wide default files. Otherwise you have to "     \
                "explicitly assign a value to them each time as a "     \
                "command-line option. See `%s --help` or `info %s` for "\
                "more information.\n\n", SPACK, SPACK);                 \
        userconfig_file=                                                \
          gal_configfiles_add_home_dir(USERCONFIG_FILEEND);             \
        fprintf(stderr, "Default files checked (existing or not):\n");  \
        fprintf(stderr, "   %s\n", CURDIRCONFIG_FILE);                  \
        if(p->cp.onlydirconf==0)                                        \
          fprintf(stderr, "   %s\n   %s\n", userconfig_file,            \
                  SYSCONFIG_FILE);                                      \
        free(userconfig_file);                                          \
        exit(EXIT_FAILURE);                                             \
      }                                                                 \
  }





#define GAL_CONFIGFILES_REPORT_PARAMETERS_SET {                         \
    fprintf(stdout, "# "SPACK_STRING"\n");                              \
    fprintf(stdout, "# Written on %s", ctime(&p->rawtime));             \
    printvalues(stdout, p);                                             \
    exit(EXIT_SUCCESS);                                                 \
  }





/* Read the options that are common to all programs from the
   configuration file. Since these two checks are within an if-else
   structure, they should not be placed within an `{' and `}'. */
#define GAL_CONFIGFILES_READ_COMMONOPTIONS_FROM_CONF                    \
    else if(strcmp(name, "quiet")==0)                                   \
      {                                                                 \
        int tint;                                                       \
        if(cp->quietset) continue;                                      \
        gal_checkset_int_zero_or_one(value, &tint, name, key, SPACK,    \
                                     filename, lineno);                 \
        cp->verb=!tint;                                                 \
        cp->quietset=1;                                                 \
      }                                                                 \
    else if(strcmp(name, "numthreads")==0)                              \
      {                                                                 \
        if(cp->numthreadsset) continue;                                 \
        gal_checkset_sizet_l_zero(value, &cp->numthreads, name, key,    \
                                  SPACK, filename, lineno);             \
        cp->numthreadsset=1;                                            \
      }                                                                 \
    else if(strcmp(name, "onlydirconf")==0)                             \
      {                                                                 \
        if(cp->onlydirconfset) continue;                                \
        gal_checkset_int_zero_or_one(value, &cp->onlydirconf, name,     \
                                     key, SPACK, filename, lineno);     \
        cp->onlydirconfset=1;                                           \
      }                                                                 \
    else if(strcmp(name, "onlyversion")==0)                             \
      {                                                                 \
        if(cp->onlyversionset) continue;                                \
        gal_checkset_allocate_copy_set(value, &cp->onlyversion,         \
                                       &cp->onlyversionset);            \
        cp->onlyversionset=1;                                           \
      }                                                                 \
    else if(strcmp(name, "nolog")==0)                                   \
      {                                                                 \
        if(cp->nologset) continue;                                      \
        gal_checkset_int_zero_or_one(value, &cp->nolog, name, key,      \
                                     SPACK, filename, lineno);          \
        cp->nologset=1;                                                 \
      }                                                                 \
    else if(strcmp(name, "minmapsize")==0)                              \
      {                                                                 \
        if(cp->minmapsizeset) continue;                                 \
        gal_checkset_sizet_l_zero(value, &cp->minmapsize, name, key,    \
                                  SPACK, filename, lineno);             \
        cp->minmapsizeset=1;                                            \
      }                                                                 \
                                                                        \
                                                                        \
    else if(strcmp(name, "dontdelete")==0)                              \
      {                                                                 \
        if(cp->dontdeleteset) continue;                                 \
        gal_checkset_int_zero_or_one(value, &cp->dontdelete, name, key, \
                                     SPACK, filename, lineno);          \
        cp->dontdeleteset=1;                                            \
      }                                                                 \
    else if(strcmp(name, "keepinputdir")==0)                            \
      {                                                                 \
        int tint;                                                       \
        if(cp->removedirinfoset) continue;                              \
        gal_checkset_int_zero_or_one(value, &tint, name, key,           \
                                     SPACK, filename, lineno);          \
        cp->removedirinfo=!tint;                                        \
        cp->removedirinfoset=1;                                         \
      }                                                                 \







/* Write common options: */
#define GAL_CONFIGFILES_PRINT_COMMONOPTIONS {                           \
    if(cp->quietset)                                                    \
      fprintf(fp, CONF_SHOWFMT"%d\n", "quiet", !p->cp.verb);            \
    if(cp->numthreadsset)                                               \
      fprintf(fp, CONF_SHOWFMT"%zu\n", "numthreads", p->cp.numthreads); \
    if(cp->onlydirconfset)                                              \
      fprintf(fp, CONF_SHOWFMT"%d\n", "onlydirconf", p->cp.onlydirconf);\
    if(cp->onlyversionset)                                              \
      GAL_CHECKSET_PRINT_STRING_MAYBE_WITH_SPACE("onlyversion",         \
                                                 cp->onlyversion);      \
    if(cp->nologset)                                                    \
      fprintf(fp, CONF_SHOWFMT"%d\n", "nolog", p->cp.nolog);            \
    if(cp->dontdeleteset)                                               \
      fprintf(fp, CONF_SHOWFMT"%d\n", "dontdelete", p->cp.dontdelete);  \
    if(cp->removedirinfoset)                                            \
      fprintf(fp, CONF_SHOWFMT"%d\n", "keepinputdir", !p->cp.removedirinfo); \
  }








/**************************************************************/
/************       Function declarations         *************/
/**************************************************************/
char *
gal_configfiles_add_home_dir(char *dir);

void
gal_configfiles_read_name_value(char *line, char *filename, size_t lineno,
                                char **name, char **value);

FILE *
gal_configfiles_write_local_config_stop(char *indir, char *filename,
                                        char *spack, char *spack_name,
                                        char **outfilename);
void
gal_configfiles_print_type(FILE *fp, int bitpix);



__END_C_DECLS    /* From C++ preparations */

#endif           /* __GAL_CONFIGFILES_H__ */
