/*********************************************************************
Function to parse options and configuration file values.

Original author:
     Mohammad Akhlaghi <akhlaghi@gnu.org>
Contributing author(s):
Copyright (C) 2017, Free Software Foundation, Inc.

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

#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/txt.h>
#include <gnuastro/data.h>
#include <gnuastro/linkedlist.h>

#include <options.h>










/**********************************************************************/
/************             Option utilities              ***************/
/**********************************************************************/
int
gal_options_is_last(struct argp_option *option)
{
  return ( option->key==0 && option->name==NULL
           && option->doc==NULL && option->group==0 );
}




int
gal_options_is_category_title(struct argp_option *option)
{
  return ( option->key==0 && option->name==NULL );
}




/* The modified `argp_option' structure contains a void pointer, depending
   on the type of value inside it, you need to the pointer differently. */
void
gal_options_free(struct argp_option *options)
{
  size_t i;

  /* Go over all the options and free those that are allocated. After
     freeing, set them to NULL, so they are not mistakenly used later. */
  for(i=0; !gal_options_is_last(&options[i]); ++i)
    if(options[i].value)
      {
        if(options[i].type==GAL_DATA_TYPE_STRLL)
          gal_linkedlist_free_stll(options[i].value, 1);
        else
          free(options[i].value);
        options[i].value=NULL;
      }
}




















/**********************************************************************/
/************            Command-line options           ***************/
/**********************************************************************/
/* Set the value given to the command-line, where we have the integer `key'
   of the option, not its long name as a string. */
error_t
gal_options_set_from_key(int key, char *arg, struct argp_option *options)
{
  size_t i;
  char **strarr=NULL;

  /* Go through all the options and find the one that should keep this
     value, then put its value into the appropriate key. Note that the
     options array finishs with an all zero element, so we don't need to
     know the number before hand.*/
  for(i=0;1;++i)
    {
      /* Check if the key corresponds to this option. */
      if(options[i].key==key)
        {
          /* When options are read from keys, they are read from the
             command-line. On the commandline, the last invokation of the
             option is important. Especially in contexts like scripts, this
             is important because you can change a given command-line
             option (that is not a linked list) by calling it a second
             time, instead of going back and changing the first value.

             As a result, only when searching for options on the
             command-line, a second value to the same option will replace
             the first one. This will not happen in configuration files. */
          if(options[i].value && gal_data_is_linked_list(options[i].type)==0)
            {
              free(options[i].value);
              options[i].value=NULL;
            }

          /* We have two types of options: those which need an argument and
             those that don't. When they need an argument the `value'
             pointer will be NULL. In such cases, the option must be in
             integer type, so just set it to 1. When an argument is given,
             convert the given value to the appropriate type and put it in
             the `value' element of options[i]. */
          if(arg)
            {
              /* For strings, `gal_data_string_to_type' is going to return
                 an allocated pointer to an allocated string (`char
                 **'). In this context, we are just dealing with one
                 string, and arrays of strings are never used (a linked
                 list is defined when multiple strings must be read). So
                 only keep the actual string and free the one that kept
                 it. */
              if(options[i].type==GAL_DATA_TYPE_STRING)
                {
                  gal_data_string_to_type((void **)(&strarr), arg,
                                                 options[i].type);
                  options[i].value=strarr[0];
                  free(strarr);
                }
              else
                {
                  if( gal_data_string_to_type(&options[i].value, arg,
                                              options[i].type) )
                    error(EXIT_FAILURE, 0, "`%s' (value to option `%s') "
                          "couldn't be read as a number", arg,
                          options[i].name);

                }
            }
          else
            {
              /* Make sure the option has the type set for options with no
                 argument. So, give it a value of 1. */
              if(options[i].type==GAL_OPTIONS_NO_ARG_TYPE)
                gal_data_string_to_type(&(options[i].value), "1",
                                        options[i].type);
              else
                error(EXIT_FAILURE, 0, "A bug! Please contact us at %s to "
                      "correct it. Options with no arguments, must have "
                      "type `%s' to be read in `gal_options_read_from_key'. "
                      "However, the option with key `%d' has type %s",
                      PACKAGE_BUGREPORT,
                      gal_data_type_as_string(GAL_OPTIONS_NO_ARG_TYPE, 1),
                      key, gal_data_type_as_string(options[i].type, 1));
            }


          /* We have found and set the value given to this option, so just
             return success (an error_t of 0 means success). */
          return 0;
        }
      else
        {
          /* The last option has all its values set to zero. */
          if(gal_options_is_last(&options[i]))
            return ARGP_ERR_UNKNOWN;
        }
    }
}





error_t
gal_options_common_argp_parse(int key, char *arg, struct argp_state *state)
{
  /* In case the user incorrectly uses the equal sign (for example
     with a short format or with space in the long format, then `arg`
     start with (if the short version was called) or be (if the long
     version was called with a space) the equal sign. So, here we
     check if the first character of arg is the equal sign, then the
     user is warned and the program is stopped: */
  if(arg && arg[0]=='=')
    argp_error(state, "incorrect use of the equal sign (`=`). For short "
               "options, `=` should not be used and for long options, "
               "there should be no space between the option, equal sign "
               "and value");

  /*********************************************************************/
  /*********************************************************************/
  /*********************************************************************


   -----> Do this check before passing control to the program. <-----


  if( key==ARGP_KEY_END )
    {
      if(cp->setdirconf && cp->setusrconf)
        error(EXIT_FAILURE, 0, "only one of `--setusrconf` or "
              "`--setdirconf` may be set in each run. You have asked "
              "for both");
    }
   *********************************************************************/
  /*********************************************************************/
  /*********************************************************************/

  return gal_options_set_from_key(key, arg, gal_commonopts_options);
}




















/**********************************************************************/
/************            Configuration files            ***************/
/**********************************************************************/

/* Read the option and the argument from the line and return.*/
static void
options_read_name_arg(char *line, char *filename, size_t lineno,
                        char **name, char **arg)
{
  int notyetfinished=1, inword=0, inquote=0;

  /* Initialize name and value: */
  *arg=NULL;
  *name=NULL;

  /* Go through the characters and set the values: */
  do
    switch(*line)
      {
      case ' ': case '\t': case '\v': case '\n': case '\r':
        if(inword) /* Only considered in a word, not in a quote*/
          {
            inword=0;
            *line='\0';
            if(*arg && inquote==0)
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
                            "option name should not start with "
                            "double quotes (\").");
            inquote=1;
            *arg=line+1;
          }
        break;
      default:
        if(inword==0 && inquote==0)
          {
            if(*name==NULL)
              *name=line;
            else  /* *name is set, now assign *arg. */
              *arg=line;
            inword=1;
          }
        break;
      }
  while(*(++line)!='\0' && notyetfinished);

  /* In the last line of the file, there is no new line to be
     converted to a '\0' character! So if value has been assigned, we
     are not in a quote and the line has finished, it means the given
     value has also finished. */
  if(*line=='\0' && *arg && inquote==0)
    notyetfinished=0;

  /* This was a blank line: */
  if(*name==NULL && *arg==NULL)
    return;

  /* Name or value were set but not yet finished. */
  if(notyetfinished)
    error_at_line(EXIT_FAILURE, 0, filename, lineno,
                  "line finished before option name and value could "
                  "be read.");
}





static int
options_set_from_name(char *name, char *arg, struct argp_option *options,
                      char *filename, size_t lineno)
{
  size_t i;
  char **strarr=NULL;

  /* Go through all the options and find the one that should keep this
     value, then put its value into the appropriate key. Note that the
     options array finishs with an all zero element, so we don't need to
     know the number before hand.*/
  for(i=0;1;++i)
    {
      /* Check if the key corresponds to this option. */
      if(options[i].name && !strcmp(options[i].name, name))
        {
          /* If the option already has a value and it isn't a linked
             list, then ignore it. */
          if(options[i].value && !gal_data_is_linked_list(options[i].type))
            return 0;

          /* For strings, `gal_data_string_to_type' is going to return an
             allocated pointer to an allocated string (`char **'). In this
             context, we are just dealing with one string, and arrays of
             strings are never used (a linked list is defined when multiple
             strings must be read). So only keep the actual string and free
             the one that kept it. */
          if(options[i].type==GAL_DATA_TYPE_STRING)
            {
              gal_data_string_to_type((void **)(&strarr), arg,
                                      options[i].type);
              options[i].value=strarr[0];
              free(strarr);
            }
          else
            {
              if( gal_data_string_to_type(&options[i].value, arg,
                                          options[i].type) )
                error_at_line(EXIT_FAILURE, 0, filename, lineno,
                              "`%s' (value to option `%s') couldn't be "
                              "read as a number", arg, options[i].name);
            }


          /* If this is an on/off option (with no argument), then check if
             the given value is 0 or 1. */
          if( options[i].type==GAL_OPTIONS_NO_ARG_TYPE
              && *(unsigned char *)(options[i].value) > 1 )
            error_at_line(EXIT_FAILURE, 0, filename, lineno,
                          "`%s' is an on/off option, so its value can only "
                          "be 1 (for `on'), or 0 (for `off'), it was given "
                          "a value of `%s'", options[i].name, arg);


          /* We have found and set the value given to this option, so just
             return success (an error_t of 0 means success). */
          return 0;
        }
      else
        {
          /* The last option has all its values set to zero. */
          if(gal_options_is_last(&options[i]))
            return 1;
        }
    }
}





static void
options_parse_file(char *filename,  struct argp_option *poptions,
                   struct argp_option *coptions)
{
  FILE *fp;
  char *line, *name, *arg;
  size_t linelen=10, lineno=0;


  /* Open the file. If the file doesn't exist, then just ignore the
     configuration file and return. */
  errno=0;
  fp=fopen(filename, "r");
  if(fp==NULL)
    {
      if(errno==ENOENT)
        return;
      else
        error(EXIT_FAILURE, errno, "%s: to read as a configuration file",
              filename);
    }


  /* Allocate the space necessary to keep a copy of each line as we parse
     it. Note that `getline' is going to later `realloc' this space to fit
     the line length. */
  errno=0;
  line=malloc(linelen*sizeof *line);
  if(line==NULL)
    error(EXIT_FAILURE, errno, "%zu bytes for `line' in `gal_txt_table_read'",
          linelen*sizeof *line);


  /* Read the parameters line by line. */
  while( getline(&line, &linelen, fp) != -1 )
    {
      ++lineno;
      if( gal_txt_line_stat(line) == GAL_TXT_LINESTAT_DATAROW )
        {
          /* Get the option name and argument/value. */
          options_read_name_arg(line, filename, lineno, &name, &arg);

          /* First look into this program's options, if the option isn't
             found there, `options_set_from_name' will return 1. So the
             condition will succeed and we will start looking into the
             common options, if it isn't found there either, then report an
             error.*/
          if( options_set_from_name(name, arg, poptions, filename, lineno) )
            if( options_set_from_name(name, arg, coptions, filename, lineno) )
              error_at_line(EXIT_FAILURE, 0, filename, lineno,
                            "unrecognized option `%s', for the full list of "
                            "options, please run with `--help'", name);
        }
    }


  /* Close the file. */
  errno=0;
  if(fclose(fp))
    error(EXIT_FAILURE, errno, "%s: couldn't close after reading as "
          "a configuration file", filename);


  /* Clean up and return. */
  free(line);
}




/* Macro to easily specify if we should continue with the rest of the
   configuration files. */
#define OPTIONS_THIS_IS_LASTCONFIG coptions[last_config_index].value     \
  && *((unsigned char *)(coptions[last_config_index].value))





/* Read the configuration files and put the values of the options not given
   into it. The directories containing the configuration files are fixed
   for all the programs.

    - `SYSCONFIG_DIR' is passed onto the library functions at compile time
      from the command-line. You can search for it in the outputs of
      `make'. The main reason is that we want the the user still has the
      chance to change the installation directory after `configure'.

    - `USERCONFIG_DIR' is defined in `config.h'.

    - `CURDIRCONFIG_DIR' is defined in `config.h'. */
static void
gal_options_parse_config_files(char *progexec, struct argp_option *poptions,
                               struct argp_option *coptions)
{
  char *filename;
  size_t i, last_config_index;
  struct gal_linkedlist_stll *tmp;



  /* A small sanity check because in multiple places, we have assumed the
     on/off options have a type of `unsigned char'. */
  if(GAL_OPTIONS_NO_ARG_TYPE != GAL_DATA_TYPE_UCHAR)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix the "
          "problem. The `GAL_OPTIONS_NO_ARG_TYPE' must have the "
          "`unsigned char' type. But",
          PACKAGE_BUGREPORT);



  /* Set the index of the `lastconfig' option. It is important to check
     this option after every configuration file. Thus, to avoid having to
     loop through all the options on every check, we'll find the index and
     just check it each time.

     Although very unlikely, it might happen that the user calls this on
     the command-line (to avoid any configuration files). So, also check if
     we should continue with this function or not. */
  for(i=0; gal_options_is_last(&coptions[i])==0; ++i)
    if(!strcmp(coptions[i].name, "lastconfig"))
      {
        last_config_index=i;
        if(OPTIONS_THIS_IS_LASTCONFIG) return;
        break;
      }



  /* Before doing anything, see if the user specified any configuration
     files to parse before the standard ones. */
  for(i=0; gal_options_is_last(&coptions[i])==0; ++i)
    if(!strcmp(coptions[i].name, "config"))
      {
        /* When more than one configuration file is given, reverse the list
           so it comes into the same order that the user specified on the
           command-line. */
        gal_linkedlist_reverse_stll(
                (struct gal_linkedlist_stll **)(&coptions[i].value) );

        /* Go through all the configuration files and fill in the values. */
        for(tmp=coptions[i].value; tmp!=NULL; tmp=tmp->next)
          {
            options_parse_file(tmp->v, poptions, coptions);
            if( OPTIONS_THIS_IS_LASTCONFIG ) return;
          }
        break;
      }



  /* The program's current directory configuration file. */
  asprintf(&filename, ".%s/%s.conf", PACKAGE, progexec);
  options_parse_file(filename, poptions, coptions);
  if( OPTIONS_THIS_IS_LASTCONFIG ) return;
  free(filename);

  /* General Gnuastro configuration file. */
  asprintf(&filename, ".%s/%s.conf", PACKAGE, PACKAGE);
  printf("%s\n", filename);
  options_parse_file(filename, poptions, coptions);
  if( OPTIONS_THIS_IS_LASTCONFIG ) return;
  free(filename);

  /* User level configuration files. */
  asprintf(&filename, "%s/%s.conf", USERCONFIG_DIR, progexec);
  printf("%s\n", filename);
  options_parse_file(filename, poptions, coptions);
  if( OPTIONS_THIS_IS_LASTCONFIG ) return;
  free(filename);

  /* User level general Gnuastro configuration file. */
  asprintf(&filename, "%s/%s.conf", USERCONFIG_DIR, PACKAGE);
  printf("%s\n", filename);
  options_parse_file(filename, poptions, coptions);
  if( OPTIONS_THIS_IS_LASTCONFIG ) return;
  free(filename);

  /* User level configuration files. */
  asprintf(&filename, "%s/%s.conf", SYSCONFIG_DIR, progexec);
  printf("%s\n", filename);
  options_parse_file(filename, poptions, coptions);
  if( OPTIONS_THIS_IS_LASTCONFIG ) return;
  free(filename);

  /* User level general Gnuastro configuration file. */
  asprintf(&filename, "%s/%s.conf", SYSCONFIG_DIR, PACKAGE);
  printf("%s\n", filename);
  options_parse_file(filename, poptions, coptions);
  if( OPTIONS_THIS_IS_LASTCONFIG ) return;
  free(filename);
}





void
gal_options_config_files(char *progexec, struct argp_option *poptions,
                         struct argp_option *coptions)
{
  /* Parse all the configuration files. */
  gal_options_parse_config_files(progexec, poptions, coptions);

  /* Do the necessary checks (printing, saving and etc). */
}
