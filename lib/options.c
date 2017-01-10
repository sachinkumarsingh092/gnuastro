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

#include <gnuastro/data.h>

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




















/**********************************************************************/
/************            Command-line options           ***************/
/**********************************************************************/
/* Set the value given to the command-line, where we have the integer `key'
   of the option, not its long name as a string. */
error_t
gal_options_set_from_key(int key, char *arg, struct argp_option *options)
{
  size_t i;

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
            gal_data_string_to_type(&options[i].value, arg,
                                    options[i].type);
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
/*
error_t
gal_options_set_from_name(char *name, char *arg, struct argp_option *options)
{

}
*/




/* Read the configuration files and put the values of the options not given
   into it. The directories containing the configuration files are fixed
   for all the programs.

    - `SYSCONFIG_DIR' is passed onto each program at compile time from the
      command-line. You can search for it in the outputs of `make'.

    - `'
*/
void
gal_options_parse_configs(char *progname, struct argp_option *progopts,
                          struct argp_option *commopts)
{
  printf("\nhere\n");
}
