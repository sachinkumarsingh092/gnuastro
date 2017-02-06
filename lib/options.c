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

#include <time.h>
#include <argp.h>
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <string.h>

#include <gnuastro/git.h>
#include <gnuastro/txt.h>
#include <gnuastro/data.h>
#include <gnuastro/table.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/linkedlist.h>

#include <timing.h>
#include <options.h>
#include <checkset.h>









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





void
gal_options_add_to_not_given(struct gal_options_common_params *cp,
                             struct argp_option *option)
{
  gal_linkedlist_add_to_stll(&cp->novalue_doc, (char *)option->doc, 0);
  gal_linkedlist_add_to_stll(&cp->novalue_name, (char *)option->name, 0);
}





void
gal_options_abort_if_mandatory_missing(struct gal_options_common_params *cp)
{
  int namewidth=0;
  char info[5000], *name, *doc;
  struct gal_linkedlist_stll *tmp;

  /* If there is no mandatory options, then just return. */
  if(cp->novalue_name==NULL)
    return;

  /* Get the maximum width of the given names: */
  for(tmp=cp->novalue_name; tmp!=NULL; tmp=tmp->next)
    if( strlen(tmp->v) > namewidth ) namewidth=strlen(tmp->v);

  /* Print the introductory information. */
  sprintf(info, "to continue, the following options need a value ");
  sprintf(info+strlen(info), "(parenthesis after option name contain its "
          "description):\n\n");

  /* Print the list of options along with their description. */
  while(cp->novalue_name!=NULL)
    {
      gal_linkedlist_pop_from_stll(&cp->novalue_doc, &doc);
      gal_linkedlist_pop_from_stll(&cp->novalue_name, &name);
      sprintf(info+strlen(info), "  %-*s (%s\b)\n", namewidth+4, name, doc);
    }
  sprintf(info+strlen(info), "\n");

  /* Print suggestions, way to solve it. */
  sprintf(info+strlen(info), "Use the command-line or a configuration file "
          "to set value(s).\n\nFor a complete description of command-line "
          "options and configuration files, please see the \"Options\" and "
          "\"Configuration files\" section of the Gnuastro book "
          "respectively. You can read them on the command-line by running "
          "the following commands (type `SPACE' to flip through pages, type "
          "`Q' to return to the command-line):\n\n"
          "  info gnuastro Options\n"
          "  info gnuastro \"Configuration files\"\n");

  error(EXIT_FAILURE, 0, "%s", info);
}





static char *
options_get_home()
{
  char *home;
  home=getenv("HOME");
  if(home==NULL)
    error(EXIT_FAILURE, 0, "the HOME environment variable "
          "is not defined");
  return home;
}





















/**********************************************************************/
/************              Option actions               ***************/
/**********************************************************************/

static void
options_check_version(char *version_string)
{
  if( strcmp(version_string, PACKAGE_VERSION) )
    error(EXIT_FAILURE, 0, "version mis-match: you are running GNU "
          "Astronomy Utilities version %s, but this program was configured "
          "to run with version %s (value to the `onlyversion' option, "
          "either in a configuration file or on the command-line). This was "
          "probably done for reproducibility. Therefore, removing, or "
          "changing, the option value might produce errors or unexpected "
          "results. It is hence strongly advised to build GNU Astronomy "
          "Utilities version %s and re-run this command/script",
          PACKAGE_VERSION, version_string, version_string);
}





static void
options_print_citation_exit(struct gal_options_common_params *cp)
{
  char *gnuastro_bibtex=
    "Gnuastro package/infrastructure\n"
    "-------------------------------\n"
    "  @ARTICLE{2015ApJS..220....1A,\n"
    "     author = {{Akhlaghi}, M. and {Ichikawa}, T.},\n"
    "      title = \"{Noise-based Detection and Segmentation of Nebulous "
    "Objects}\",\n"
    "    journal = {\\apjs},\n"
    "  archivePrefix = \"arXiv\",\n"
    "     eprint = {1505.01664},\n"
    "   primaryClass = \"astro-ph.IM\",\n"
    "   keywords = {galaxies: irregular, galaxies: photometry, "
    "galaxies: structure, methods: data analysis, "
    "techniques: image processing, techniques: photometric},\n"
    "       year = 2015,\n"
    "      month = sep,\n"
    "     volume = 220,\n"
    "        eid = {1},\n"
    "      pages = {1},\n"
    "        doi = {10.1088/0067-0049/220/1/1},\n"
    "     adsurl = {http://adsabs.harvard.edu/abs/2015ApJS..220....1A},\n"
    "    adsnote = {Provided by the SAO/NASA Astrophysics Data System}\n"
    "  }";


  /* Print the statements. */
  printf("\nThank you for using %s (%s) %s.\n\n", cp->program_name,
         PACKAGE_NAME, PACKAGE_VERSION);
  printf("Citations are vital for the continued work on Gnuastro.\n\n"
         "Please cite these BibTeX record(s) in your paper(s).\n"
         "(don't forget to also include the version as shown above)\n\n"
         "%s\n\n",
         gnuastro_bibtex);


  /* Only print the citation for the program if one exists. */
  if(cp->program_bibtex[0]!='\0') printf("%s\n\n", cp->program_bibtex);


  /* Print a thank you message. */
  printf("                                               ,\n"
         "                                              {|'--.\n"
         "                                             {{\\    \\\n"
         "      Many thanks from all                   |/`'--./=.\n"
         "      Gnuastro developers!                   `\\.---' `\\\\\n"
         "                                                  |\\  ||\n"
         "                                                  | |//\n"
         "                                                   \\//_/|\n"
         "                                                   //\\__/\n"
         "                                                  //\n"
         "                   (http://www.chris.com/ascii/) |/\n");



  /* Exit the program. */
  exit(EXIT_SUCCESS);
}





/* Declaration of `option_parse_file' that is defined below, since it is
   also needed in `options_immediate', but it fits into context below
   better. */
static void
options_parse_file(char *filename,  struct gal_options_common_params *cp,
                   int enoent_abort);





/* Some options need immediate attention/action before continuing to read
   the rest of the options. In these cases we need to (maybe) check and
   (possibly) abort immediately. */
void
options_immediate(int key, char *arg, struct gal_options_common_params *cp)
{
  switch(key)
    {
    /* We don't want later options that were set for the given version to
       cause errors if this option was given. */
    case GAL_OPTIONS_KEY_ONLYVERSION:
      options_check_version(arg);
      break;

    /* This option is completely independent of anything and doesn't need
       out attention later. */
    case GAL_OPTIONS_KEY_CITE:
      options_print_citation_exit(cp);
      break;

    /* Read the given configuration file. */
    case GAL_OPTIONS_KEY_CONFIG:
      options_parse_file(arg, cp, 1);
      break;
    }
}





/* The option value has been read and put into the `value' field of the
   `argp_option' structure. This function will use the `range' field to
   define a check and abort with an error if the value is not in the given
   range. It also takes the `arg' so it can be used for good error message
   (showing the value that could not be read). */
static void
options_sanity_check(struct argp_option *option, char *arg,
                     char *filename, size_t lineno)
{
  size_t dsize=1;
  char *message=NULL;
  int operator1=GAL_ARITHMETIC_OP_INVALID;
  int operator2=GAL_ARITHMETIC_OP_INVALID;
  int multicheckop=GAL_ARITHMETIC_OP_INVALID;
  gal_data_t *value, *ref1=NULL, *ref2=NULL, *check1, *check2;
  int mcflag = ( GAL_ARITHMETIC_NUMOK
                 | GAL_ARITHMETIC_FREE
                 | GAL_ARITHMETIC_INPLACE );

  /* Currently, this function is only for numeric types, so if the value is
     string type, or its `range' field is `GAL_OPTIONS_RANGE_ANY', then
     just return without any checks. */
  if( option->type==GAL_DATA_TYPE_STRING
      || option->type==GAL_DATA_TYPE_STRLL
      || option->range==GAL_OPTIONS_RANGE_ANY )
    return;

  /* Put the option value into a data structure. */
  value=gal_data_alloc(option->value, option->type, 1, &dsize, NULL,
                       0, -1, NULL, NULL, NULL);

  /* Set the operator(s) and operands: */
  switch(option->range)
    {

    case GAL_OPTIONS_RANGE_GT_0:
      message="greater than zero";
      ref1=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      *(unsigned char *)(ref1->array)=0;
      operator1=GAL_ARITHMETIC_OP_GT;
      ref2=NULL;
      break;


    case GAL_OPTIONS_RANGE_GE_0:
      message="greater or equal to zero";
      ref1=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      *(unsigned char *)(ref1->array)=0;
      operator1=GAL_ARITHMETIC_OP_GE;
      ref2=NULL;
      break;


    case GAL_OPTIONS_RANGE_0_OR_1:
      message="either 0 or 1";
      ref1=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      ref2=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      *(unsigned char *)(ref1->array)=0;
      *(unsigned char *)(ref2->array)=1;

      operator1=GAL_ARITHMETIC_OP_EQ;
      operator2=GAL_ARITHMETIC_OP_EQ;
      multicheckop=GAL_ARITHMETIC_OP_OR;
      break;


    case GAL_OPTIONS_RANGE_GE_0_LE_1:
      message="between zero and one";
      ref1=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      ref2=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      *(unsigned char *)(ref1->array)=0;
      *(unsigned char *)(ref2->array)=1;

      operator1=GAL_ARITHMETIC_OP_GE;
      operator2=GAL_ARITHMETIC_OP_LE;
      multicheckop=GAL_ARITHMETIC_OP_AND;
      break;


    case GAL_OPTIONS_RANGE_GT_0_ODD:
      message="greater than zero and odd";
      ref1=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      ref2=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      *(unsigned char *)(ref1->array)=0;
      *(unsigned char *)(ref2->array)=2;

      operator1=GAL_ARITHMETIC_OP_GT;
      operator2=GAL_ARITHMETIC_OP_MODULO;
      multicheckop=GAL_ARITHMETIC_OP_AND;
      break;

    case GAL_OPTIONS_RANGE_0_OR_ODD:
      message="greater than, or equal to, zero and odd";
      ref1=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      ref2=gal_data_alloc(NULL, GAL_DATA_TYPE_UCHAR, 1, &dsize, NULL,
                          0, -1, NULL, NULL, NULL);
      *(unsigned char *)(ref1->array)=0;
      *(unsigned char *)(ref2->array)=2;

      operator1=GAL_ARITHMETIC_OP_EQ;
      operator2=GAL_ARITHMETIC_OP_MODULO;
      multicheckop=GAL_ARITHMETIC_OP_OR;
      break;

    default:
      error(EXIT_FAILURE, 0, "range code %d not recognized in "
            "`gal_options_check_set'", option->range);
    }


  /* Use the arithmetic library to check for the condition. We don't want
     to free the value or change its value, so when dealing with the value
     directly, we won't use the `GAL_ARITHMETIC_FREE', or
     `GAL_ARITHMETIC_INPLACE' flags. But we will do this when there are
     multiple checks so from the two check data structures, we only have
     one remaining. */
  check1=gal_arithmetic(operator1, GAL_ARITHMETIC_NUMOK, value, ref1);
  if(ref2)
    {
      check2=gal_arithmetic(operator2, GAL_ARITHMETIC_NUMOK, value, ref2);
      check1=gal_arithmetic(multicheckop, mcflag, check1, check2);
    }


  /* If the final check is not successful, then print an error. */
  if( *(unsigned char *)(check1->array)==0 )
    error_at_line(EXIT_FAILURE, 0, filename, lineno,
                  "value to option `%s' must be %s, but the given value "
                  "is `%s'. Recall that `%s' is \"%s\"", option->name,
                  message, arg, option->name, option->doc);


  /* Clean up and finish. Note that we used the actual value pointer in the
     data structure, so first we need to set it to NULL, so `gal_data_free'
     doesn't free it, we need it for later (for example to print the option
     values). */
  value->array=NULL;
  gal_data_free(ref1);
  gal_data_free(ref2);
  gal_data_free(value);
  gal_data_free(check1);
}





static void
gal_options_read_check(struct argp_option *option, char *arg, char *filename,
                       size_t lineno)
{
  if(option->type==GAL_DATA_TYPE_STRLL)
    gal_linkedlist_add_to_stll(option->value, arg, 1);
  else
    {
      /* Read the string argument into the value. */
      if( gal_data_string_to_type(&option->value, arg, option->type) )
        /* Fortunately `error_at_line' will behave like `error' when the
           filename is NULL (the option was read from a command-line). */
        error_at_line(EXIT_FAILURE, 0, filename, lineno,
                      "`%s' (value to option `--%s') couldn't be read into "
                      "the proper numerical type. Common causes for this "
                      "error are:\n"
                      "  - It contains non-numerical characters.\n"
                      "  - It is negative, but the expected value is "
                      "positive.\n"
                      "  - It is floating point, but the expected value "
                      "is an integer.\n"
                      "  - The previous option required a value, but you "
                      "forgot to give it one, so the next option's "
                      "name(+value, if there are no spaces between them) "
                      "is read as the value of the previous option.", arg,
                      option->name);
    }

  /* Do a sanity check. */
  options_sanity_check(option, arg, filename, lineno);


  /* Flip the `set' flag to `GAL_OPTIONS_SET'. */
  option->set=GAL_OPTIONS_SET;
}



















/**********************************************************************/
/************            Command-line options           ***************/
/**********************************************************************/
/* Set the value given to the command-line, where we have the integer `key'
   of the option, not its long name as a string. */
error_t
gal_options_set_from_key(int key, char *arg, struct argp_option *options,
                         struct gal_options_common_params *cp)
{
  size_t i;

  /* Go through all the options and find the one that should keep this
     value, then put its value into the appropriate key. Note that the
     options array finishs with an all zero element, so we don't need to
     know the number before hand.*/
  for(i=0;1;++i)
    {
      /* Check if the key corresponds to this option. */
      if( options[i].key==key )
        {
          /* For options that need immediate attention. */
          options_immediate(key, arg, cp);

          /* When options are read from keys (by this function), they are
             read from the command-line. On the commandline, the last
             invokation of the option is important. Especially in contexts
             like scripts, this is important because you can change a given
             command-line option (that is not a linked list) by calling it
             a second time, instead of going back and changing the first
             value.

             As a result, only when searching for options on the
             command-line, a second value to the same option will replace
             the first one. This will not happen in configuration files. */
          if(options[i].set && gal_data_is_linked_list(options[i].type)==0)
            options[i].set=GAL_OPTIONS_NOT_SET;

          /* We have two types of options: those which need an argument and
             those that don't. For those that don't `arg' will be
             NULL. When it accepts an argument then read itinto the option
             structure and do a sanity check.*/
          if(arg)
            gal_options_read_check(&options[i], arg, NULL, 0);
          else
            {
              /* Make sure the option has the type set for options with no
                 argument. So, give it a value of 1. */
              if(options[i].type==GAL_OPTIONS_NO_ARG_TYPE)
                {
                  *(unsigned char *)(options[i].value)=1;
                  options[i].set=GAL_OPTIONS_SET;
                }
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
  struct gal_options_common_params *cp=state->input;

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

  /* Read the options. */
  return gal_options_set_from_key(key, arg, cp->coptions, cp);
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
options_set_from_name(char *name, char *arg,  struct argp_option *options,
                      struct gal_options_common_params *cp, char *filename,
                      size_t lineno)
{
  size_t i;

  /* Go through all the options and find the one that should keep this
     value, then put its value into the appropriate key. Note that the
     options array finishs with an all zero element, so we don't need to
     know the number before hand.*/
  for(i=0;1;++i)
    {
      /* Check if the key corresponds to this option. */
      if( options[i].name && !strcmp(options[i].name, name) )
        {
          /* If the option already has a value and it isn't a linked list,
             or it is not relevant to this program, then ignore it. */
          if( options[i].set && !gal_data_is_linked_list(options[i].type ) )
            return 0;

          /* For options that need immediate attention. */
          options_immediate(options[i].key, arg, cp);

          /* Read the value into the option and do a sanity check. */
          gal_options_read_check(&options[i], arg, filename, lineno);

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





/* If the last config option has a value which is 1, then in some previous
   configuration file the user has asked to stop parsing configuration
   files. In that case, don't read this configuration file. */
static int
options_lastconfig_has_been_called(struct argp_option *coptions)
{
  size_t i;

  for(i=0; !gal_options_is_last(&coptions[i]); ++i)
    if( coptions[i].key == GAL_OPTIONS_KEY_LASTCONFIG
        && coptions[i].set
        && *((unsigned char *)(coptions[i].value)) )
      return 1;
  return 0;
}





static void
options_parse_file(char *filename,  struct gal_options_common_params *cp,
                   int enoent_abort)
{
  FILE *fp;
  char *line, *name, *arg;
  size_t linelen=10, lineno=0;


  /* If `lastconfig' was called prior to this file, then just return and
     ignore this configuration file. */
  if( options_lastconfig_has_been_called(cp->coptions) )
    return;


  /* Open the file. If the file doesn't exist, then just ignore the
     configuration file and return. */
  errno=0;
  fp=fopen(filename, "r");
  if(fp==NULL)
    {
      if(errno==ENOENT && enoent_abort==0)
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
          if( options_set_from_name(name, arg, cp->poptions, cp,
                                    filename, lineno) )
            if( options_set_from_name(name, arg, cp->coptions, cp,
                                      filename, lineno) )
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
gal_options_parse_config_files(struct gal_options_common_params *cp)
{
  char *home;
  char *filename;

  /* A small sanity check because in multiple places, we have assumed the
     on/off options have a type of `unsigned char'. */
  if(GAL_OPTIONS_NO_ARG_TYPE != GAL_DATA_TYPE_UCHAR)
    error(EXIT_FAILURE, 0, "A bug! Please contact us at %s so we can fix the "
          "problem. The `GAL_OPTIONS_NO_ARG_TYPE' must have the "
          "`unsigned char' type. But",
          PACKAGE_BUGREPORT);

  /* Common options configuration file. */
  asprintf(&filename, ".%s/%s.conf", PACKAGE, PACKAGE);
  options_parse_file(filename, cp, 0);
  free(filename);

  /* The program's current directory configuration file. */
  asprintf(&filename, ".%s/%s.conf", PACKAGE, cp->program_exec);
  options_parse_file(filename, cp, 0);
  free(filename);

  /* Read the home environment variable. */
  home=options_get_home();

  /* Common options user-wide configuration file. */
  asprintf(&filename, "%s/%s/%s.conf", home, USERCONFIG_DIR, PACKAGE);
  options_parse_file(filename, cp, 0);
  free(filename);

  /* The program's user-wide configuration file. */
  asprintf(&filename, "%s/%s/%s.conf", home, USERCONFIG_DIR,
           cp->program_exec);
  options_parse_file(filename, cp, 0);
  free(filename);

  /* Common options system-wide configuration file. */
  asprintf(&filename, "%s/%s.conf", SYSCONFIG_DIR, PACKAGE);
  options_parse_file(filename, cp, 0);
  free(filename);

  /* The program's system-wide configuration file. */
  asprintf(&filename, "%s/%s.conf", SYSCONFIG_DIR, cp->program_exec);
  options_parse_file(filename, cp, 0);
  free(filename);
}





static void
options_reverse_lists_check_mandatory(struct gal_options_common_params *cp,
                                      struct argp_option *options)
{
  size_t i;

  for(i=0; !gal_options_is_last(&options[i]); ++i)
    {
      if(options[i].set)
        {
          if( gal_data_is_linked_list(options[i].type) )
            switch(options[i].type)
              {
              case GAL_DATA_TYPE_STRLL:
                gal_linkedlist_reverse_stll(
                    (struct gal_linkedlist_stll **)(options[i].value) );
                break;
              }
        }
      else if(options[i].mandatory==GAL_OPTIONS_MANDATORY)
        gal_options_add_to_not_given(cp, &options[i]);
    }
}





/* Read all configuration files and set common options */
void
gal_options_read_config_set(struct gal_options_common_params *cp)
{
  size_t i;
  struct argp_option *coptions=cp->coptions;

  /* Parse all the configuration files. */
  gal_options_parse_config_files(cp);

  /* Reverse the order of all linked list type options so the popping order
     is the same as the user's input order. We need to do this here because
     when printing those options, their order matters.*/
  options_reverse_lists_check_mandatory(cp, cp->poptions);
  options_reverse_lists_check_mandatory(cp, cp->coptions);

  /* Some of the values that the user gives as strings should be stored
     internally as integers. For the program-specific options, this is done
     in the program's `ui.c', here, we do it for the common options. */
  for(i=0; !gal_options_is_last(&coptions[i]); ++i)
    if(coptions[i].set)
      switch(coptions[i].key)
      {
      case GAL_OPTIONS_KEY_SEARCHIN:
        cp->searchin=gal_table_string_to_searchin(cp->searchinstr);
        break;
      case GAL_OPTIONS_KEY_TABLEFORMAT:
        cp->tableformat=gal_table_string_to_format(cp->tableformatstr);
        break;
      }

  /* Abort if any of the mandatory options are not set. */
  gal_options_abort_if_mandatory_missing(cp);
}




















/**********************************************************************/
/************              Printing/Writing             ***************/
/**********************************************************************/
/* We don't want to print the values of configuration specific options and
   the output option. The output value is assumed to be specific to each
   input, and the configuration options are for reading the configuration,
   not writing it. */
static int
option_is_printable(struct argp_option *option)
{
  /* First check if option is hidden (not relevant to this program). */
  if(option->flags & OPTION_HIDDEN)
    return 0;

  /* Then check if it is a pre-program option. */
  switch(option->key)
    {
    case GAL_OPTIONS_KEY_OUTPUT:
    case GAL_OPTIONS_KEY_CITE:
    case GAL_OPTIONS_KEY_PRINTPARAMS:
    case GAL_OPTIONS_KEY_CONFIG:
    case GAL_OPTIONS_KEY_SETDIRCONF:
    case GAL_OPTIONS_KEY_SETUSRCONF:
    case GAL_OPTIONS_KEY_LASTCONFIG:
      return 0;
    }
  return 1;
}





/* For a given type, print the value in `ptr' in a space of `width'
   elements. If `width==0', then return the width necessary to print the
   value. */
static int
options_print_any_type(void *ptr, int type, int width, FILE *fp)
{
  char *str;

  /* Write the value into a string. */
  str=gal_data_write_to_string(ptr, type, 1);

  /* If only the width was desired, don't actually print the string, just
     return its length. Otherwise, print it. */
  if(width)
    fprintf(fp, "%-*s ", width, str);
  else
    width=strlen(str);

  /* Free the allocated space and return. */
  free(str);
  return width;
}





/* An option structure is given, return its name and value print
   lengths. */
static void
options_correct_max_lengths(struct argp_option *option, int *max_nlen,
                            int *max_vlen)
{
  int vlen;
  struct gal_linkedlist_stll *tmp;

  /* Get the length of the value and save its length length if its
     larger than the widest value. */
  if(gal_data_is_linked_list(option->type))
    {
      /* A small sanity check. */
      if(option->type!=GAL_DATA_TYPE_STRLL)
        error(EXIT_FAILURE, 0, "currently only string linked lists "
              "are acceptable for printing");

      /* Check each node, one by one. */
      for(tmp=*(struct gal_linkedlist_stll **)(option->value);
          tmp!=NULL; tmp=tmp->next)
        {
          /* Get the length of this node: */
          vlen=options_print_any_type(&tmp->v, GAL_DATA_TYPE_STRING, 0, NULL);

          /* If its larger than the maximum length, then put it in. */
          if( vlen > *max_vlen )
            *max_vlen=vlen;
        }
    }
  else
    {
      vlen=options_print_any_type(option->value, option->type, 0, NULL);
      if( vlen > *max_vlen )
        *max_vlen=vlen;
    }

  /* If the name of this option is larger than all existing, set its
     length as the largest name length. */
  if( strlen(option->name) > *max_nlen )
    *max_nlen = strlen(option->name);
}





/* To print the options nicely, we need the maximum lengths of the options
   and their values. */
static void
options_set_lengths(struct argp_option *poptions,
                    struct argp_option *coptions, int *namelen, int *valuelen)
{
  int i, max_nlen=0, max_vlen=0;

  /* For program specific options. */
  for(i=0; !gal_options_is_last(&poptions[i]); ++i)
    if(poptions[i].name && poptions[i].set)
      options_correct_max_lengths(&poptions[i], &max_nlen, &max_vlen);

  /* For common options. Note that the options that will not be printed are
     in this category, so we also need to check them. The detailed steps
     are the same as before. */
  for(i=0; !gal_options_is_last(&coptions[i]); ++i)
    if( coptions[i].name && coptions[i].set
        && option_is_printable(&coptions[i]) )
      options_correct_max_lengths(&coptions[i], &max_nlen, &max_vlen);

  /* Save the final values in the output pointers. */
  *namelen=max_nlen;
  *valuelen=max_vlen;
}





/* The `#' before the `doc' string are not required by the configuration
   file parser when the documentation string fits in a line. However, when
   the `doc' string is longer than 80 characters, it will be cut between
   multiple lines and without the `#', the start of the line will be read
   as an option. */
static void
options_print_doc(FILE *fp, const char *doc, int nvwidth)
{
  size_t len=strlen(doc);

  /* The `+3' is because of the three extra spaces in this line: one before
     the variable name, one after it and one after the value. */
  int i, prewidth=nvwidth+3, width=78-prewidth, cwidth;

  /* We only want the formatting when writing to stdout. */
  if(len<width)
    fprintf(fp, "# %s\n", doc);
  else
    {
      /* If the break is in the middle of a word, then pull set it before
         the word starts.*/
      cwidth=width; while( doc[cwidth]!=' ' ) --cwidth;
      fprintf(fp, "# %.*s\n", cwidth, doc);
      i=cwidth;

      /* Go over the rest of the line */
      while(i<len)
        {
          /* Remove any possible space before the first word. */
          while( doc[i]==' ' ) ++i;

          /* Check if the line break won't fall in the middle of a word. */
          cwidth=width;
          if( i+cwidth<len) while( doc[i+cwidth]!=' ' ) --cwidth;
          fprintf(fp, "%*s# %.*s\n", prewidth, "", cwidth, &doc[i]);
          i+=cwidth;
        }
    }
}





static void
options_print_all_in_group(struct argp_option *options, int groupint,
                           int namelen, int valuelen, FILE *fp)
{
  size_t i;
  struct gal_linkedlist_stll *tmp;
  int namewidth=namelen+1, valuewidth=valuelen+1;

  /* Go over all the options. */
  for(i=0; !gal_options_is_last(&options[i]); ++i)
    if( options[i].group == groupint          /* Is in this group.        */
        && options[i].set                     /* Has been given a value.  */
        && option_is_printable(&options[i]) ) /* Is relevant for printing.*/
      {
        /* Linked lists */
        if(gal_data_is_linked_list(options[i].type))
          for(tmp=*(struct gal_linkedlist_stll **)(options[i].value);
              tmp!=NULL; tmp=tmp->next)
            {
              fprintf(fp, " %-*s ", namewidth, options[i].name);
              options_print_any_type(&tmp->v, GAL_DATA_TYPE_STRING,
                                     valuewidth, fp);
              options_print_doc(fp, options[i].doc, namewidth+valuewidth);
            }

        /* Normal types. */
        else
          {
            fprintf(fp, " %-*s ", namewidth, options[i].name);
            options_print_any_type(options[i].value, options[i].type,
                                   valuewidth, fp);
            options_print_doc(fp, options[i].doc, namewidth+valuewidth);
          }
      }
}





static void
options_print_all(struct gal_options_common_params *cp, char *dirname,
                  const char *optionname)
{
  size_t i;
  FILE *fp;
  time_t rawtime;
  char *topicstr, *filename;
  int groupint, namelen, valuelen;
  struct gal_linkedlist_ill *group=NULL;
  struct gal_linkedlist_stll *topic=NULL;
  struct argp_option *coptions=cp->coptions, *poptions=cp->poptions;

  /* If the configurations are to be written to a file, then do the
     preparations. */
  if(dirname)
    {
      /* Make the host directory if it doesn't already exist. */
      gal_checkset_mkdir(dirname);

      /* Prepare the full filename: */
      asprintf(&filename, "%s/%s.conf", dirname, cp->program_exec);

      /* Remove the file if it already exists. */
      gal_checkset_check_remove_file(filename, 0);

      /* Open the file for writing */
      errno=0;
      fp=fopen(filename, "w");
      if(fp==NULL)
        error(EXIT_FAILURE, errno, "%s: could't open to write "
              "configuration file", dirname);

      /* Print the basic information as comments in the file first. */
      time(&rawtime);
      fprintf(fp,
              "# %s (%s) %s.\n"
              "# Written at %s#\n"
              "#  - Empty lines are ignored.\n"
              "#  - Lines starting with `#` are ignored.\n"
              "#  - The long option name is followed by a value.\n"
              "#  - The name and value should be separated by atleast\n"
              "#    one white space character (for example space or tab).\n"
              "#  - If the value has space, enclose the whole value in\n"
              "#    double quotation (\") signs.\n"
              "#  - After the value, the rest of the line is ignored.\n"
              "#\n# Run `info %s' for a more elaborate description of each "
              "option.\n",
              cp->program_name, PACKAGE_NAME, PACKAGE_VERSION,
              ctime(&rawtime), cp->program_exec);
    }
  else fp=stdout;

  /* Parse all the options with a title, note that the `Input', `Output'
     and `Operating mode' options are defined in the common options, while
     the (possible) other groups are in the program specific options. We
     will only be dealing with the `topics' linked list in this function
     and the strings in `poption' are statically allocated, so its fine to
     not waste CPU cycles allocating and freeing.*/
  for(i=0; !gal_options_is_last(&coptions[i]); ++i)
    if(coptions[i].name==NULL && coptions[i].key==0 && coptions[i].doc)
      {
        /* The `(char *)' is because `.doc' is a constant and this helps
           remove the compiler warning. */
        gal_linkedlist_add_to_ill(&group, coptions[i].group);
        gal_linkedlist_add_to_stll(&topic, (char *)coptions[i].doc, 0);
      }
  for(i=0; !gal_options_is_last(&poptions[i]); ++i)
    if(poptions[i].name==NULL && poptions[i].key==0 && poptions[i].doc)
      {
        gal_linkedlist_add_to_ill(&group, poptions[i].group);
        gal_linkedlist_add_to_stll(&topic, (char *)poptions[i].doc, 0);
      }

  /* Reverse the linked lists to get the same input order. */
  gal_linkedlist_reverse_stll(&topic);
  gal_linkedlist_reverse_ill(&group);

  /* Get the maximum width of names and values. */
  options_set_lengths(poptions, coptions, &namelen, &valuelen);

  /* Go over each topic and print every option that is in this group. */
  while(topic)
    {
      /* Pop the nodes from the linked list. */
      gal_linkedlist_pop_from_ill(&group, &groupint);
      gal_linkedlist_pop_from_stll(&topic, &topicstr);

      /* First print the topic, */
      fprintf(fp, "\n# %s\n", topicstr);
      fprintf(fp, "# ");
      i=0; while(i++<strlen(topicstr)) fprintf(fp, "%c", '-');
      fprintf(fp, "\n");

      /* Then, print all the options that are in this group. */
      options_print_all_in_group(coptions, groupint, namelen, valuelen, fp);
      options_print_all_in_group(poptions, groupint, namelen, valuelen, fp);
    }

  /* Let the user know. */
  if(dirname)
    {
      printf("\nNew/updated configuration file:\n\n  %s\n\n"
             "You may inspect it with `cat %s'.\n"
             "You may use your favorite text editor to modify it later.\n"
             "Or, run %s again with new values for the options and `--%s'.\n",
             filename, filename, cp->program_name, optionname);
      free(filename);
    }

  /* Exit the program successfully */
  exit(EXIT_SUCCESS);
}





#define OPTIONS_UCHARVAL *(unsigned char *)(cp->coptions[i].value)
void
gal_options_print_state(struct gal_options_common_params *cp)
{
  size_t i;
  unsigned char sum=0;
  char *home, *dirname;


  /* A sanity check is necessary first. We want to make sure that the user
     hasn't called more than one of these options. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    if(cp->coptions[i].set)
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_PRINTPARAMS:
        case GAL_OPTIONS_KEY_SETDIRCONF:
        case GAL_OPTIONS_KEY_SETUSRCONF:
          sum += OPTIONS_UCHARVAL;
        }
  if(sum>1)
    error(EXIT_FAILURE, 0, "only one of the `printparams', `setdirconf' "
          "and `setusrconf' options can be called in each run");


  /* Print the required configuration files. Note that simply having a
     non-NULL value is not enough. They can have a value of 1 or 0, and the
     respective file should only be created if we have a value of 1. */
  for(i=0; !gal_options_is_last(&cp->coptions[i]); ++i)
    if(cp->coptions[i].set && OPTIONS_UCHARVAL)
      switch(cp->coptions[i].key)
        {
        case GAL_OPTIONS_KEY_PRINTPARAMS:
          options_print_all(cp, NULL, NULL);
          break;

        case GAL_OPTIONS_KEY_SETDIRCONF:
          asprintf(&dirname, ".%s", PACKAGE);
          options_print_all(cp, dirname, cp->coptions[i].name);
          free(dirname);
          break;

        case GAL_OPTIONS_KEY_SETUSRCONF:
          home=options_get_home();
          asprintf(&dirname, "%s/%s", home, USERCONFIG_DIR);
          options_print_all(cp, dirname, cp->coptions[i].name);
          free(dirname);
          break;
        }
}





void
gal_options_print_log(gal_data_t *logll, char *program_string,
                      time_t *rawtime, char *comments, char *filename,
                      struct gal_options_common_params *cp)
{
  char *finalcomment;
  char *msg, gitdescribe[100], *gd;

  /* Get the Git description in the running folder. */
  gd=gal_git_describe();
  if(gd) sprintf(gitdescribe, " from %s,", gd);
  else   gitdescribe[0]='\0';

  /* Write the top level comment. */
  asprintf(&finalcomment, "# %s\n# Created%s on %s%s",
           program_string, gitdescribe, ctime(rawtime),
           comments ? comments : "#");

  /* Write the log file to disk */
  gal_table_write(logll, finalcomment, GAL_TABLE_FORMAT_TXT, filename,
                  cp->dontdelete);

  /* In verbose mode, print the information. */
  if(!cp->quiet)
    {
      asprintf(&msg, "%s created.", filename);
      gal_timing_report(NULL, msg, 1);
      free(msg);
    }

  /* Clean up. */
  free(finalcomment);
}
