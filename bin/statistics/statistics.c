/*********************************************************************
Statistics - Statistical analysis on input dataset.
Statistics is part of GNU Astronomy Utilities (Gnuastro) package.

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
#include <config.h>

#include <math.h>
#include <errno.h>
#include <error.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <gnuastro/fits.h>
#include <gnuastro/arithmetic.h>
#include <gnuastro/statistics.h>

#include <timing.h>
#include <checkset.h>

#include "main.h"

#include "ui.h"
#include "statistics.h"





/*******************************************************************/
/**************           Print in one row           ***************/
/*******************************************************************/
static void
ui_print_one_number(gal_data_t *out)
{
  char *toprint=gal_data_write_to_string(out->array, out->type, 0);
  printf("%s ", toprint);
  gal_data_free(out);
  free(toprint);
}





static void
ui_print_one_row(struct statisticsparams *p)
{
  struct gal_linkedlist_ill *tmp;

  /* Print the numbers. */
  for(tmp=p->toprint; tmp!=NULL; tmp=tmp->next)
    switch(tmp->v)
      {
      case ARGS_OPTION_KEY_NUMBER:
        ui_print_one_number( gal_statistics_number(p->input) );      break;
      case ARGS_OPTION_KEY_MINIMUM:
        ui_print_one_number( gal_statistics_minimum(p->input) );     break;
      case ARGS_OPTION_KEY_MAXIMUM:
        ui_print_one_number( gal_statistics_maximum(p->input) );     break;
      case ARGS_OPTION_KEY_SUM:
        ui_print_one_number( gal_statistics_sum(p->input) );         break;
      case ARGS_OPTION_KEY_MEAN:
        ui_print_one_number( gal_statistics_mean(p->input) );        break;
      case ARGS_OPTION_KEY_STD:
        ui_print_one_number( gal_statistics_std(p->input) );         break;
      case ARGS_OPTION_KEY_MEDIAN:
        ui_print_one_number( gal_statistics_median(p->sorted, 0) );  break;
      case ARGS_OPTION_KEY_MODE:
        error(EXIT_FAILURE, 0, "mode isn't implemented yet!");
        break;
      default:
        error(EXIT_FAILURE, 0, "A bug! Operation code %d not recognized in "
              "`ui_print_row'. Please contact us at %s so we can address "
              "the problem", tmp->v, PACKAGE_BUGREPORT);
      }

  /* Print a new line. */
  printf("\n");
}




















/*******************************************************************/
/**************             ASCII plots              ***************/
/*******************************************************************/
static void
print_ascii_plot(struct statisticsparams *p, gal_data_t *plot,
                 gal_data_t *bins, int h1_c0, int printinfo)
{
  int i, j;
  size_t *s, *sf, max=0;
  double *b, v, halfbinwidth, correction;

  /* Find the maximum of the plot. */
  sf=(s=plot->array)+plot->size; do max = *s>max ? *s : max; while(++s<sf);

  /* Print the range so the user knows. */
  if(printinfo)
    {
      b=bins->array;
      halfbinwidth = (b[1]-b[0])/2;
      printf("\nASCII %s:\n", ( h1_c0 ? "Histogram" :
                                "Cumulative frequency plot") );
      if(h1_c0) printf("Number: %zu\n", p->input->size);
      printf("Y: (linear: 0 to %zu)\n", max);
      printf("X: (linear: %g -- %g, in %zu bins)\n", b[0]-halfbinwidth,
             b[ bins->size - 1 ] + halfbinwidth, bins->size);
    }

  /* Print the ASCII plot: */
  s=plot->array;
  correction = (double)(p->asciiheight) / (double)max;
  for(i=p->asciiheight;i>=0;--i)
    {
      printf(" |");
      for(j=0;j<plot->size;++j)
        {
          v = (double)s[j] * correction;
          if( v >= ((double)i-0.5f) && v > 0.0f ) printf("*");
          else printf(" ");
        }
      printf("\n");
    }
  printf(" |");
  for(j=0;j<plot->size;++j) printf("-");
  printf("\n\n");
}




static void
ascii_plots(struct statisticsparams *p)
{
  gal_data_t *bins, *hist, *cfp;

  /* Make the bins and the respective plot. */
  bins=gal_statistics_regular_bins(p->input, NULL, p->numasciibins, NAN);
  hist=gal_statistics_histogram(p->input, bins, 0, 0);
  if(p->asciicfp)
    {
      bins->next=hist;
      cfp=gal_statistics_cfp(p->input, bins, 0);
    }

  /* Print the plots. */
  if(p->asciihist)  print_ascii_plot(p, hist, bins, 1, 1);
  if(p->asciicfp)   print_ascii_plot(p, cfp,  bins, 0, 1);

  /* Clean up.*/
  gal_data_free(bins);
  gal_data_free(hist);
  if(p->asciicfp) gal_data_free(cfp);
}




















/*******************************************************************/
/*******    Histogram and cumulative frequency tables    ***********/
/*******************************************************************/
static void
save_hist_and_or_cfp(struct statisticsparams *p)
{
  gal_data_t *bins, *hist, *cfp=NULL;
  struct gal_linkedlist_stll *comments=NULL;
  char *tmp, *contents, *output, *suf, *fix, *suffix;


  /* Set the bins and make the histogram, this is necessary for both the
     histogram and CFP (recall that the CFP is built from the
     histogram). */
  bins=gal_statistics_regular_bins(p->input, NULL, p->numbins,
                                   p->onebinstart);
  hist=gal_statistics_histogram(p->input, bins, p->normalize, p->maxbinone);


  /* Set the histogram as the next pointer of bins. This is again necessary
     in both cases: when only a histogram is requested, it is used for the
     plotting. When only a CFP is desired, it is used as input into
     `gal_statistics_cfp'. */
  bins->next=hist;


  /* Make the cumulative frequency plot if the user wanted it. Make the
     CFP, note that for the CFP, `maxbinone' and `normalize' are the same:
     the last bin (largest value) must be one. So if any of them are given,
     then set the last argument to 1.*/
  if(p->cumulative)
    cfp=gal_statistics_cfp(p->input, bins, p->normalize || p->maxbinone);


  /* FITS tables don't accept `uint64_t', so to be consistent, we'll conver
     the histogram and CFP to `uint32_t'.*/
  if(hist->type==GAL_DATA_TYPE_UINT64)
    hist=gal_data_copy_to_new_type_free(hist, GAL_DATA_TYPE_UINT32);
  if(cfp && cfp->type==GAL_DATA_TYPE_UINT64)
    cfp=gal_data_copy_to_new_type_free(cfp, GAL_DATA_TYPE_UINT32);


  /* Finalize the next pointers. */
  bins->next=hist;
  hist->next=cfp;


  /* Prepare the contents. */
  if(p->histogram && p->cumulative)
    { suf="_hist_cfp"; contents="Histogram and cumulative frequency plot"; }
  else if(p->histogram)
    { suf="_hist";     contents="Histogram"; }
  else
    { suf="_cfp";      contents="Cumulative frequency plot"; }


  /* Set the basename and output suffix. */
  fix = ( p->cp.output
          ? gal_fits_name_is_fits(p->cp.output) ? "fits" : "txt"
          : "txt" );
  asprintf(&suffix, "%s.%s", suf, fix);


  /* If another file is to be created or output name isn't set, we'll need
     to use Automatic output. */
  output = ( p->cp.output
             ? p->cp.output
             : gal_checkset_automatic_output(&p->cp, p->inputname, suffix) );


  /* Write the comments, NOTE: we are writing the first two in reverse of
     the order we want them. They will later be freed as part of the list's
     freeing.*/
  tmp=gal_fits_name_save_as_string(p->inputname, p->cp.hdu);
  gal_linkedlist_add_to_stll(&comments, tmp, 0);

  asprintf(&tmp, "%s created from:", contents);
  gal_linkedlist_add_to_stll(&comments, tmp, 0);

  if(strcmp(fix, "fits"))  /* The intro info will be in FITS files anyway.*/
    gal_table_comments_add_intro(&comments, PROGRAM_STRING, &p->rawtime);


  /* Write the table. */
  gal_table_write(bins, comments, p->cp.tableformat, output,
                  p->cp.dontdelete);


  /* Let the user know, if we aren't in quiet mode. */
  if(!p->cp.quiet)
    printf("%s created.\n", output);


  /* Clean up. */
  free(suffix);
  if(output!=p->cp.output) free(output);
  gal_linkedlist_free_stll(comments, 1);
}



















/*******************************************************************/
/**************           Basic information          ***************/
/*******************************************************************/
/* To keep things in `print_basics' clean, we'll define the input data
   here, then only print the values there. */
void
print_input_info(struct statisticsparams *p)
{
  char *str, *name, *col=NULL;

  /* Print the program name and version. */
  printf("%s\n", PROGRAM_STRING);

  /* Print the input information, if the input was a table, we also need to
     give the column information. When the column has a name, it will be
     printed, when it doesn't, we'll use the same string the user gave. */
  printf("-------\n");
  name=gal_fits_name_save_as_string(p->inputname, p->cp.hdu);
  printf("Input: %s\n", name);

  /* If a table was given, print the column. */
  if(p->column) printf("Column: %s\n", p->column);

  /* Range. */
  str=NULL;
  if( !isnan(p->greaterequal) && !isnan(p->lessthan) )
    asprintf(&str, "from (inclusive) %g, upto (exclusive) %g",
             p->greaterequal, p->lessthan);
  else if( !isnan(p->greaterequal) )
    asprintf(&str, "from (inclusive) %g", p->greaterequal);
  else if( !isnan(p->lessthan) )
    asprintf(&str, "upto (exclusive) %g", p->lessthan);
  if(str)
    {
      printf("Range: %s.\n", str);
      free(str);
    }

  /* Units. */
  if(p->input->unit) printf("Unit: %s\n", p->input->unit);

  /* Clean up. */
  if(col) free(col);
  free(name);
  printf("-------\n");
}





/* This function will report the simple immediate statistics of the
   data. For the average and standard deviation, the unsorted data is
   used so we don't suddenly encounter rounding errors. */
void
print_basics(struct statisticsparams *p)
{
  char *str;
  double mean, std;
  int namewidth=40;
  gal_data_t *tmp, *bins, *hist;

  /* Define the input dataset. */
  print_input_info(p);

  /* Print the number: */
  printf("  %-*s %zu\n", namewidth, "Number of elements:", p->input->size);

  /* Minimum: */
  tmp=gal_statistics_minimum(p->input);
  str=gal_data_write_to_string(tmp->array, tmp->type, 0);
  printf("  %-*s %s\n", namewidth, "Minimum:", str);
  gal_data_free(tmp);
  free(str);

  /* Maximum: */
  tmp=gal_statistics_maximum(p->input);
  str=gal_data_write_to_string(tmp->array, tmp->type, 0);
  printf("  %-*s %s\n", namewidth, "Maximum:", str);
  gal_data_free(tmp);
  free(str);

  /* Find the mean and standard deviation, but don't print them, see
     explanations under median. */
  tmp=gal_statistics_mean_std(p->input);
  mean = ((double *)(tmp->array))[0];
  std  = ((double *)(tmp->array))[1];
  gal_data_free(tmp);

  /* Find and print the median: we want the median to be found in place to
     save time/memory. But having a sorted array can decrease the floating
     point accuracy of the standard deviation. So we'll do the median
     calculation in the end.*/
  tmp=gal_statistics_median(p->input, 1);
  str=gal_data_write_to_string(tmp->array, tmp->type, 0);
  printf("  %-*s %s\n", namewidth, "Median:", str);
  gal_data_free(tmp);
  free(str);

  /* Print the mean and standard deviation. */
  printf("  %-*s %g\n", namewidth, "Mean:", mean);
  printf("  %-*s %g\n", namewidth, "Standard deviation:", std);

  /* Ascii histogram. Note that we don't want to force the user to have the
     plotting parameters. */
  printf("-------\nHistogram:\n");
  p->asciiheight = p->asciiheight ? p->asciiheight : 10;
  p->numasciibins = p->numasciibins ? p->numasciibins : 70;
  bins=gal_statistics_regular_bins(p->input, NULL, p->numasciibins, NAN);
  hist=gal_statistics_histogram(p->input, bins, 0, 0);
  print_ascii_plot(p, hist, bins, 1, 0);
  gal_data_free(bins);
  gal_data_free(hist);
}




















/*******************************************************************/
/**************            Sigma clipping            ***************/
/*******************************************************************/
void
print_sigma_clip(struct statisticsparams *p)
{
  float *a;
  char *mode;
  int namewidth=40;
  gal_data_t *sigclip;

  /* Set the mode for printing: */
  if( p->sigclipparam>=1.0f )
    asprintf(&mode, "for %g clips", p->sigclipparam);
  else
    asprintf(&mode, "until relative change in STD is less than %g",
             p->sigclipparam);

  /* Report the status */
  if(!p->cp.quiet)
    {
      print_input_info(p);
      printf("%g-sigma clipping steps %s:\n\n", p->sigclipmultip, mode);
    }

  /* Do the Sigma clipping: */
  sigclip=gal_statistics_sigma_clip(p->sorted, p->sigclipmultip,
                                    p->sigclipparam, p->cp.quiet);
  a=sigclip->array;

  /* Finish the introduction. */
  if(!p->cp.quiet)
    printf("-------\nSummary:\n");
  else
    printf("%g-sigma clipped %s:\n", p->sigclipmultip, mode);

  /* Print the final results: */
  printf("  %-*s %zu\n", namewidth, "Number of input elements:",
         p->input->size);
  if( p->sigclipparam < 1.0f )
    printf("  %-*s %d\n", namewidth, "Number of clips:",     sigclip->status);
  printf("  %-*s %g\n", namewidth, "Final number of elements:", a[0]);
  printf("  %-*s %g\n", namewidth, "Median:",                a[1]);
  printf("  %-*s %g\n", namewidth, "Mean:",                  a[2]);
  printf("  %-*s %g\n", namewidth, "Standard deviation:",    a[3]);

  /* Clean up. */
  free(mode);
}



















/*******************************************************************/
/**************             Main function            ***************/
/*******************************************************************/
void
statistics(struct statisticsparams *p)
{
  int print_basic_info=1;

  /* Print the one-row numbers if the user asked for them. */
  if(p->toprint)
    {
      print_basic_info=0;
      ui_print_one_row(p);
    }

  /* Print the ASCII plots if requested. */
  if(p->asciihist || p->asciicfp)
    {
      ascii_plots(p);
      print_basic_info=0;
    }

  /* Save the histogram and CFP as tables if requested. */
  if(p->histogram || p->cumulative)
    {
      print_basic_info=0;
      save_hist_and_or_cfp(p);
    }

  /* Print the sigma-clipped results. */
  if( !isnan(p->sigclipmultip ) )
    {
      print_basic_info=0;
      print_sigma_clip(p);
    }

  /* If nothing was requested print the simple statistics. */
  if(print_basic_info)
    print_basics(p);

}
