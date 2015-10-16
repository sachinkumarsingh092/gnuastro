/*********************************************************************
Functions for linked lists.
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
#ifndef SLL_H
#define SLL_H



/******************* Two doubles (for coordinates) */
struct tdll
{
    double a;
    double b;
    struct tdll *next;
};

void
add_to_tdll(struct tdll **list, double a, double b);

void
pop_from_tdll(struct tdll **list, double *a, double *b);

size_t
numintdll(struct tdll *list);

void
tdlltoarrayinv(struct tdll *list, double **d, size_t *num);

void
freetdll(struct tdll *list);




















/******************* float: */
struct fll
{
    float v;
    struct fll *next;
};

void
printarrayoffll(struct fll **afll, size_t num);

void
add_to_fll(struct fll **list, float value);

void
pop_from_fll(struct fll **list, float *value);

size_t
numinfll(struct fll *list);

void
flltoarray(struct fll *list, float **f, size_t *num);

void
freefll(struct fll *list);

void
freearrayoffll(struct fll **afll, size_t num);




















/******************* String: */
struct stll
{
    char *v;
    struct stll *next;
};
void
add_to_stll(struct stll **list, char *value);

void
pop_from_stll(struct stll **list, char **value);

void
print_stll(struct stll *list);

size_t
numinstll(struct stll *list);

















/******************* size_t: */
struct sll
{
    size_t v;
    struct sll *next;
};

void
add_to_sll(struct sll **list, size_t value);

void
pop_from_sll(struct sll **list, size_t *value);

size_t
numinsll(struct sll *list);

void
printsll(struct sll *list);

void
slltoarray(struct sll *list, size_t **f, size_t *num);

void
freesll(struct sll *list);




















/******************* Two way size_t: */
struct tsll
{
  size_t v;
  struct tsll *next;
  struct tsll *prev;
};

void
add_to_tsll_end(struct tsll **last, size_t value);

void
pop_from_tsll_start(struct tsll **first,  size_t *value);




















/******************* Ordered size_t: */
struct osll
{
  size_t v;			/* The actual value. */
  float s;			/* The parameter to sort by. */
  struct osll *next;
};

void
add_to_osll(struct osll **list, size_t value, float tosort);

void
pop_from_osll(struct osll **list,  size_t *value, float *sortvalue);

void
osll_into_sll(struct osll *in, struct sll **out);




















/******************* Two way ordered size_t: */
struct tosll
{
  size_t v;
  float s;
  struct tosll *prev;
  struct tosll *next;
};

void
print_tosll(struct tosll *l, struct tosll *s);

void
add_to_tosll_end(struct tosll **largest, struct tosll **smallest,
		  size_t value, float tosort);

void
pop_from_tosll_start(struct tosll **lartest, struct tosll **smallest,
		      size_t *value, float *tosort);

void
smallest_tosll(struct tosll *largest, struct tosll **smallest);

void
tosll_into_sll(struct tosll *in, struct sll **out);

void
tosll_free(struct tosll *largest);

#endif
