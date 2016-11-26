/*********************************************************************
Arithmetic - Do arithmetic operations on images.
Arithmetic is part of GNU Astronomy Utilities (Gnuastro) package.

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
#ifndef OPERATORS_H
#define OPERATORS_H


void
sum(struct imgarithparams *p);

void
subtract(struct imgarithparams *p);

void
multiply(struct imgarithparams *p);

void
divide(struct imgarithparams *p);

void
topower(struct imgarithparams *p, char *op);

void
alloppixs(struct imgarithparams *p, char *operator);

void
takesqrt(struct imgarithparams *p);

void
takelog(struct imgarithparams *p);

void
takelog10(struct imgarithparams *p);

void
takeabs(struct imgarithparams *p);

void
findmin(struct imgarithparams *p);

void
findmax(struct imgarithparams *p);

int
lessthan(double left, double right);

int
lessequal(double left, double right);

int
greaterthan(double left, double right);

int
greaterequal(double left, double right);

int
equal(double left, double right);

int
notequal(double left, double right);

void
conditionals(struct imgarithparams *p, char *operator);

void
andor(struct imgarithparams *p, char *operator);

void
notfunc(struct imgarithparams *p);

void
opisblank(struct imgarithparams *p);

void
where(struct imgarithparams *p);

#endif
