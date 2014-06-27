/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef PNUM_ROUTINES_H
#define PNUM_ROUTINES_H

#include "PBasic.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#define NR_END 1
class PNumRoutines{
  public:

   static void nr_svd(double **a, int m, int n, double w[], double **v);
   static void nr_multimin(double p[], int n, double ftol, int *iter, double *fret, FunctFunctor *Funct, DerivFunctor *Deriv);
   static void nr_inverse(double **a, int N, double **y);
   static void nr_inverse( double** a, int N, double ** y, int& status);

  private:
};
#endif
