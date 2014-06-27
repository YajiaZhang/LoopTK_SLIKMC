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

#ifndef PPCA_H
#define PPCA_H
#include "PLibraries.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_eigen.h"
class PPCA{
    public:
    static Vector3 getMean(Vector3 *points, int number);
    static void normalize( Vector3 mean, Vector3 *data, Vector3 **target, int number);
    static void getCovariance(gsl_matrix *covm, Vector3 *data, int number);
    static void getEigenvector(gsl_matrix *matrix, gsl_vector *eigenvalues, gsl_matrix *eigenvector);
    static void getPrincipalComponent(gsl_vector *pc, Vector3 *points, int number);
    static void getPrincipalComponent(gsl_vector *pc, vector<Vector3*> *points);
};
#endif
