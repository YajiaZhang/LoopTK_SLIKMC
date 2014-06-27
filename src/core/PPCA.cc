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

#include "PPCA.h"

Vector3 PPCA::getMean(Vector3 *points, int number){
      Vector3 mean = Vector3(0,0,0);
      for(int i=0; i<number; i++){
          mean.x += points[i].x/(double)number;
          mean.y += points[i].y/(double)number;
          mean.z += points[i].z/(double)number;
      }
      return mean;    
}

void PPCA::normalize( Vector3 mean, Vector3 *data, Vector3 **target, int number){
    for(int i=0; i<number; i++){
        (*target)[i].x = data[i].x-mean.x;
        (*target)[i].y = data[i].y-mean.y;
        (*target)[i].z = data[i].z-mean.z;
    }
}

void PPCA::getCovariance(gsl_matrix *covm, Vector3 *data, int number){
    

    double xx=0, xy=0, xz=0, yy=0, yz=0, zz=0;
    for(int i=0; i<number; i++){
         xx += data[i].x*data[i].x;
         xy += data[i].x*data[i].y;
         xz += data[i].x*data[i].z;
         yy += data[i].y*data[i].y;
         yz += data[i].y*data[i].z;
         zz += data[i].z*data[i].z;
    }   
    gsl_matrix_set(covm, 0, 0, xx/(number-1));
    gsl_matrix_set(covm, 0, 1, xy/(number-1));
    gsl_matrix_set(covm, 1, 0, xy/(number-1));
    gsl_matrix_set(covm, 0, 2, xz/(number-1));
    gsl_matrix_set(covm, 2, 0, xz/(number-1));
    gsl_matrix_set(covm, 1, 1, yy/(number-1));
    gsl_matrix_set(covm, 1, 2, yz/(number-1));
    gsl_matrix_set(covm, 2, 1, yz/(number-1));
    gsl_matrix_set(covm, 2, 2, zz/(number-1));
}

void PPCA::getEigenvector(gsl_matrix *matrix, gsl_vector *eigenvalues, gsl_matrix *eigenvector){
       gsl_eigen_symmv_workspace * w =
         gsl_eigen_symmv_alloc (3);
       gsl_eigen_symmv (matrix,eigenvalues, eigenvector, w);
       gsl_eigen_symmv_free (w);
       gsl_eigen_symmv_sort (eigenvalues, eigenvector,
                             GSL_EIGEN_SORT_ABS_DESC);
}

void PPCA::getPrincipalComponent(gsl_vector *pc, vector<Vector3*>*points){
    Vector3 *pts = (Vector3*)malloc(sizeof(Vector3)*points->size());
    for(int i=0; i<points->size(); i++) pts[i] = *(*points)[i];
    PPCA::getPrincipalComponent(pc, pts, points->size()); 
    delete pts;
}

void PPCA::getPrincipalComponent(gsl_vector *pc, Vector3 *points, int number){
    
    Vector3 mean = PPCA::getMean(points, number);
    Vector3 *np = (Vector3*)malloc(sizeof(Vector3)*number);
    PPCA::normalize( mean, points, &np, number);
    gsl_matrix *cov = gsl_matrix_alloc(3,3);
    PPCA::getCovariance(cov, np, 10);
    gsl_vector *eigenvalue = gsl_vector_alloc(3);
    gsl_matrix *eigenvector = gsl_matrix_alloc(3, 3);
    PPCA::getEigenvector(cov,eigenvalue, eigenvector);
    gsl_vector_view evec_i = gsl_matrix_column (eigenvector, 0);
    *pc = evec_i.vector;
    delete np;
    delete cov;
    delete eigenvalue;
    delete eigenvector;
}

