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

#include "GramSchmidt.h"
#include "complex.h"
#include <iostream>
using namespace std;

namespace Math {

template <class T>
int OrthonormalBasis(const VectorTemplate<T>* x, VectorTemplate<T>* basis, int n)
{
  int k=0;			//counter to number of nonzero output vectors
  VectorTemplate<T> tmp;
  for(int i=0; i<n; i++) {
    if(basis != x)
      basis[k] = x[i];
    for(int j=0; j<k; j++) {
      tmp.mul(basis[j], basis[j].dot(x[i]));
      basis[k] -= tmp;
    }
    T lensquared = basis[k].normSquared();
    if(lensquared != Zero) {
      basis[k].inplaceMul(Inv(Sqrt(lensquared)));
      k++;
    }
    else {
      cout<<"Redundant vector "<<i<<endl;
    }
  }
  return k;
}

template <class T>
int OrthogonalBasis(const VectorTemplate<T>* x, VectorTemplate<T>* basis, int n)
{
  int k=0;			//counter to number of nonzero output vectors
  T* basisSquared = new T[n];
  VectorTemplate<T> tmp;
  for(int i=0; i<n; i++) {
    if(basis != x)
      basis[k] = x[i];
    for(int j=0; j<k; j++) {
      tmp.mul(basis[j], basis[j].dot(x[i])/basisSquared[j]);
      basis[k] -= tmp;
    }
    basisSquared[k] = basis[k].normSquared();
    if(basisSquared[k] != Zero) {
      k++;
    }
    else {
      cout<<"Redundant vector "<<i<<endl;
    }
  }
  delete [] basisSquared;
  return k;
}

template <class T>
void Orthogonalize(VectorTemplate<T>& x,const VectorTemplate<T>* basis, int n)
{
  for(int i=0;i<n;i++)
    x.madd(basis[i],-basis[i].dot(x)/basis[i].normSquared());
}



#define DEFINEGRAMSCHMIDT(T) \
  template int OrthonormalBasis<T>(const VectorTemplate<T>* x, VectorTemplate<T>* basis, int n); \
  template int OrthogonalBasis<T>(const VectorTemplate<T>* x, VectorTemplate<T>* basis, int n); \
  template void Orthogonalize<T>(VectorTemplate<T>& x,const VectorTemplate<T>* basis, int n);
DEFINEGRAMSCHMIDT(float);
DEFINEGRAMSCHMIDT(double);
DEFINEGRAMSCHMIDT(Complex);


} //namespace Math
