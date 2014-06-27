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

#include "backsubstitute.h"
#include "complex.h"
//#include "MatrixPrinter.h"
#include <errors.h>
using namespace std;

namespace Math {

const static Real kBackSubZeroTolerance = (Real)1e-4;

template <class T>
bool UBackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x)
{
  Assert(a.isSquare());
  Assert(a.n == b.n);
  Assert(a.n == x.n);
  int n=a.n;
	T aii,sum;
	for(int i=n-1; i>=0; i--) {
		aii=a(i,i);
		sum=b[i];
		for(int j=i+1; j<n; j++)
			sum-=a(i,j)*x[j];
		if(aii == 0) {
		  if(!FuzzyZero(sum,(T)kBackSubZeroTolerance)) {
		    cerr<<"UBackSubstitute: dividing by zero: "<<sum<<"/"<<aii<<endl;
		    return false;
		  }
		  x[i]=0;
		}
		else
		  x[i]=sum/aii;
	}
	return true;
}

// If A is lower triangular nxn, solves Ax=b
template <class T>
bool LBackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x)
{
  Assert(a.isSquare());
  Assert(a.n == b.n);
  Assert(a.n == x.n);
  int n=a.n;
	T aii,sum;
	for(int i=0; i<n; i++) {
		aii=a(i,i);
		sum=b[i];
		for(int j=0; j<i; j++)
			sum-=a(i,j)*x[j];
		if(aii == 0) {
		  if(!FuzzyZero(sum,(T)kBackSubZeroTolerance)) {
		    cerr<<"LBackSubstitute: dividing by zero: "<<sum<<"/"<<aii<<endl;
//		    cerr<<MatrixPrinter(a)<<endl;
		    return false;
		  }
		  x[i]=0;
		}
		else
		  x[i]=sum/aii;
	}
	return true;
}

// If A is lower triangular nxn, solves A^t*x=b
template <class T>
bool LtBackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x)
{ 
  Assert(a.isSquare());
  Assert(a.n == b.n);
  Assert(a.n == x.n);
  int n=a.n;
 	T aii,sum;
	for(int i=n-1; i>=0; i--) {
		aii=a(i,i);
		sum=b[i];
		for(int j=i+1; j<n; j++)
			sum-=a(j,i)*x[j];
		if(aii == 0) {
		  if(!FuzzyZero(sum,(T)kBackSubZeroTolerance)) {
		    cerr<<"LtBackSubstitute: dividing by zero: "<<sum<<"/"<<aii<<endl;
		    return false;
		  }
		  x[i]=0;
		}
		else
		  x[i]=sum/aii;
	}
	return true;
}

template <class T>
void U1BackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x)
{ 
  Assert(a.isSquare());
  Assert(a.n == b.n);
  Assert(a.n == x.n);
  int n=a.n;
  T sum;
	for(int i=n-1; i>=0; i--) {
		sum=b(i);
		for(int j=i+1; j<n; j++)
			sum-=a(i,j)*x[j];
    x[i]=sum;
	}
}

template <class T>
void L1BackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x)
{
  Assert(a.isSquare());
  Assert(a.n == b.n);
  Assert(a.n == x.n);
  int n=a.n;
	T sum;
	for(int i=0; i<n; i++) {
		sum=b[i];
		for(int j=0; j<i; j++)
			sum-=a(i,j)*x[j];
		x[i]=sum;
	}
}

template <class T>
void Lt1BackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x)
{
  Assert(a.isSquare());
  Assert(a.n == b.n);
  Assert(a.n == x.n);
  int n=a.n;
  T sum;
	for(int i=n-1; i>=0; i--) {
		sum=b[i];
		for(int j=i+1; j<n; j++)
			sum-=a(j,i)*x[j];
		x[i]=sum;
	}
}

#define DEFINEBACKSUBSTITUTE(T) \
template bool UBackSubstitute<T>(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x); \
template bool LBackSubstitute<T>(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x); \
template bool LtBackSubstitute<T>(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x); \
template void U1BackSubstitute<T>(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x); \
template void L1BackSubstitute<T>(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x); \
template void Lt1BackSubstitute<T>(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x);

DEFINEBACKSUBSTITUTE(float);
DEFINEBACKSUBSTITUTE(double);
DEFINEBACKSUBSTITUTE(Complex);

} //namespace Math
