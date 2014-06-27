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

#include "LDL.h"
#include "backsubstitute.h"
#include "DiagonalMatrix.h"
#include <errors.h>
using namespace std;

namespace Math {

template<class T>
LDLDecomposition<T>::LDLDecomposition()
  :zeroTolerance((T)1e-8)
{}

template<class T>
LDLDecomposition<T>::LDLDecomposition(const MatrixT& A)
  :zeroTolerance((T)1e-8)
{
  set(A);
}

template<class T>
void LDLDecomposition<T>::set(const MatrixT& A)
{
  Assert(A.m == A.n);
  LDL.resize(A.n,A.n);
  int i,j,k;
  T sum;
  for(i=0;i<A.n;i++) {
    sum = A(i,i);
    for(k=0;k<i;k++) sum -= LDL(k,k)*Sqr(LDL(i,k));
    LDL(i,i) = sum;
    if(FuzzyZero(sum,zeroTolerance)) {
      cerr<<"Warning: LDLt decomposition has a zero element on diagonal "<<i<<endl;
    }

    for(j=i+1;j<A.n;j++) {
      sum = A(i,j);
      for(int k=0;k<i;k++) sum -= LDL(k,k)*LDL(i,k)*LDL(j,k);
      if(LDL(i,i) == 0) {
	if(sum != 0) {
	  cerr<<"Zero diagonal, what to do with sum "<<sum<<endl;
	  sum = 0;
	}
      }
      else 
	sum /= LDL(i,i);
      LDL(j,i) = LDL(i,j) = sum;
    }
  }


  /*
  MatrixT L,LD,LDLt;
  VectorT D;
  getL(L);
  getD(D);
  //cout<<"A: "; A.print();
  //cout<<"L: "; L.print();
  //cout<<"D: "; D.print();
  LD = L;
  for(int i=0;i<A.n;i++)
    LD.scaleCol(i,LDL(i,i));
  LDLt.mulTransposeB(LD,L);
  //cout<<"LDLt: "; LDLt.print();
  LDLt -= A;
  cout<<"Max error in LDL "<<LDLt.maxAbsElement()<<endl;
  */
}

template<class T>
void LDLDecomposition<T>::backSub(const VectorT& b, VectorT& x) const
{
  //LDLt*x=b
  //DLt*x=L^-1*b=y
  //Lt*x=D^-1*y=y'
  //x=(Lt^-1)y
  VectorT y;
  LBackSub(b,y);
  DBackSub(y,y);
  LTBackSub(y,x);
}

template<class T>
void LDLDecomposition<T>::LBackSub(const VectorT& b, VectorT& x) const
{
  Assert(b.n == LDL.n);
  x.resize(LDL.n);
  L1BackSubstitute(LDL,b,x);
}

template<class T>
void LDLDecomposition<T>::DBackSub(const VectorT& b, VectorT& x) const
{
  x.resize(b.n);
  Assert(b.n==x.n);
  for(int i=0;i<x.n;i++) {
    if(!FuzzyZero(LDL(i,i),zeroTolerance))
      x(i) = b(i)/LDL(i,i);
    else {
      if(!FuzzyZero(b(i),zeroTolerance))
	cerr<<"LDLDecomposition::DBackSub(): Warning, zero on the diagonal"<<endl;
      x(i) = 0;
    }
  }
}

template<class T>
void LDLDecomposition<T>::LTBackSub(const VectorT& b, VectorT& x) const
{
  Assert(b.n == LDL.n);
  x.resize(LDL.n);
  Lt1BackSubstitute(LDL,b,x);
}

template<class T>
void LDLDecomposition<T>::getInverse(MatrixT& Ainv) const
{
  Ainv.resize(LDL.n,LDL.n);
  VectorT temp(LDL.n,Zero),y,x;
  for(int i=0;i<LDL.n;i++) {
    temp(i)=One;
    LBackSub(temp,y);
    DBackSub(y,y);
    LTBackSub(y,x);
    //fill in a column
    for(int j=0;j<LDL.n;j++)
      Ainv(j,i)=x(j);
    temp(i)=Zero;
  }
}

template<class T>
void LDLDecomposition<T>::getL(MatrixT& L) const
{
  Assert(LDL.m == LDL.n);
  L.resize(LDL.m,LDL.n);
  for(int i=0;i<LDL.n;i++) {
    L(i,i) = One;
    for(int j=0;j<i;j++)
      L(i,j) = LDL(i,j);
    for(int j=i+1;j<LDL.n;j++)
      L(i,j) = Zero;
  }
}

template<class T>
void LDLDecomposition<T>::getD(VectorT& d) const
{
  Assert(LDL.m == LDL.n);
  d.resize(LDL.n);
  LDL.getDiagCopy(0,d);
}

template <class T>
void LDLDecomposition<T>::mulL(const Vector& x,Vector& y) const
{
  int n=LDL.n;
  Assert(x.n == n);
  y.resize(n);
  for(int i=0;i<n;i++) {
    Real sum = x(i);  //Lii = 1
    for(int j=0;j<i;j++)
      sum += LDL(i,j)*x(j);
    y(i) = sum;
  }
}

template <class T>
void LDLDecomposition<T>::mulLT(const Vector& x,Vector& y) const
{
  int n=LDL.n;
  Assert(x.n == n);
  y.resize(n);
  for(int i=0;i<n;i++) {
    Real sum = x(i);  //Lii = 1
    for(int j=i+1;j<n;j++)
      sum += LDL(j,i)*x(j);
    y(i) = sum;
  }
}

template <class T>
void LDLDecomposition<T>::mulD(const Vector& x,Vector& y) const
{
  int n=LDL.n;
  Assert(x.n == n);
  y.resize(n);
  for(int i=0;i<n;i++) y(i) = x(i)*LDL(i,i);
}

template <class T>
void LDLDecomposition<T>::getA(MatrixT& A) const
{
  MatrixT L,temp;
  DiagonalMatrixTemplate<T> D;
  getL(L);
  getD(D);
  D.postMultiply(L,temp);
  A.mulTransposeB(temp,L);
}

template<class T>
void LDLDecomposition<T>::update(const VectorT& _x)
{
  VectorT x = _x;  //make a copy, we'll change it
  int n=LDL.n;
  Assert(x.n == n);

  T alpha=1;
  for(int i=0;i<n;i++) {
    T deltai = LDL(i,i);
    T temp = alpha + Sqr(x(i))/deltai;
    deltai = deltai*temp;
    T gamma = x(i)/deltai;
    deltai = deltai / alpha;
    alpha = temp;
    LDL(i,i) = deltai;
    for(int k=i+1;k<n;k++) {
      x(k) -= x(i)*LDL(k,i);
      LDL(k,i) += gamma*x(k);
    }
  }
}

template <class T>
bool LDLDecomposition<T>::downdate(const VectorT& _x)
{
  VectorT x = _x;  //make a copy, we'll change it
  int n=LDL.n;
  Assert(x.n == n);

  T alpha=1;
  for(int i=0;i<n;i++) {
    T deltai = LDL(i,i);
    T temp = alpha - Sqr(x(i))/deltai;
    deltai = deltai*temp;
    if(deltai == 0) return false;
    T gamma = x(i)/deltai;
    deltai = deltai / alpha;
    alpha = temp;
    LDL(i,i) = deltai;
    for(int k=i+1;k<n;k++) {
      x(k) -= x(i)*LDL(k,i);
      LDL(k,i) -= gamma*x(k);
    }
  }
  return true;
}

template class LDLDecomposition<float>;
template class LDLDecomposition<double>;
//template class LDLDecomposition<Complex>;

} //namespace Math
