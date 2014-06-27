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

#include "metric.h"
#include "complex.h"
#include <errors.h>

namespace Math {

template <class T>
inline T DotSelf(const T& a) { return dot(a,a); }

template <class T>
T Norm_L1(const VectorTemplate<T>& x)
{
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum += Abs(x[i]);
  return sum;
}

template <class T>
T Norm_L2(const VectorTemplate<T>& x)
{
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum += DotSelf(x[i]);
  return Sqrt(sum);
}

template <class T>
T Norm_L2_Safe(const VectorTemplate<T>& x)
{
  //factor out maximum element to avoid overflow
  T xmax = x.maxAbsElement();
  if(xmax == 0) return 0;
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum += DotSelf(x[i]/xmax);
  return xmax*Sqrt(sum);
}

template <class T>
T Norm_LInf(const VectorTemplate<T>& x)
{
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum = Max(Abs(x[i]),sum);
  return sum;
}

template <class T>
T Norm_Mahalanobis(const VectorTemplate<T>& x,T k)
{
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum += Pow(x[i],k);
  return Pow(sum,Inv(k));
}


template<class T>
T Distance_L1(const VectorTemplate<T>& x,const VectorTemplate<T>& y)
{
  Assert(x.n == y.n);
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum += Abs(x[i]-y[i]);
  return sum;
}

template<class T>
T Distance_L2(const VectorTemplate<T>& x,const VectorTemplate<T>& y)
{
  Assert(x.n == y.n);
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum += DotSelf(x[i]-y[i]);
  return Sqrt(sum);
}

template<class T>
T Distance_L2_Safe(const VectorTemplate<T>& x,const VectorTemplate<T>& y)
{
  Assert(x.n == y.n);
  T dmax = 0;
  for(int i=0;i<x.n;i++) dmax = Max(dmax,Abs(x[i]-y[i]));
  if(dmax == 0) return 0;
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum += DotSelf((x[i]-y[i])/dmax);
  return dmax*Sqrt(sum);
}

template<class T>
T Distance_LInf(const VectorTemplate<T>& x,const VectorTemplate<T>& y)
{
  Assert(x.n == y.n);
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum = Max(Abs(x[i]-y[i]),sum);
  return sum;
}

template<class T>
T Distance_Mahalanobis(const VectorTemplate<T>& x,const VectorTemplate<T>& y,T k)
{
  Assert(x.n == y.n);
  T sum=0;
  for(int i=0;i<x.n;i++)
    sum += Pow(x[i]-y[i],k);
  return Pow(sum,Inv(k));
}



template <class T>
T Norm_L1(const MatrixTemplate<T>& A)
{
  T v=0;
  for(int j=0;j<A.n;j++) {
    T sum=0;
    for(int i=0;i<A.m;i++) sum += Abs(A(i,j));
    v = Max(v,sum);
  }
  return v;
}

template <class T>
T Norm_LInf(const MatrixTemplate<T>& A)
{
  T v=0;
  for(int i=0;i<A.m;i++) {
    T sum=0;
    for(int j=0;j<A.m;j++) sum += Abs(A(i,j));
    v = Max(v,sum);
  }
  return v;
}

template <class T>
T Norm_Frobenius(const MatrixTemplate<T>& A)
{
  T sum=0;
  MatrixIterator<T> a = A.begin();
  for(int i=0;i<A.m;i++,a.nextRow())
    for(int j=0;j<A.n;j++,a.nextCol())
      sum += DotSelf(*a);
  return Sqrt(sum);
}

template <class T>
T Norm_Frobenius_Safe(const MatrixTemplate<T>& A)
{
  MatrixIterator<T> a = A.begin();
  T vmax=A.maxAbsElement();
  if(vmax == 0) return 0;
  T sum=0;
  a=A.begin();
  for(int i=0;i<A.m;i++,a.nextRow())
    for(int j=0;j<A.n;j++,a.nextCol())
      sum += DotSelf(*a/vmax);
  return vmax*Sqrt(sum);
}


template<class T>
T Distance_L1(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B)
{
  Assert(A.hasDims(B.m,B.n));
  T v=0;
  for(int j=0;j<A.n;j++) {
    T sum=0;
    for(int i=0;i<A.m;i++) sum += Abs(A(i,j)-B(i,j));
    v = Max(v,sum);
  }
  return v;
}


template<class T>
T Distance_LInf(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B)
{
  T v=0;
  for(int i=0;i<A.m;i++) {
    T sum=0;
    for(int j=0;j<A.m;j++) sum += Abs(A(i,j)-B(i,j));
    v = Max(v,sum);
  }
  return v;
}

template<class T>
T Distance_Frobenius(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B)
{
  Assert(A.hasDims(B.m,B.n));
  MatrixIterator<T> a = A.begin();
  MatrixIterator<T> b = B.begin();
  T sum=0;
  for(int i=0;i<A.m;i++,a.nextRow(),b.nextRow())
    for(int j=0;j<A.n;j++,a.nextCol(),b.nextCol())
      sum += DotSelf(*a-*b);
  return Sqrt(sum);
}

template<class T>
T Distance_Frobenius_Safe(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B)
{
  MatrixIterator<T> a = A.begin();
  MatrixIterator<T> b = B.begin();
  T vmax=0;
  for(int i=0;i<A.m;i++,a.nextRow(),b.nextRow())
    for(int j=0;j<A.n;j++,a.nextCol(),b.nextCol())
      vmax = Max(vmax,Abs(*a-*b));
  if(vmax == 0) return 0;
  T sum=0;
  for(int i=0;i<A.m;i++,a.nextRow(),b.nextRow())
    for(int j=0;j<A.n;j++,a.nextCol(),b.nextCol())
      sum += Sqr((*a-*b)/vmax);
  return vmax*Sqrt(sum);
}



//instantiations for complex
template <> Complex Norm_LInf(const VectorTemplate<Complex>& x)
{
  Real sum=Zero;
  for(int i=0;i<x.n;i++)
    sum = Max(Abs(x[i]),sum);
  return sum;
}

template <> Complex Norm_Mahalanobis(const VectorTemplate<Complex>& x,Complex k)
{
  Assert(k.y == Zero);
  Complex sum=Zero;
  for(int i=0;i<x.n;i++)
    sum += Pow(x[i],k.x);
  return Pow(sum,Inv(k.x));
}

template <> Complex Distance_LInf(const VectorTemplate<Complex>& x,const VectorTemplate<Complex>& y)
{
  Assert(x.n == y.n);
  Real sum=Zero;
  for(int i=0;i<x.n;i++)
    sum = Max(Abs(x[i]-y[i]),sum);
  return sum;
}

template <> Complex Distance_Mahalanobis(const VectorTemplate<Complex>& x,const VectorTemplate<Complex>& y,Complex k)
{
  Assert(x.n == y.n);
  Assert(k.y == Zero);
  Complex sum=Zero;
  for(int i=0;i<x.n;i++)
    sum += Pow(x[i]-y[i],k.x);
  return Pow(sum,Inv(k.x));
}

template <> Complex Distance_L2_Safe(const VectorTemplate<Complex>& x,const VectorTemplate<Complex>& y)
{
  Assert(x.n == y.n);
  Real dmax = Zero;
  for(int i=0;i<x.n;i++) dmax = Max(dmax,Abs(x[i]-y[i]));
  if(dmax == Zero) return Zero;
  Complex sum=Zero;
  for(int i=0;i<x.n;i++)
    sum += DotSelf((x[i]-y[i])/dmax);
  return dmax*Sqrt(sum);
}

template <> Complex Norm_L1(const MatrixTemplate<Complex>& A)
{
  Real v=Zero;
  for(int j=0;j<A.n;j++) {
    Real sum=Zero;
    for(int i=0;i<A.m;i++) sum += Abs(A(i,j));
    v = Max(v,sum);
  }
  return v;
}

template <> Complex Norm_LInf(const MatrixTemplate<Complex>& A)
{
  Real v=Zero;
  for(int i=0;i<A.m;i++) {
    Real sum=Zero;
    for(int j=0;j<A.m;j++) sum += Abs(A(i,j));
    v = Max(v,sum);
  }
  return v;
}

template <> Complex Distance_L1(const MatrixTemplate<Complex>& A,const MatrixTemplate<Complex>& B)
{
  Assert(A.hasDims(B.m,B.n));
  Real v=Zero;
  for(int j=0;j<A.n;j++) {
    Real sum=Zero;
    for(int i=0;i<A.m;i++) sum += Abs(A(i,j)-B(i,j));
    v = Max(v,sum);
  }
  return v;
}


template <> Complex Distance_LInf(const MatrixTemplate<Complex>& A,const MatrixTemplate<Complex>& B)
{
  Real v=Zero;
  for(int i=0;i<A.m;i++) {
    Real sum=Zero;
    for(int j=0;j<A.m;j++) sum += Abs(A(i,j)-B(i,j));
    v = Max(v,sum);
  }
  return v;
}

template<> Complex Distance_Frobenius_Safe(const MatrixTemplate<Complex>& A,const MatrixTemplate<Complex>& B)
{
  MatrixIterator<Complex> a = A.begin();
  MatrixIterator<Complex> b = B.begin();
  Real vmax=0;
  for(int i=0;i<A.m;i++,a.nextRow(),b.nextRow())
    for(int j=0;j<A.n;j++,a.nextCol(),b.nextCol())
      vmax = Max(vmax,Abs(*a-*b));
  if(vmax == Zero) return Zero;
  Complex sum=Zero;
  for(int i=0;i<A.m;i++,a.nextRow(),b.nextRow())
    for(int j=0;j<A.n;j++,a.nextCol(),b.nextCol())
      sum += Sqr((*a-*b)/vmax);
  return vmax*Sqrt(sum);
}


#define DEFINEMETRIC(T) \
  template T Norm_L1(const VectorTemplate<T>& x); \
  template T Norm_L2(const VectorTemplate<T>& x); \
  template T Norm_L2_Safe(const VectorTemplate<T>& x); \
  template T Norm_LInf(const VectorTemplate<T>& x); \
  template T Norm_Mahalanobis(const VectorTemplate<T>& x,T k); \
  template T Distance_L1(const VectorTemplate<T>& x,const VectorTemplate<T>& y); \
  template T Distance_L2(const VectorTemplate<T>& x,const VectorTemplate<T>& y); \
  template T Distance_L2_Safe(const VectorTemplate<T>& x,const VectorTemplate<T>& y); \
  template T Distance_LInf(const VectorTemplate<T>& x,const VectorTemplate<T>& y); \
  template T Distance_Mahalanobis(const VectorTemplate<T>& x,const VectorTemplate<T>& y,T k); \
  template T Norm_L1(const MatrixTemplate<T>& A); \
  template T Norm_LInf(const MatrixTemplate<T>& A); \
  template T Norm_Frobenius(const MatrixTemplate<T>& A); \
  template T Norm_Frobenius_Safe(const MatrixTemplate<T>& A); \
  template T Distance_L1(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B); \
  template T Distance_LInf(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B); \
  template T Distance_Frobenius(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B); \
  template T Distance_Frobenius_Safe(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B);

DEFINEMETRIC(float);
DEFINEMETRIC(double);
DEFINEMETRIC(Complex);


} //namespace Math
