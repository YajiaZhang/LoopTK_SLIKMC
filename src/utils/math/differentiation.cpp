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

#include "differentiation.h"
#include <vector>
#include <errors.h>
using namespace std;

namespace Math {

Real dfCenteredDifference(RealFunction& f, Real x, Real h)
{
  return (f(x+h)-f(x-h))/(Two*h);
}

void dfCenteredDifference(VectorFunction& f, Real x, Real h, Vector& df)
{
  Vector temp;
  f(x+h,temp);
  df = temp;
  f(x-h,temp);
  df -= temp;
  df /= Two*h;
}

//centered differences for 2nd derivative
//f''(x) ~= 1/h^2*(f(x+h) - 2f(x) + f(x-h)) + O(h^2)
Real ddfCenteredDifference(RealFunction& f, Real x, Real h)
{
  return (f(x+h)-Two*f(x)+f(x-h))/Sqr(h);
}

void ddfCenteredDifference(VectorFunction& f, Real x, Real h, Vector& ddf)
{
  Vector temp;
  f(x+h,temp);
  ddf = temp;
  f(x,temp);
  ddf.madd(temp,-Two);
  f(x-h,temp);
  ddf += temp;
  ddf /= Sqr(h);
}


void GradientForwardDifference(ScalarFieldFunction& f,Vector& x,Real h,Vector& g)
{
  Assert(g.n == x.n);
  Real f0 = f(x),f1;
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    x(i) += h; f1=f(x);
    x(i) = xi;

    g(i) = (f1-f0);
  }
  g /= h;
}

void JacobianForwardDifference(VectorFieldFunction& f,Vector& x,Real h,Matrix& J)
{
  int nd = f.NumDimensions();
  Vector f0(nd),f1(nd);
  Assert(J.m == f.NumDimensions() && J.n == x.n);
  f(x,f0);
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    x(i) += h;  f(x,f1);
    x(i) = xi;

    Vector Ji; J.getColRef(i,Ji);
    Ji.sub(f1,f0);
  }
  J /= h;
}

void HessianForwardDifference(ScalarFieldFunction& f,Vector& x,Real h,Matrix& H)
{
  Assert(H.hasDims(x.n,x.n));
  Real f00,f01,f10,f11;
  f00 = f(x);
  Real scale = Inv(Sqr(h));
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    for(int j=i;j<x.n;j++) {
      Real xj=x(j);
      x(j) += h;   f01 = f(x);
      x(i) += h;   f11 = f(x);
      x(j) = xj;   f10 = f(x);
      x(i) = xi; 

      H(i,j) = H(j,i) = (f11-f10-f01+f00)*scale;
    }
  }
}

void HessianForwardDifference_Grad(ScalarFieldFunction& f,Vector& x,Real h,Matrix& H)
{
  Assert(H.hasDims(x.n,x.n));
  Vector g0(x.n);
  vector<Vector> g(x.n);
  f.PreEval(x);
  f.Gradient(x,g0);
  for(int i=0;i<x.n;i++) {
    g[i].resize(x.n);
    Real xi=x(i);
    x(i)+=h;
    f.PreEval(x);
    f.Gradient(x,g[i]);
    x(i) = xi;
  }
  Real scale = Inv(Two*h);
  for(int i=0;i<x.n;i++) {
    for(int j=i;j<x.n;j++) {
      H(i,j) = H(j,i) = (g[j](i) - g0(i) + g[i](j)-g0(j))*scale;
    }
  }
}


void GradientCenteredDifference(ScalarFieldFunction& f,Vector& x,Real h,Vector& g)
{
  Assert(g.n == x.n);
  Real twoh = h+h;
  Real f0,f1;
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    x(i) += h;   f1=f(x);
    x(i) -= twoh;  f0=f(x);
    x(i) = xi;

    g(i) = (f1-f0);
  }
  g /= twoh;
}

void JacobianCenteredDifference(VectorFieldFunction& f,Vector& x,Real h,Matrix& J)
{
  int nd = f.NumDimensions();
  Vector f0(nd),f1(nd);
  Assert(J.m == f.NumDimensions() && J.n == x.n);
  Real twoh = h+h;
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    x(i) += h;  f(x,f1);
    x(i) -= twoh;  f(x,f0);
    x(i) = xi;

    Vector Ji; J.getColRef(i,Ji);
    Ji.sub(f1,f0);
  }
  J /= twoh;
}

void HessianCenteredDifference(ScalarFieldFunction& f,Vector& x,Real h,Matrix& H)
{
  Assert(H.hasDims(x.n,x.n));
  Real f0,f_1,f_2,f1,f2;
  Real scale=Sqr(Inv(h))*Sqr(Half);
  Real iiscale=Sqr(Inv(h))/12;
  f0 = f(x);
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    //get diagonal entry
    x(i)+=h;   f1 = f(x); 
    x(i)+=h;   f2 = f(x); 
    x(i)=xi-h; f_1 = f(x); 
    x(i)-=h;   f_2 = f(x); 
    x(i) = xi;
    H(i,i) = (Real)(-f2+16.0*f1-30.0*f0+16.0*f_1-f_2)*iiscale;
    /*  This is a test of this 4th order approx vs the std 2nd order
    Real temp = (f1 - Two*f0 + f_1)/Sqr(h);
    printf("f's: %f %f %f %f %f => H: %f... or %f\n",f_2,f_1,f0,f1,f2,H(i,i),temp);
    getchar();
    */

    for(int j=i+1;j<x.n;j++) {
      Real xj = x(j);
      //get dfj,df_j
      x(j) += h;  x(i) += h;   f2 = f(x);
      x(i) -= Two*h;  f1 = f(x);

      x(j) -= Two*h;   f_2 = f(x);
      x(i) += Two*h;   f_1 = f(x);
      x(i) = xi; x(j) = xj;

      H(i,j) = H(j,i) = (f2-f1-f_1+f_2)*scale;
    }
  }
}

void HessianCenteredDifference_Grad(ScalarFieldFunction& f,Vector& x,Real h,Matrix& H)
{
  Assert(H.hasDims(x.n,x.n));
  vector<Vector> g0(x.n);
  vector<Vector> g1(x.n);
  for(int i=0;i<x.n;i++) {
    g0[i].resize(x.n);
    g1[i].resize(x.n);
    Real xi=x(i);
    x(i)+=h;
    f.PreEval(x);
    f.Gradient(x,g1[i]);
    x(i)-=Two*h;
    f.PreEval(x);
    f.Gradient(x,g0[i]);
    x(i) = xi;
  }
  Real scale = Inv(h)*(Real)0.25;
  for(int i=0;i<x.n;i++) {
    for(int j=i;j<x.n;j++) {
      H(i,j) = H(j,i) = (g1[j](i) - g0[j](i) + g1[i](j)-g0[i](j))*scale;
    }
  }
}





void GradientForwardDifference(ScalarFieldFunction& f,Vector& x,const Vector& h,Vector& g)
{
  Assert(g.n == x.n);
  Assert(h.n == x.n);
  Real f0 = f(x),f1;
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    x(i) += h(i); f1=f(x);
    x(i) = xi;

    g(i) = (f1-f0)/h(i);
  }
}

void JacobianForwardDifference(VectorFieldFunction& f,Vector& x,const Vector& h,Matrix& J)
{
  Assert(h.n == x.n);
  int nd = f.NumDimensions();
  Vector f0(nd),f1(nd);
  Assert(J.m == f.NumDimensions() && J.n == x.n);
  f(x,f0);
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    x(i) += h(i);  f(x,f1);
    x(i) = xi;

    Vector Ji; J.getColRef(i,Ji);
    Ji.sub(f1,f0);
    Ji /= h(i);
  }
}

void HessianForwardDifference(ScalarFieldFunction& f,Vector& x,const Vector& h,Matrix& H)
{
  Assert(h.n == x.n);
  Assert(H.hasDims(x.n,x.n));
  Real f00,f01,f10,f11;
  f00 = f(x);
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    for(int j=i;j<x.n;j++) {
      Real xj=x(j);
      x(j) += h(j);   f01 = f(x);
      x(i) += h(i);   f11 = f(x);
      x(j) = xj;      f10 = f(x);
      x(i) = xi;

      H(i,j) = H(j,i) = (f11-f10-f01+f00)/(h(i)*h(j));
    }
  }
}

void HessianForwardDifference_Grad(ScalarFieldFunction& f,Vector& x,const Vector& h,Matrix& H)
{
  Assert(h.n == x.n);
  Assert(H.hasDims(x.n,x.n));
  Vector g0(x.n);
  vector<Vector> g(x.n);
  f.PreEval(x);
  f.Gradient(x,g0);
  for(int i=0;i<x.n;i++) {
    g[i].resize(x.n);
    Real xi=x(i);
    x(i)+=h(i);
    f.PreEval(x);
    f.Gradient(x,g[i]);
    x(i) = xi;
  }
  for(int i=0;i<x.n;i++) {
    for(int j=i;j<x.n;j++) {
      H(i,j) = H(j,i) = Half*((g[j](i) - g0(i))/h(j) + (g[i](j)-g0(j))/h(i));
    }
  }
}


void GradientCenteredDifference(ScalarFieldFunction& f,Vector& x,const Vector& h,Vector& g)
{
  Assert(h.n == x.n);
  Assert(g.n == x.n);
  Real f0,f1;
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    Real hi=h(i);
    Real twoh=hi+hi;
    x(i) += hi;   f1=f(x);
    x(i) -= twoh; f0=f(x);
    x(i) = xi;

    g(i) = (f1-f0)/twoh;
  }
}

void JacobianCenteredDifference(VectorFieldFunction& f,Vector& x,const Vector& h,Matrix& J)
{
  Assert(h.n == x.n);
  int nd = f.NumDimensions();
  Vector f0(nd),f1(nd);
  Assert(J.m == f.NumDimensions() && J.n == x.n);
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    Real hi=h(i);
    Real twoh=hi+hi;
    x(i) += hi;   f(x,f1);
    x(i) -= twoh; f(x,f0);
    x(i) = xi;

    Vector Ji; J.getColRef(i,Ji);
    Ji.sub(f1,f0);
    Ji /= twoh;
  }
}

void HessianCenteredDifference(ScalarFieldFunction& f,Vector& x,const Vector& h,Matrix& H)
{
  Assert(h.n == x.n);
  Assert(H.hasDims(x.n,x.n));
  Real f0,f_1,f_2,f1,f2;
  f0 = f(x); 
  for(int i=0;i<x.n;i++) {
    Real xi=x(i);
    Real hi=h(i);
    //get diagonal entry
    x(i)+=hi;   f1 = f(x); 
    x(i)+=hi;   f2 = f(x); 
    x(i)=xi-hi; f_1 = f(x); 
    x(i)-=hi;   f_2 = f(x); 
    x(i)=xi;
    H(i,i) = (-f2+16*f1-30*f0+16*f_1-f_2)/(12*Sqr(hi));

    for(int j=i+1;j<x.n;j++) {
      Real xj = x(j);
      Real hj = h(j);
      //get dfj,df_j
      x(j) += hj;  x(i) += hi;   f2 = f(x);
      x(i) -= Two*hi;  f1 = f(x);

      x(j) -= Two*hj;   f_2 = f(x);
      x(i) += Two*hi;   f_1 = f(x);
      x(i) = xi; x(j) = xj;

      H(i,j) = H(j,i) = (f2-f1-f_1+f_2)/(hi*hj)*(Real)0.25;
    }
  }
}

void HessianCenteredDifference_Grad(ScalarFieldFunction& f,Vector& x,const Vector& h,Matrix& H)
{
  Assert(h.n == x.n);
  Assert(H.hasDims(x.n,x.n));
  vector<Vector> g0(x.n);
  vector<Vector> g1(x.n);
  for(int i=0;i<x.n;i++) {
    g0[i].resize(x.n);
    g1[i].resize(x.n);
    Real xi=x(i);
    x(i)+=h(i);
    f.PreEval(x);
    f.Gradient(x,g1[i]);
    x(i)-=Two*h(i);
    f.PreEval(x);
    f.Gradient(x,g0[i]);
    x(i)=xi;
  }
  for(int i=0;i<x.n;i++) {
    for(int j=i;j<x.n;j++) {
      H(i,j) = H(j,i) = ((g1[j](i) - g0[j](i))/h(j) + (g1[i](j)-g0[i](j))/h(i))*(Real)0.25;
    }
  }
}

} //namespace math
