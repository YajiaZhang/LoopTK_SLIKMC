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

#include "function.h"
#include "vectorfunction.h"
#include <errors.h>
using namespace Math;
using namespace std;


Real ScalarFieldFunction::DirectionalDeriv(const Vector& x,const Vector& h)
{
  Vector grad;
  Gradient(x,grad);
  return grad.dot(h);
}

Real ScalarFieldFunction::DirectionalDeriv2(const Vector& x,const Vector& h)
{
  cerr<<"ScalarFieldFunction::DirectionalDeriv2: Warning, possibly inefficient evaluation"<<endl;
  Matrix H(x.n,x.n);
  Hessian(x,H);
  //calc h^t H h
  Real d=Zero;
  for(int i=0;i<x.n;i++) {
    d += h(i)*H.dotRow(i,h);
  }
  return d;
}





std::string VectorFieldFunction::Label() const
{
  return "<unknown Rm->Rn>";
}

std::string VectorFieldFunction::Label(int i) const
{
  std::string str = Label();
  char buf[32];
  sprintf(buf,"[%d]",i);
  str += buf;
  return str;
}

Real VectorFieldFunction::Eval_i(const Vector& x,int i)
{
  cout<<"Warning: really inefficient Eval_i"<<endl;
  Vector v(NumDimensions());
  Eval(x,v);
  return v(i);
}

void VectorFieldFunction::Jacobian_i(const Vector& x,int i,Vector& Ji)
{
  Ji.resize(x.n);
  for(int j=0;j<Ji.n;j++)
    Ji(j) = Jacobian_ij(x,i,j);
}

void VectorFieldFunction::Jacobian_j(const Vector& x,int j,Vector& Jj)
{
  Jj.resize(NumDimensions());
  for(int i=0;i<Jj.n;i++)
    Jj(i) = Jacobian_ij(x,i,j);
}

void VectorFieldFunction::Jacobian(const Vector& x,Matrix& J)
{
  J.resize(NumDimensions(),x.n);
  for(int i=0;i<J.m;i++)
    for(int j=0;j<J.n;j++)
      J(i,j) = Jacobian_ij(x,i,j);
}

void VectorFieldFunction::DirectionalDeriv(const Vector& x,const Vector& h,Vector& v)
{
  Matrix J;
  Jacobian(x,J);
  J.mul(h,v);
}

Real VectorFieldFunction::Divergence(const Vector& x)
{
  Assert(x.n==NumDimensions());
  Real sum=Zero;
  for(int i=0;i<x.n;i++) sum+=Jacobian_ij(x,i,i);
  return sum;
}





