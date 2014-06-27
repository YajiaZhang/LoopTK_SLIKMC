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

#include "sample.h"
#include <math/random.h>
#include <iostream>
#include <errors.h>
using namespace std;

namespace Math {

int WeightedSample(const vector<Real>& weights)
{
  Real totWeight=Zero;
  for(size_t i=0;i<weights.size();i++)  totWeight+=weights[i];
  return WeightedSample(weights,totWeight);
}

int WeightedSample(const vector<Real>& weights,Real totWeight)
{
  Real tmp = Rand()*totWeight;
  for(size_t i=0;i<weights.size();i++) {
    tmp -= weights[i];
    if(tmp <= 0) {
      return (int) i;
      break;
    }
  }
  AssertNotReached();
  return -1;
}

Real Sample(const Interval& s)
{
  return Rand(s.a,s.b);
}

Real Sample(const ClosedIntervalSet& s)
{
  Assert(!s.empty());
  Real max=0;
  for(size_t i=0;i<s.size();i++) {
    if(s[i].b < s[i].a) {
      cerr<<"Sample(ClosedInterval): interval "<<i<<" is invalid!"<<endl;
      cerr<<"["<<s[i].a<<","<<s[i].b<<"]"<<endl;
      Abort();
    }
    max += s[i].b-s[i].a;
  }
  if(isinf(max)) {
    cerr<<"Sample(ClosedInterval): interval has infinite size! aborting"<<endl;
    Abort();
  }
  if(max==0) {  //a bunch of single points
    int n=rng.randInt(s.size());
    return s[n].a;
  }
  Real u=Rand()*max;
  for(size_t i=0;i<s.size();i++) {
    u -= (s[i].b-s[i].a);
    if(u < 0) return s[i].b+u;
  }
  cerr<<"Shouldn't get here!"<<endl;
  AssertNotReached();
  return 0;
}


void SampleCircle(Real r, Real& x, Real& y)
{
  Real angle = Rand()*TwoPi;
  x = Cos(angle)*r;
  y = Sin(angle)*r;
}

void SampleDisk(Real r, Real& x, Real& y)
{
  Real angle = Rand()*TwoPi;
  Real dist = Sqrt(Rand());
  x = Cos(angle)*dist*r;
  y = Sin(angle)*dist*r;
}

void SampleTriangle(Real& x, Real& y)
{
  Real u = Rand();
  x = One - Sqrt(One-u);
  y = Rand()*(One-x);
}

void SampleSphere(Real r, Real& x, Real& y, Real& z)
{
  Real theta=Rand()*TwoPi;
  z = Rand()*Two-One;
  Real rad=r*Sqrt(One-Sqr(z));
  x = rad*Cos(theta);
  y = rad*Sin(theta);
  z *= r;
}

void SampleBall(Real r, Real& x, Real& y, Real& z)
{
  SampleSphere(r,x,y,z);
  Real dist = Pow(Rand(),(Real)(1.0/3.0));
  x*=dist;
  y*=dist;
  z*=dist;
}

} //namespace Math
