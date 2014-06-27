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

#include "LinearlyDependent.h"

namespace Math {

template <class T>
bool LinearlyDependent_Robust(const VectorTemplate<T>& a, const VectorTemplate<T>& b, T& c, bool& cb, T eps) 
{
  assert(a.n == b.n);
  T aDotB = a.dot(b);
  T aNorm2 = a.normSquared();
  if(aNorm2 > Abs(aDotB)) {
    cb = false;
    c = aDotB/aNorm2;
    T aNorm = Sqrt(aNorm2);
    T relEps = eps*aNorm;
    for(int i=0;i<a.n;i++) {
      if(!FuzzyEquals(c*a[i],b[i],relEps)) return false;
    }
    return true;
  }
  else {
    T bNorm2 = b.normSquared();
    cb = true;
    if(bNorm2 == Zero) {  //both a and b are 0
      c = One;
      return true;
    }
    c = aDotB/bNorm2;
    T bNorm = Sqrt(bNorm2);
    T relEps = eps*bNorm;
    for(int i=0;i<a.n;i++) {
      if(!FuzzyEquals(a[i],c*b[i],relEps)) return false;
    }
    return true;
  }
}


} //namespace Math
