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

namespace Math3D {

template <class VT,int n>
bool LinearlyDependent_Robust_Template(const VT& a, const VT& b, Real& c, bool& cb, Real eps) {
  Real aDotB = a.dot(b);
  Real aNorm2 = a.normSquared();
  if(aNorm2 > Abs(aDotB)) {
    cb = false;
    c = aDotB/aNorm2;
    Real aNorm = Sqrt(aNorm2);
    Real relEps = eps*aNorm;
    for(int i=0;i<n;i++) {
      if(!FuzzyEquals(c*a[i],b[i],relEps)) return false;
    }
    return true;
  }
  else {
    Real bNorm2 = b.normSquared();
    cb = true;
    if(bNorm2 == Zero) {  //both a and b are 0
      c = One;
      return true;
    }
    c = aDotB/bNorm2;
    Real bNorm = Sqrt(bNorm2);
    Real relEps = eps*bNorm;
    for(int i=0;i<n;i++) {
      if(!FuzzyEquals(a[i],c*b[i],relEps)) return false;
    }
    return true;
  }
}

bool LinearlyDependent_Robust(const Vector2& a, const Vector2& b, Real& c, bool& cb, Real eps) {
  return LinearlyDependent_Robust_Template<Vector2,2>(a,b,c,cb,eps);
}

bool LinearlyDependent_Robust(const Vector3& a, const Vector3& b, Real& c, bool& cb, Real eps) {
  return LinearlyDependent_Robust_Template<Vector3,3>(a,b,c,cb,eps);
}

bool LinearlyDependent_Robust(const Vector4& a, const Vector4& b, Real& c, bool& cb, Real eps) {
  return LinearlyDependent_Robust_Template<Vector4,4>(a,b,c,cb,eps);
}

} //namespace Math
