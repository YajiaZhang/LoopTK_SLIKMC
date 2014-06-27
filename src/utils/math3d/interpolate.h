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

#ifndef MATH3D_INTERPOLATE_H
#define MATH3D_INTERPOLATE_H

#include "primitives.h"

namespace Math3D {

template <class type>
void interpolate(const type& a, const type& b, Real u, type& x)
{
  x=(One-u)*a + u*b;
}

inline void interpolate(const Vector2& a, const Vector2& b, Real u, Vector2& x)
{
  x.mul(a,One-u);
  x.madd(b,u);
}

inline void interpolate(const Vector3& a, const Vector3& b, Real u, Vector3& x)
{
  x.mul(a,One-u);
  x.madd(b,u);
}

inline void interpolate(const Vector4& a, const Vector4& b, Real u, Vector4& x)
{
  x.mul(a,One-u);
  x.madd(b,u);
}

/*

inline void interpolate(const Matrix2& a, const Matrix2& b, Real u, Matrix2& x)
{
  x.mul(a,One-u);
  x.madd(b,u);
}

inline void interpolate(const Matrix3& a, const Matrix3& b, Real u, Matrix3& x)
{
  x.mul(a,One-u);
  x.madd(b,u);
}

inline void interpolate(const Matrix4& a, const Matrix4& b, Real u, Matrix4& x)
{
  x.mul(a,One-u);
  x.madd(b,u);
}

*/

///if segment is x = a + u*(b-a), returns the value u s.t. x=0
inline Real SegmentZeroCrossing(Real a,Real b)
{
  if(a == b) return Zero;
  return a/(a-b);
}

///if segment is x = a + u*(b-a), returns the value u s.t. x=x0
inline Real SegmentCrossing(Real a,Real b,Real x0)
{
  if(a == b) return Zero;
  return (a-x0)/(a-b);
}

}  //namespace Math3D

#endif
