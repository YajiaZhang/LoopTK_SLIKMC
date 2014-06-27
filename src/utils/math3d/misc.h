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

#ifndef MATH3D_MISC_H
#define MATH3D_MISC_H

#include "primitives.h"
#include <assert.h>

namespace Math3D {

//these utils save a sqrt
template<class V>
inline bool NormLess(const V& x, Real d) { return x.normSquared() < Sqr(d); }
template<class V>
inline bool NormLEQ(const V& x, Real d) { return x.normSquared() <= Sqr(d); }
template<class V>
inline bool NormGreater(const V& x, Real d) { return x.normSquared() > Sqr(d); }
template<class V>
inline bool NormGEQ(const V& x, Real d) { return x.normSquared() <= Sqr(d); }

template<class V>
inline bool DistanceLess(const V& x, const V& y, Real d) { return x.distanceSquared(y) < Sqr(d); }
template<class V>
inline bool DistanceLEQ(const V& x, const V& y, Real d) { return x.distanceSquared(y) <= Sqr(d); }
template<class V>
inline bool DistanceGreater(const V& x, const V& y, Real d) { return x.distanceSquared(y) > Sqr(d); }
template<class V>
inline bool DistanceGEQ(const V& x, const V& y, Real d) { return x.distanceSquared(y) <= Sqr(d); }

/// @brief Orientation test for 2d points
/// @return
/// >0 for p2 left of the line through p0 and p1 <br>
/// =0 for p2 on the line <br>
/// <0 for p2 right of the line
inline Real Orient2D(const Vector2& p0, const Vector2& p1, const Vector2& p2)
{
  return (p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y);
}

inline Real DistanceSquared2D(Real x1,Real y1,Real x2,Real y2)
{
  return Sqr(x1-x2)+Sqr(y1-y2);
}

inline Real Distance2D(Real x1,Real y1,Real x2,Real y2)
{
  return Sqrt(DistanceSquared2D(x1,y1,x2,y2));
}

} //namespace Math3D

#endif
