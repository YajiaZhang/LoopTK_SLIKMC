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

#ifndef MATH_AABB_H
#define MATH_AABB_H

#include "vector.h"
#include <errors.h>

/** @file math/AABB.h
 * @ingroup Math
 * \brief Functions defining axis-aligned bounding boxes (AABBs).
 * 
 * An AABB consists of two Vectors (bmin,bmax), such that a Vector x
 * is inside the AABB if \f$ bmin_i \leq x_i \leq bmax_i \f$ for all i.
 */

namespace Math {

/** @addtogroup Math */
/*@{*/

/// Clamps x to be contained in the AABB (bmin,bmax)
inline void AABBClamp(Vector& x,const Vector& bmin,const Vector& bmax,Real d=Zero)
{
  Assert(x.n==bmin.n);
  Assert(x.n==bmax.n);
  for(int i=0;i<x.n;i++)
    x(i) = Clamp(x(i),bmin(i)+d,bmax(i)-d);
}

/// Returns the minimum distance from x to the boundaries of the
/// AABB (bmin,bmax).  Also returns the element index that contains the
/// minimum margin.
inline Real AABBMargin(const Vector& x,const Vector& bmin,const Vector& bmax,int& index)
{
  Assert(x.n==bmin.n);
  Assert(x.n==bmax.n);
  Real margin=Inf;
  index=-1;
  for(int i=0;i<x.n;i++) {
    if(x(i)-bmin(i) < margin) {
      margin = x(i)-bmin(i);
      index=i;
    }
    if(bmax(i)-x(i) < margin) {
      margin = bmax(i)-x(i);
      index=i;
    } 
  }
  return margin;
}

/// Returns the minimum distance from x to the boundaries of the
/// AABB (bmin,bmax)
inline Real AABBMargin(const Vector& x,const Vector& bmin,const Vector& bmax)
{
  Assert(x.n==bmin.n);
  Assert(x.n==bmax.n);
  Real margin=Inf;
  for(int i=0;i<x.n;i++) {
    margin = Min(margin,x(i)-bmin(i));
    margin = Min(margin,bmax(i)-x(i));
  }
  return margin;
}

/// Returns true if the AABB (bmin,bmax) contains x
inline bool AABBContains(const Vector& x,const Vector& bmin,const Vector& bmax,Real d=Zero)
{
  Assert(x.n==bmin.n);
  Assert(x.n==bmax.n);
  for(int i=0;i<x.n;i++)
    if(x(i) < bmin(i)+d || x(i) > bmax(i)+d) return false;
  return true;
}

/// For a point inside the AABB (bmin,bmax), limits the desired step t
/// x = x0+t*dx such that x is inside the AABB
inline void AABBLineSearch(const Vector& x0,const Vector& dx,const Vector& bmin,const Vector& bmax,Real& t)
{
  Assert(x0.n == dx.n);
  Assert(x0.n == bmin.n);
  Assert(x0.n == bmax.n);
  for(int i=0;i<bmax.n;i++) {
    Assert(x0(i) <= bmax(i));
    Assert(x0(i) >= bmin(i));
    if(x0(i) + t*dx(i) > bmax(i)) {
      Assert(dx(i) > Zero);
      t = (bmax(i)-x0(i))/(dx(i)+Epsilon);
    }
    if(x0(i) + t*dx(i) < bmin(i)) {
      Assert(dx(i) < Zero);
      t = (bmin(i)-x0(i))/(dx(i)-Epsilon);
    }
    Assert(x0(i)+t*dx(i) <= bmax(i));
    Assert(x0(i)+t*dx(i) >= bmin(i));
  }
}

/// Clips the line segment x=x0+u*dx, with t in [u0,u1]
/// to the AABB (bmin,bmax).  Returns false if there is no intersection.
inline bool AABBClipLine(const Vector& x0,const Vector& dx,
			 const Vector& bmin,const Vector& bmax,
			 Real& u0,Real& u1);

/*@}*/
} //namespace Math

#endif

