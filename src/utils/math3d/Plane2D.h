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

#ifndef MATH3D_PLANE2D_H
#define MATH3D_PLANE2D_H

#include "Segment2D.h"
#include "Ray2D.h"

namespace Math3D {

/** @brief A 2D plane class
 * @ingroup Math3D
 *
 * Represents plane with a normal and offset such that x on the plane 
 * satisfy dot(normal,x) = offset.
 */
struct Plane2D
{
  void setPointNormal(const Point2D& a, const Vector2& n);
  void setLine(const Line2D& l);
  void setPoints(const Point2D& a, const Point2D& b);
  void setTransformed(const Plane2D& pin, const Matrix3& xform);
  
  Real distance(const Point2D& v) const;
  void project(const Point2D& in, Point2D& out) const;		///<projects onto the plane
  void getBasis(Vector2& xb) const;            ///<returns a plane basis vector
  
  bool intersectsSegment(const Segment2D&, Real* t);
  bool intersectsLine(const Line2D&, Real* t);
  bool intersectsRay(const Ray2D&, Real* t);
  bool intersects(const AABB2D&) const;
  
  ///returns the dimension of the intersection of the 2 planes. <br>
  ///if 0, they don't intersect, pt is intersection "point at infinty". <br>
  ///1, it's a point returned in pt <br>
  ///2, the planes are identical
  int allIntersections(const Plane2D& p,Vector2& pt) const;

  bool Read(File& f);
  bool Write(File& f) const;

  Vector2 normal;
  Real offset;
};

} //namespace Math3D

#endif
