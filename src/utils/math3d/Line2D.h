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

#ifndef MATH3D_LINE2D_H
#define MATH3D_LINE2D_H

#include "primitives.h"
#include "AABB2D.h"

namespace Math3D {

using namespace Math; 

struct Segment2D;

/** @ingroup Math3D
 * @brief A 2D line class
 *
 * A redundant representation, using a point s on the line and a direction d.
 * Is parameterized by x = s+t*d for all real t.
 */
struct Line2D
{
  void setPoints(const Point2D& a, const Point2D& b);
  void setSegment(const Segment2D& s);
  void setTransformed(const Line2D&, const Matrix3& xform);
  Real closestPointParameter(const Point2D& in) const;
  Real closestPoint(const Point2D& in, Point2D& out) const;  ///<returns the parameter value of the point
  Real closestPoint(const Point2D& in, Point2D& out, Real tmin, Real tmax) const;  ///<tmin,tmax limit the range of the parameter t
  Real distance(const Point2D& pt) const;
  void eval(Real t, Point2D& out) const;
  Real orientation(const Point2D& p) const;
  bool isLeft(const Point2D& p) const { return orientation(p) > Zero; }
  bool isRight(const Point2D& p) const { return orientation(p) < Zero; }
  bool Read(File& f);
  bool Write(File& f) const;

  void getAABB(AABB2D&, Real tmin=-Inf, Real tmax=Inf) const;
  bool lineIntersects(const AABB2D&) const;
  bool rayIntersects(const AABB2D&) const;
  ///given bounds [tmin,tmax] of the line, returns the clipping min/max
  bool intersects(const AABB2D&, Real& tmin, Real& tmax) const;

  Point2D source;
  Vector2 direction;
};

} //namespace Math3D

#endif
