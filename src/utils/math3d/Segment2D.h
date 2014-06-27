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

#ifndef MATH3D_SEGMENT2D_H
#define MATH3D_SEGMENT2D_H

#include "Line2D.h"

namespace Math3D {

/** @ingroup Math3D
 * @brief A 2D segment class
 *
 * Represented by the endpoints A and B.  Is parameterized by t in [0,1]
 * going from A->B as t approaches 1.
 */
struct Segment2D
{
  void setTransformed(const Segment2D&, const Matrix3& xform);
  void getLine(Line2D& l) const;
  Real closestPointParameter(const Point2D& in) const;
  Real closestPoint(const Point2D& in, Point2D& out) const;  //returns the parameter value of the point
  Real distance(const Point2D& pt) const;
  void eval(Real t, Point2D& out) const;
  Real orientation(const Point2D& p) const;  ///< >0 for left, <0 for right
  bool isLeft(const Point2D& p) const;  ///< p left of A->B
  bool isRight(const Point2D& p) const;  ///< p right of A->B
  bool intersects(const Segment2D& S) const;
  bool intersects(const Vector2& a,const Vector2& b) const;
  bool intersects(const Segment2D& S,Vector2& p) const;
  bool intersects(const Vector2& a,const Vector2& b,Vector2& p) const;
  bool Read(File& f);
  bool Write(File& f) const;

  void getAABB(AABB2D&) const;
  bool intersects(const AABB2D&) const;
  ///given bounds [tmin,tmax] of the segment, returns the clipping min/max
  bool intersects(const AABB2D&, Real& tmin, Real& tmax) const;

  ///orientation of x vs a,b
  static Real orientation(const Point2D& a,const Point2D& b,const Point2D& x);

  Point2D A;
  Point2D B;
};

} //namespace Math3D

#endif
