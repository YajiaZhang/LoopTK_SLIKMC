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

#ifndef MATH3D_LINE3D_H
#define MATH3D_LINE3D_H

#include "AABB3D.h"

namespace Math3D {

using namespace Math; 

struct Segment3D;

/** @ingroup Math3D
 * @brief A 3D line class
 *
 * A redundant representation, using a point s on the line and a direction d.
 * Is parameterized by x = s+t*d for all real t.
 */
struct Line3D
{
  void setPoints(const Point3D& a, const Point3D& b);
  void setSegment(const Segment3D& s);
  void setTransformed(const Line3D&, const Matrix4& xform);
  Real closestPointParameter(const Point3D& in) const;
  Real closestPoint(const Point3D& in, Point3D& out) const;  ///<returns the parameter value of the point
  Real closestPoint(const Point3D& in, Point3D& out, Real tmin, Real tmax) const;  //tmin,tmax limit the range of the parameter t
  Real distance(const Point3D& pt) const;
  void eval(Real t, Point3D& out) const;
  bool intersects(const Line3D&, Real* t=NULL, Real* u=NULL, Real epsilon=0) const;	///<t is the parameter of this line, u is the other line 
  void closestPoint(const Line3D&,Real& t,Real& u) const;  ///<same comment as above
  void getAABB(AABB3D&, Real tmin=-Inf, Real tmax=Inf) const;
  bool lineIntersects(const AABB3D&) const;
  bool rayIntersects(const AABB3D&) const;
  ///given bounds [tmin,tmax] of the line, returns the clipping min/max
  bool intersects(const AABB3D&, Real& tmin, Real& tmax) const;
  bool Read(File& f);
  bool Write(File& f) const;

  Point3D source;
  Vector3 direction;
};

} //namespace Math3D

#endif
