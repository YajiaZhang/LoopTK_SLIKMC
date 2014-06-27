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

#ifndef MATH3D_PLANE3D_H
#define MATH3D_PLANE3D_H

#include "Line3D.h"
#include "Ray3D.h"
#include "Segment3D.h"

namespace Math3D {

/** @brief A 3D plane class
 * @ingroup Math3D
 *
 * Represents plane with a normal and offset such that x on the plane 
 * satisfy dot(normal,x) = offset.
 */
struct Plane3D
{
  void setPointNormal(const Point3D& a, const Vector3& n);
  void setPointBases(const Point3D& a, const Vector3& b1, const Vector3& b2);
  void setPoints(const Point3D& a, const Point3D& b, const Point3D& c);
  void setTransformed(const Plane3D& pin, const Matrix4& xform);
  
  Real distance(const Point3D& v) const;
  void project(const Point3D& in, Point3D& out) const;		///<projects onto the plane
  void getBasis(Vector3& xb, Vector3& yb) const;		///<returns a basis of the plane
  
  bool intersectsSegment(const Segment3D&, Real* t);
  bool intersectsLine(const Line3D&, Real* t);
  bool intersectsRay(const Ray3D&, Real* t);
  bool intersects(const AABB3D&) const;
  bool intersects(const AABB3D&,Real& dmin,Real& dmax) const;   ///<calculates min/max bbox distances as well
  
  ///returns the dimension of the intersection of the 2 planes. <br>
  ///if 0, they don't intersect. <br>
  ///1, it's a line returned in l <br>
  ///2, the planes are identitcal
  int allIntersections(const Plane3D& p,Line3D& l) const;

  bool Read(File& f);
  bool Write(File& f) const;
  
  Vector3 normal;
  Real offset;
};

} //namespace Math

#endif
