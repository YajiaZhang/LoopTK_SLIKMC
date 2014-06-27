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

#ifndef MATH3D_TRIANGLE2D_H
#define MATH3D_TRIANGLE2D_H

#include "Plane2D.h"

namespace Math3D {

/** @ingroup Math3D
 * @brief A 2D triangle class
 *
 * Represented by its vertices a,b,c.
 *
 * Barycentric coordinates (u,v,w) are such that 0 <= u,v,w <= 1
 * and u+v+w = 1.  They parameterize the triangle as x = u*a+v*b+w*c.
 *
 * "Plane" coordinates (p,q) are such that 0 <= p,q and p+q<= 1.
 * They parameterize the triangle as x = a + p*(b-a) + q*(c-a).
 * Barycentric coordinates (u,v,w) = (1-p-q,p,q).
 */
struct Triangle2D
{
  Triangle2D();
  Triangle2D(const Vector2& a,const Vector2& b,const Vector2& c);
  void set(const Vector2& a,const Vector2& b,const Vector2& c);
  void setTransformed(const Triangle2D& t, const Matrix3& xform);
  
  Real orientation() const;
  Real area() const;
  void getAABB(AABB2D&) const;
  
  Vector3 barycentricCoords(const Point2D& x) const;
  Point2D barycentricCoordsToPoint(const Vector3& bc) const;
  Vector2 planeCoords(const Point2D& x) const;
  Point2D planeCoordsToPoint(const Vector2& pc) const;
  
  Vector2 closestPointCoords(const Point2D& in) const;  ///<returns the plane-coords of the point
  Point2D closestPoint(const Point2D& in) const;
  bool contains(const Point2D& x) const;
  
  bool intersect(const Plane2D&) const;
  bool intersect(const Plane2D&, Segment2D& S) const;
  //bool intersect(const Triangle2D&) const;
  //bool intersect(const Triangle2D&, Segment2D& S) const;
  //edges are (a,b) (b,c) (c,a), u gives the interpolation parameter 
  //of the plane intersection along each (or -1 if no intersection exists)
  //void edgeIntersections(const Plane2D&, Real u[2]) const;
  //void edgeIntersections(const Triangle2D&, Real u[2]) const;

  bool Read(File& f);
  bool Write(File& f) const;
 
  static Real orientation(const Point2D& a, const Point2D& b, const Point2D& c);
  static Real area(const Point2D& a, const Point2D& b, const Point2D& c);
  static Vector3 barycentricCoords(const Vector2& x, const Point2D& a, const Point2D& b, const Point2D& c);
  static Point2D barycentricCoordsToPoint(const Vector3& bc, const Point2D& a, const Point2D& b, const Point2D& c);
  static bool containsBarycentricCoords(const Vector3& bc);
  static Point2D planeCoordsToPoint(const Vector2& pc, const Point2D& a, const Point2D& b, const Point2D& c);
  static bool containsPlaneCoords(const Vector2& pc);
  
  Point2D a,b,c;
};

} //namespace Math3D

#endif
