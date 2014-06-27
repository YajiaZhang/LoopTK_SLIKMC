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

#ifndef MATH3D_SPHERE3D_H
#define MATH3D_SPHERE3D_H

#include "Line3D.h"

namespace Math3D {

/** @ingroup Math3D
 * @brief A 3D sphere class
 *
 * Represented by a center and a radius.
 * 
 * Most methods consider the sphere as the solid ball.
 * Methods that use the sphere boundary have the prefix "boundary".
 */
struct Sphere3D
{
  Real distance(const Point3D& v) const;
  bool contains(const Point3D& v) const;
  bool contains(const Sphere3D& s) const;
  bool withinDistance(const Point3D& v, Real dist) const;
  bool boundaryWithinDistance(const Point3D& v, Real dist) const;
  bool intersects(const Line3D&, Real* t1=NULL, Real* t2=NULL) const;
  bool intersects(const Ray3D&, Real* t1=NULL, Real* t2=NULL) const;
  bool intersects(const Segment3D&, Real* t1=NULL, Real* t2=NULL) const;
  bool intersects(const Plane3D& p) const;
  bool intersects(const Sphere3D& s) const;
  bool boundaryIntersects(const Sphere3D& s) const;
  bool boundaryIntersectsBoundary(const Sphere3D& s) const;
  bool Read(File& f);
  bool Write(File& f) const;

  void getAABB(AABB3D& bb) const;
  bool intersects(const AABB3D& bb) const;
  bool intersectsApprox(const AABB3D& bb) const;
  
  static bool ballsIntersect(const Point3D& ca,Real ra,const Point3D& cb,Real rb);
  static bool ballSphereIntersect(const Point3D& ca,Real ra,const Point3D& cb,Real rb);
  static bool spheresIntersect(const Point3D& ca,Real ra,const Point3D& cb,Real rb);

  Point3D center;
  Real radius;
};

} //namespace Math3D

#endif
