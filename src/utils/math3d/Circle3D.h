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

#ifndef MATH3D_CIRCLE3D_H
#define MATH3D_CIRCLE3D_H

#include "Line3D.h"
#include "Plane3D.h"
#include "Sphere3D.h"

namespace Math3D {

/** @ingroup Math3D
 * @brief A 2D circle in 3D space class
 *
 * Represented by a center, axis, and radius.
 * 
 * Most methods consider the circle as the solid disk.
 * Methods that use the circle boundary have the prefix "boundary".
 */
struct Circle3D
{
  bool setIntersection(const Sphere3D& s,const Plane3D& p);  //false if they don't intersect
  Real distance(const Point3D& v) const;
  Real boundaryDistance(const Point3D& v) const;
  bool intersects(const Circle3D& c) const;
  bool intersects(const Sphere3D& s) const;
  bool intersects(const Line3D& l,Real* t=NULL) const;
  bool intersects(const Plane3D& p) const;
  bool boundaryIntersects(const Sphere3D& s) const;
  void getPlane(Plane3D& p) const;
  void getAABB(AABB3D&) const;

  bool Read(File& f);
  bool Write(File& f) const;

  Point3D center;
  Vector3 axis;
  Real radius;
};

} //namespace Math3D

#endif
