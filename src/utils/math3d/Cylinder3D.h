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

#ifndef MATH3D_CYLINDER3D_H
#define MATH3D_CYLINDER3D_H

#include "Circle3D.h"

namespace Math3D {

/** @ingroup Math3D
 * @brief A 3D cylinder.
 *
 * The base is centered at center, extruding a circle of the
 * given radius along the axis, for the given height.
 */
struct Cylinder3D
{
  void getBase(Circle3D& c) const;
  void getCap(Circle3D& c) const;
  void getAABB(AABB3D&) const;

  /// Returns the entry/exit points if they intersect
  bool intersects(const Line3D& line,Real* tmin=NULL,Real* tmax=NULL) const;

  Point3D center;
  Vector3 axis;
  Real radius;
  Real height;
};

} //namespace Math3D

#endif
