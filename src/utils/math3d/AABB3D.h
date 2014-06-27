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

#ifndef MATH3D_AABB3D_H
#define MATH3D_AABB3D_H

#include "AABBTemplate.h"
#include "primitives.h"

namespace Math3D {

typedef Vector3 Point3D;

/** @brief A 3D axis-aligned bounding box
 * @ingroup Math3D
 */
struct AABB3D : public AABBTemplate<Vector3>
{
  AABB3D();
  AABB3D(const Vector3& bmin,const Vector3& bmax);
  AABB3D(const AABB3D&);
  void justify();  ///<swaps negative sized entries (where min<max)
  void setTransform(const AABB3D&,const Matrix4& mat);
  void inplaceTransform(const Matrix4& mat);
  
  bool contains(const Point3D&) const;
  bool withinDistance(const Point3D&,Real d) const;
  bool intersects(const AABB3D&) const;
};

} //namespace Math3D

#endif

