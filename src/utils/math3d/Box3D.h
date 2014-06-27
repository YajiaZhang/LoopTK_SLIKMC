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

#ifndef MATH3D_BOX3D_H
#define MATH3D_BOX3D_H

#include "LocalCoordinates3D.h"

namespace Math3D {

/** @brief A 3D box
 * @ingroup Math3D
 *
 * The box is the unit cube [0,1]^3 set in the scaled local coordinate
 * system.  That is, one corner is at the origin, and it has dimensions
 * [dims.x,dims.y,dims.z] in the coordinates given by {xbasis,ybasis,zbasis}. 
 */
struct Box3D : public ScaledLocalCoordinates3D
{
  void set(const AABB3D& bb);
  bool contains(const Point3D& pt) const;
  bool withinDistance(const Point3D& pt, Real dist) const;
  void getAABB(AABB3D& bb) const;
  bool intersects(const Box3D& b) const;
  bool intersectsApprox(const Box3D& b) const;  ///<faster, approximate version
};

} //namespace Math3D

#endif
