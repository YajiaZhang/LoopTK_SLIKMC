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

#ifndef MATH3D_ELLIPSOID_H
#define MATH3D_ELLIPSOID_H

#include "LocalCoordinates3D.h"
#include "Line3D.h"

namespace Math3D {

/** @brief A 3D ellipsoid
 * @ingroup Math3D
 *
 * The box is the centered unit cube [-1,1]^3 set in the scaled local
 * coordinate system.  That is, the center is at the origin, and its
 * major axes are {xbasis,ybasis,zbasis} with radii {dims.x,dims.y,dims.z}. 
 */
struct Ellipsoid3D : public ScaledLocalCoordinates3D
{
  bool contains(const Point3D& pt) const;
  bool intersects(const Line3D& l, Real* t1=NULL, Real* t2=NULL) const;
  void getAABB(AABB3D& bb) const;
};

} //namespace Math3D

#endif

