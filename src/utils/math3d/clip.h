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

#ifndef MATH3D_CLIP_H
#define MATH3D_CLIP_H

#include "Plane2D.h"
#include "Polygon2D.h"
#include "AABB2D.h"
#include "Plane3D.h"
#include "Polyhedron3D.h"
#include "AABB3D.h"
#include "Box3D.h"

/** @file math3d/clip.h
 * @ingroup Math3D
 * @brief Functions for clipping lines against geometric primitives
 */

namespace Math3D {
  /** @addtogroup Math3D */
  /*@{*/

using namespace Math;

/** @brief Clipping primitive in 1D. 
 *
 * Given line segment x = q+p*u, u in [umin,umax],
 * clips the range to the constraint x <= 0.
 * 
 * Returns true the resulting range is non-empty.
 */
bool ClipLine1D(Real q, Real p, Real& umin, Real& umax);

/// Given line segment x = x0+u, u in [umin,umax],
/// clips the range to the constraint a*x <= b
inline bool ClipLine1D(Real x0, Real a, Real b, Real& umin, Real& umax)
{
  return ClipLine1D(a*x0-b,a,umin,umax);
}

///Clip a line (x,v) with range [u1,u2] to the plane's negative halfspace
///i.e. the normal points outward.
bool ClipLine(const Vector2& x, const Vector2& v, const Plane2D& b, Real& u1, Real& u2);
///Clip the line to the interior of the bbox
bool ClipLine(const Vector2& x, const Vector2& v, const AABB2D& b, Real& u1, Real& u2);
///Clip the line to the interior of the convex polygon
bool ClipLine(const Vector2& x, const Vector2& v, const ConvexPolygon2D& b, Real& u1, Real& u2);
bool ClipLine(const Vector3& x, const Vector3& v, const Plane3D& b, Real& u1, Real& u2);
bool ClipLine(const Vector3& x, const Vector3& v, const Box3D& b, Real& u1, Real& u2);
bool ClipLine(const Vector3& x, const Vector3& v, const AABB3D& b, Real& u1, Real& u2);
bool ClipLine(const Vector3& x, const Vector3& v, const ConvexPolyhedron3D& b, Real& u1, Real& u2);

  /*@}*/
} 

#endif
