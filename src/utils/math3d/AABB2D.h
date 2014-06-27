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

#ifndef MATH3D_AABB2D_H
#define MATH3D_AABB2D_H

#include "primitives.h"
#include "AABBTemplate.h"

namespace Math3D {

typedef Vector2 Point2D;

/** @ingroup Math3D
 * @brief A 2D axis-aligned bounding box
 */
struct AABB2D : public AABBTemplate<Vector2>
{
	AABB2D();
	AABB2D(const Vector2& bmin,const Vector2& bmax);
	AABB2D(const AABB2D&);
	void justify();  ///<swaps negative sized entries (where min<max)
	void setTransform(const AABB2D&,const Matrix3& mat);
	void inplaceTransform(const Matrix3& mat);

	bool contains(const Point2D&) const;
	bool intersects(const AABB2D&) const;
};

} //namespace Math3D

#endif
