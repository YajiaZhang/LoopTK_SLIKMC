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

#ifndef MATH3D_LOCAL_COORDINATES3D_H
#define MATH3D_LOCAL_COORDINATES3D_H

#include "Segment3D.h"
#include "Plane3D.h"

namespace Math3D {

//a rigid-body transformation
struct LocalCoordinates3D
{
	void getBasis (Matrix4&) const;
	void getBasisInv (Matrix4&) const;
	void toLocal (const Point3D&, Point3D&) const;
	void toLocalReorient (const Vector3&, Vector3&) const;
	void fromLocal (const Point3D&, Point3D&) const;
	void fromLocalReorient(const Vector3&, Vector3&) const;

	void toLocal (const Line3D&, Line3D&) const;
	void toLocal (const Segment3D&, Segment3D&) const;
	void toLocal (const Plane3D&, Plane3D&) const;
	void fromLocal (const Line3D&, Line3D&) const;
	void fromLocal (const Segment3D&, Segment3D&) const;
	void fromLocal (const Plane3D&, Plane3D&) const;

	bool Read(File& f);
	bool Write(File& f) const;

	Vector3 origin;
	Vector3 xbasis, ybasis, zbasis;	//orthonormal vectors representing the orientation
};

//a rigid-body transformation with a scale in each basis direction
//normalize -> from scaled coordinates to unscaled
//denormalize -> from unscaled coordinates to scaled
struct ScaledLocalCoordinates3D : public LocalCoordinates3D
{
	void getBasisScaled (Matrix4&) const;
	void getBasisScaledInv (Matrix4&) const;
	void normalize (const Vector3&, Vector3&) const;
	void denormalize (const Vector3&, Vector3&) const;

	void toLocalNormalized (const Point3D&, Point3D&) const;
	void toLocalNormalized (const Line3D&, Line3D&) const;
	void toLocalNormalized (const Segment3D&, Segment3D&) const;
	void toLocalNormalized (const Plane3D&, Plane3D&) const;
	void fromLocalNormalized (const Point3D&, Point3D&) const;
	void fromLocalNormalized (const Line3D&, Line3D&) const;
	void fromLocalNormalized (const Segment3D&, Segment3D&) const;
	void fromLocalNormalized (const Plane3D&, Plane3D&) const;

	bool Read(File& f);
	bool Write(File& f) const;

	Vector3 dims;
};

} //namespace Math3D

#endif
