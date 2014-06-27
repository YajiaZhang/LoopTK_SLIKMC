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

#ifndef CAMERA_CLIP_H
#define CAMERA_CLIP_H

// INCLUSION: The object is completely inside the volume.
// (completely visible)
#define INCLUSION -1
// INTERSECT: The object intersects with one of more of the
// volume's planes. (partially visible)
#define INTERSECT 0
// EXCLUSION: The object is completely outside the volume.
// (invisible)
#define EXCLUSION 1

#include <math3d/geometry3d.h>
#include <vector>
using namespace Math3D;

//NOTE: the planes are defined so that normals point outward
class ConvexVolume
{
public:
	int PointOverlap(const Vector3& pt) const;
	int AABBOverlap(const Vector3& bMin, const Vector3& bMax) const;
	int LineOverlap(const Line3D& pt, Real& tmin, Real& tmax) const;

	void Transform(const Matrix4& m);
	void SetTransform(const ConvexVolume&,const Matrix4& m);

	std::vector<Plane3D> planes;
};

#endif
