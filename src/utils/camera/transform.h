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

#ifndef CAMERA_TRANSFORM_H
#define CAMERA_TRANSFORM_H

#include <math3d/primitives.h>
using namespace Math3D;

//#warning "transform.h has been superceded by camera.h"

void SetOrbitTransform(const Vector3& rot, const Vector3& target, float dist, RigidTransform& xform);
void SetFreeTransform(const Vector3& pos, const Vector3& rot, RigidTransform& xform);
void SetTargetTransform(const Vector3& pos, const Vector3& tgt, const Vector3& up, RigidTransform& xform);

void GetOrbitTransform(const RigidTransform& xform, Vector3& rot, Vector3& target, float dist = One);
void GetFreeTransform(const RigidTransform& xform, Vector3& pos, Vector3& rot);
void GetTargetTransform(const RigidTransform& xform, Vector3& pos, Vector3& tgt, Vector3& up, float tgtDist = One);

#endif
