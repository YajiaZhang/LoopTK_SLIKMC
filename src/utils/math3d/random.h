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

#ifndef MATH3D_RANDOM_H
#define MATH3D_RANDOM_H

#include <math/random.h>
#include <math/sample.h>
#include <math/complex.h>
#include "primitives.h"

namespace Math3D {

inline void SampleCircle(Real r,Vector2& v)
{
  Math::SampleCircle(r,v.x,v.y);
}

inline void SampleDisk(Real r,Vector2& v)
{
  Math::SampleDisk(r,v.x,v.y);
}

inline void SampleSphere(Real r,Vector3& v)
{
  Math::SampleSphere(r,v.x,v.y,v.z);
}

inline void SampleBall(Real r,Vector3& v)
{
  Math::SampleBall(r,v.x,v.y,v.z);
}

inline void SampleSquare(Real d,Vector2& v)
{
  v.set(Rand(-d,d),Rand(-d,d));
}

inline void SampleCube(Real d,Vector3& v)
{
  v.set(Rand(-d,d),Rand(-d,d),Rand(-d,d));
}

inline void SampleHyperCube(Real d,Vector4& v)
{
  v.set(Rand(-d,d),Rand(-d,d),Rand(-d,d),Rand(-d,d));
}

void RandRotation(Quaternion& q);

} //namespace Math3D

#endif
