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

#ifndef MATH3D_POLAR_H
#define MATH3D_POLAR_H

#include "primitives.h"

namespace Math3D {

//polar coords: (r,theta) : (x,y) = (r cos theta, r sin theta)
//   r in [0,inf), theta in [0,2pi)
//cylindrical coords: (r,theta,z) : (x,y,z) = (r cos theta, r sin theta, z)
//   r in [0,inf), theta in [0,2pi)
//spherical coords: (r,theta,phi) : (x,y,z) = (r cos theta cos phi, r sin theta cos phi, r sin phi)
//   r in [0,inf), theta in [0,2pi), phi in [-pi, pi]


//conversion routines
void PolarToRectangular(const Vector2& polar, Vector2& rect);
void RectangularToPolar(const Vector2& rect, Vector2& polar);

void SphericalToRectangular(const Vector3& sphere, Vector3& rect);
void RectangularToSpherical(const Vector3& rect, Vector3& sphere);

void CylindricalToRectangular(const Vector3& cyl, Vector3& rect);
void RectangularToCylindrical(const Vector3& rect, Vector3& cyl);

//derivatives of polar coordinates in rectangular coordinates
void PolarDR(const Vector2& polar, Vector2& drect);
void PolarDTheta(const Vector2& polar, Vector2& drect);

void SphericalDR(const Vector3& sphere, Vector3& drect);
void SphericalDTheta(const Vector3& sphere, Vector3& drect);
void SphericalDPhi(const Vector3& sphere, Vector3& drect);

void CylindricalDR(const Vector3& cyl, Vector3& drect);
void CylindricalDTheta(const Vector3& cyl, Vector3& drect);
void CylindricalDZ(const Vector3& cyl, Vector3& drect);

} //namespace Math3D

#endif
