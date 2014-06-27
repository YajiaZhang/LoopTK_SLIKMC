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

#ifndef MATH3D_LINEARLY_DEPENDENT_H
#define MATH3D_LINEARLY_DEPENDENT_H

#include "primitives.h"

namespace Math3D {

//Robust determination of linear dependence.
//Return true if the vectors are dependent.
//Single-vector versions:
//  Returns a constant c s.t. they're linearly dependent within rel error eps
//    and |c| <= 1.
//  2 cases, (1) a*c = b with maxabs(a*c-b)/|a| <= eps.
//  (2) a = c*b with maxabs(a-c*b)/|b| <= eps.
//  in case 1, cb = false, in case 2, cb=true.
//  c is chosen by pseudoinverse
//Matrix versions:
//  Returns a vector c s.t. A*c = 0.  All coefficients |ci| <= 1
//  At least one ci = 1.

bool LinearlyDependent_Robust(const Vector2& a, const Vector2& b, Real& c, bool& cb, Real eps = Epsilon); 
bool LinearlyDependent_Robust(const Vector3& a, const Vector3& b, Real& c, bool& cb, Real eps = Epsilon); 
bool LinearlyDependent_Robust(const Vector4& a, const Vector4& b, Real& c, bool& cb, Real eps = Epsilon); 

} 

#endif 
