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

#ifndef MATH3D_QUATERNION_INLINE_H
#define MATH3D_QUATERNION_INLINE_H

#include "vecinline.h"

namespace Math3D {

typedef vec4_t quat_t;

#define quat_zero vec4_zero
#define quat_equal vec4_equal
#define quat_multiply vec4_multiply
#define quat_add vec4_add
#define quat_dot vec4_dot
#define quat_normalize vec4_normalize

void quat_slerp(quat_t x, const quat_t a, const quat_t b, Real t);

}

#endif
