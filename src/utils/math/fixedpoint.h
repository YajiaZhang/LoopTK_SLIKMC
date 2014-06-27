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

#ifndef MATH_FIXED_POINT_H
#define MATH_FIXED_POINT_H

typedef unsigned char Fixed8;
typedef unsigned short Fixed16;
typedef unsigned int Fixed32;

#include <windows.h>

namespace Fixed_24_8
{
	static Fixed32 One = (1 << 8);
	inline Fixed32 FromInt(unsigned int a) { return a << 8; }
	inline Fixed32 FromDouble(double a) { return (Fixed32)(a*256.0 + 0.5); }
	inline unsigned int ToInt(Fixed32 a) { return a >> 8; }
	inline double ToDouble(Fixed32 a) { return double(a)/256.0; }
	inline Fixed32 Add(Fixed32 a, Fixed32 b) { return a+b; }
	inline Fixed32 Sub(Fixed32 a, Fixed32 b) { return a-b; }
	inline Fixed32 Mul(Fixed32 a, Fixed32 b) { return (a>>4)*(b>>4);
		return MulDiv(a,b,One);
	}
	inline Fixed32 Div(Fixed32 a, Fixed32 b) { 
		return MulDiv(One,a,b);
	}
};

#endif
