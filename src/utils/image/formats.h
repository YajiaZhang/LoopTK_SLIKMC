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


typedef unsigned int COLOROPTYPE [4];
typedef void (*PIXELGETPROC) (const unsigned char* bits, COLOROPTYPE col);
typedef void (*PIXELSETPROC) (unsigned char* bits, const COLOROPTYPE col);

#include <memory.h>

struct r5g6b5
{
	unsigned int r :5;
	unsigned int g :6;
	unsigned int b :5;
};

struct x1r5g5b5
{
	unsigned int x :1;
	unsigned int r :5;
	unsigned int g :5;
	unsigned int b :5;
};


void argb_get (const unsigned char* bits, COLOROPTYPE col)
{
	col[0] = *bits;
	col[1] = *(bits+1);
	col[2] = *(bits+2);
	col[3] = *(bits+3);
}

void argb_set (unsigned char* bits, const COLOROPTYPE col)
{
	*bits     = col[0];
	*(bits+1) = col[1];
	*(bits+2) = col[2];
	*(bits+3) = col[3];
}

void rgb8_get (const unsigned char* bits, COLOROPTYPE col)
{
	col[0] = *bits;
	col[1] = *(bits+1);
	col[2] = *(bits+2);
}

void rgb8_set (unsigned char* bits, const COLOROPTYPE col)
{
	*bits     = col[0];
	*(bits+1) = col[1];
	*(bits+2) = col[2];
}

void r5g6b5_get (const unsigned char* bits, COLOROPTYPE col)
{
	r5g6b5 b;
	memcpy(&b, bits, sizeof(short));
	col[0] = b.r;
	col[1] = b.g;
	col[2] = b.b;
}

void r5g6b5_set (unsigned char* bits, const COLOROPTYPE col)
{
	r5g6b5 c;
	c.r = col[0];
	c.g = col[1];
	c.b = col[2];
	memcpy(bits, &c, sizeof(short));
}

void x1r5g5b5_get (const unsigned char* bits, COLOROPTYPE col)
{
	x1r5g5b5 c;
	memcpy(&c, bits, sizeof(short));
	col[0] = c.r;
	col[1] = c.g;
	col[2] = c.b;
}

void x1r5g5b5_set (unsigned char* bits, const COLOROPTYPE col)
{
	x1r5g5b5 c;
	c.r = col[0];
	c.g = col[1];
	c.b = col[2];
	memcpy(bits, &c, sizeof(short));
}

void a8_get (const unsigned char* bits, COLOROPTYPE col)
{
	col[0] = *bits;
}

void a8_set (unsigned char* bits, const COLOROPTYPE col)
{
	*bits     = col[0];
}

inline void zero_color(COLOROPTYPE x)
{
	x[0] = x[1] = x[2] = x[3] = 0;
}

inline void add_color(COLOROPTYPE x, const COLOROPTYPE a)
{
	x[0]+=a[0];
	x[1]+=a[1];
	x[2]+=a[2];
	x[3]+=a[3];
}

inline void shift_right_color(COLOROPTYPE x, unsigned int shift)
{
	x[0] = x[0] >> shift;
	x[1] = x[1] >> shift;
	x[2] = x[2] >> shift;
	x[3] = x[3] >> shift;
}
