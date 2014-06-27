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

#include "AABB3D.h"
using namespace Math3D;

AABB3D::AABB3D()
{
}

AABB3D::AABB3D(const Vector3& _bmin,const Vector3& _bmax)
:AABBTemplate<Vector3>(_bmin,_bmax)
{
}

AABB3D::AABB3D(const AABB3D& rhs)
:AABBTemplate<Vector3>(rhs)
{
}

bool AABB3D::contains(const Vector3& v) const
{
	return (v.x>=bmin.x && v.x<=bmax.x &&
		v.y>=bmin.y && v.y<=bmax.y &&
		v.z>=bmin.z && v.z<=bmax.z);
}

bool AABB3D::withinDistance(const Vector3& v,Real d) const
{
	return (v.x+d>=bmin.x && v.x-d<=bmax.x &&
		v.y+d>=bmin.y && v.y-d<=bmax.y &&
		v.z+d>=bmin.z && v.z-d<=bmax.z);
}

bool AABB3D::intersects(const AABB3D& a) const
{
	return (bmin.x <= a.bmax.x && bmax.x >= a.bmin.x &&
		bmin.y <= a.bmax.y && bmax.y >= a.bmin.y &&
		bmin.z <= a.bmax.z && bmax.z >= a.bmin.z);
}

void AABB3D::justify()
{
	if(bmin.x > bmax.x) std::swap(bmin.x,bmax.x);
	if(bmin.y > bmax.y) std::swap(bmin.y,bmax.y);
	if(bmin.z > bmax.z) std::swap(bmin.z,bmax.z);
}

void AABB3D::setTransform(const AABB3D& b,const Matrix4& mat)
{
	*this = b;
	inplaceTransform(mat);
}

void AABB3D::inplaceTransform(const Matrix4& mat)
{
	Vector3 dims[8];
	Vector3 dimsTransformed[8];
	dims[0].set(bmin.x,bmin.y,bmin.z);
	dims[1].set(bmin.x,bmin.y,bmax.z);
	dims[2].set(bmin.x,bmax.y,bmin.z);
	dims[3].set(bmin.x,bmax.y,bmax.z);
	dims[4].set(bmax.x,bmin.y,bmin.z);
	dims[5].set(bmax.x,bmin.y,bmax.z);
	dims[6].set(bmax.x,bmax.y,bmin.z);
	dims[7].set(bmax.x,bmax.y,bmax.z);
	for(int i=0;i<8;i++)
		mat.mulPoint(dims[i],dimsTransformed[i]);
	setPoint(dimsTransformed[0]);
	for(int i=1;i<8;i++)
		expand(dimsTransformed[i]);
}



