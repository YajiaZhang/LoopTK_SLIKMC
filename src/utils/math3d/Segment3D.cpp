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

#include "Segment3D.h"
#include "clip.h"
#include "interpolate.h"
using namespace Math3D;

bool Segment3D::Read(File& f)
{
	if(!A.Read(f)) return false;
	if(!B.Read(f)) return false;
	return true;
}

bool Segment3D::Write(File& f) const
{
	if(!A.Write(f)) return false;
	if(!B.Write(f)) return false;
	return true;
}

void Segment3D::setTransformed(const Segment3D& s, const Matrix4& xform)
{
	xform.mulPoint(s.A,A);
	xform.mulPoint(s.B,B);
}

void Segment3D::getLine(Line3D& l) const
{
	l.source=A;
	l.direction=B-A;
}

void Segment3D::getAABB(AABB3D& bb) const
{
  bb.setPoint(A);
  bb.expand(B);
}

Real Segment3D::closestPointParameter(const Point3D& in) const
{
	Vector3 dir=B-A;
	Real numer = dot(in-A,dir);
	Real denom = dot(dir,dir);
	//t = numer/denom, denom always >= 0
	if(numer <= Zero) return Zero;
	if(numer >= denom) return One;
	return numer/denom;
}

Real Segment3D::closestPoint(const Point3D& in,Point3D& out) const
{
	Real t = closestPointParameter(in);
	eval(t,out);
	return t;
}

void Segment3D::closestPoint(const Segment3D& s, Real& t, Real& u) const
{
  Line3D l1,l2;
  getLine(l1);
  s.getLine(l2);
  l1.closestPoint(l2,t,u);
  t=Clamp(t,Zero,One);
  u=Clamp(u,Zero,One);
}

void Segment3D::eval(Real t, Point3D& out) const
{
  interpolate(A,B,t,out);
}

bool Segment3D::intersects(const AABB3D& bb) const
{
	Real u1=0,u2=1;
	return intersects(bb,u1,u2);
}

bool Segment3D::intersects(const AABB3D& bb, Real& u1, Real& u2) const
{
  return ClipLine(A, B-A, bb, u1,u2);
}

