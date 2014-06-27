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

#include "Line2D.h"
#include "clip.h"
#include "misc.h"
#include <math/misc.h>
#include "interpolate.h"

using namespace Math3D;


bool Line2D::Read(File& f)
{
	if(!source.Read(f)) return false;
	if(!direction.Read(f)) return false;
	return true;
}

bool Line2D::Write(File& f) const
{
	if(!source.Write(f)) return false;
	if(!direction.Write(f)) return false;
	return true;
}

void Line2D::setPoints(const Point2D& a, const Point2D& b)
{
	source = a;
	direction.sub(b,a);
}

void Line2D::setSegment(const Segment2D& s)
{
	source = s.A;
	direction.sub(s.B,s.A);
}

void Line2D::setTransformed(const Line2D& l, const Matrix3& xform)
{
	xform.mulPoint(l.source,source);
	xform.mulVector(l.direction,direction);
}

void Line2D::eval(Real t, Point2D& out) const
{
	out = source;
	out.madd(direction,t);
}

Real Line2D::closestPointParameter(const Point2D& in) const
{
	Real denom = dot(direction,direction);
	if(denom == Zero) return Zero;
	return dot(in-source,direction)/denom;
}

Real Line2D::closestPoint(const Point2D& in, Point2D& out) const
{
	Real t=closestPointParameter(in);
	eval(t,out);
	return t;
}

Real Line2D::closestPoint(const Point2D& in, Point2D& out, Real tmin, Real tmax) const
{
	Real denom = dot(direction,direction);
	Real numer = dot(in-source,direction);
	//t = numer/denom with denom >= 0
	Real t;
	if(numer<=tmin*denom) t=tmin;
	else if(numer>=tmax*denom) t=tmax;
	else t = numer/denom;
	eval(t,out);
	return t;
}

Real Line2D::distance(const Point2D& p) const
{
  Vector2 n;
  n.setPerpendicular(direction);
  n.inplaceNormalize();
  return Abs(dot(n,p-source));
}

Real Line2D::orientation(const Vector2& x) const
{
  return cross(direction,x-source);
}


void Line2D::getAABB(AABB2D& bb, Real tmin, Real tmax) const
{
	Point2D a,b;
	eval(tmin,a);
	eval(tmax,b);
	bb.setPoint(a);
	bb.expand(b);
}

bool Line2D::lineIntersects(const AABB2D& bb) const
{
	Real u1=-Inf,u2=Inf;
	return intersects(bb,u1,u2);
}

bool Line2D::rayIntersects(const AABB2D& bb) const
{
	Real u1=0,u2=Inf;
	return intersects(bb,u1,u2);
}

bool Line2D::intersects(const AABB2D& bb, Real& u1, Real& u2) const
{
  return ClipLine(source, direction, bb, u1,u2);
}

