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

#include "geometry2d.h"
#include "clip.h"
#include "misc.h"
#include <math/misc.h>
#include <math/Interval.h>
#include "LinearlyDependent.h"
#include "interpolate.h"
using namespace std;

namespace Math3D {


bool Plane2D::Read(File& f)
{
	if(!normal.Read(f)) return false;
	if(!ReadFile(f,offset)) return false;
	return true;
}

bool Plane2D::Write(File& f) const
{
	if(!normal.Write(f)) return false;
	if(!WriteFile(f,offset)) return false;
	return true;
}

void Plane2D::setPointNormal(const Point2D& a, const Vector2& n)
{
	normal.setNormalized(n);
	offset = dot(a,normal);
}

void Plane2D::setLine(const Line2D& l)
{
	normal.setPerpendicular(l.direction);
	setPointNormal(l.source,normal);
}

void Plane2D::setPoints(const Point2D& a, const Point2D& b)
{
  Vector2 v;
  v.setPerpendicular(b-a);
  setPointNormal(a,v);
}

void Plane2D::setTransformed(const Plane2D& pin, const Matrix3& xform)
{
	xform.mulVector(pin.normal, normal);

	Vector2 v, v_out;
	v = pin.normal*pin.offset;
	xform.mulPoint(v, v_out);
	setPointNormal(v_out, normal);
}

//returns the orthogonal distance from the plane to the point
Real Plane2D::distance(const Point2D& v) const
{
	return dot(v,normal) - offset;
}

//projects a Vector2 onto this plane
void Plane2D::project(const Point2D& vin, Point2D& vout) const
{
	vout = vin - normal*distance(vin);
}

void Plane2D::getBasis(Vector2& xb) const
{
  xb.setPerpendicular(normal);
}


bool Plane2D::intersectsSegment(const Segment2D& s, Real* t)
{
  Real da = distance(s.A);
  Real db = distance(s.B);
  if(da < Zero) {
    if(db < Zero) return false;
  }
  else {
    if(db > Zero) return false;
  }
  if(t) {
    //d = da + t*(db-da)
    if(da == db) //segment on plane
      *t = Zero;
    else
      *t = da/(da-db);
  }
  return true;
}

bool Plane2D::intersectsLine(const Line2D& l, Real* t)
{
  Real ds = distance(l.source);
  if(dot(normal,l.direction) == Zero) {
    if(t) *t = Inf;
		return (ds == Zero);
  }
	if(t)
		//dot(s + t*d, n) = o   =>   dot(s,n) + t*dot(d,n) = o
		*t = -ds/dot(l.direction,normal);
	return true;
}

bool Plane2D::intersectsRay(const Ray2D& r, Real* t)
{
	//if source on - side, intersects if dir is +
	//if source on + side, intersects if dir is -
	Real src = distance(r.source);
	Real dir = dot(normal,r.direction);
	if(src < Zero) {
		if(dir > Zero) {
			if(*t) *t = -src/dir;
			return true;
		}
	}
	else if(src > Zero) {
		if(dir < Zero) {
			if(*t) *t = -src/dir;
			return true;
		}
	}
	else {
		if(*t) *t = Zero;
		return true;
	}
	return false;
}

bool Plane2D::intersects(const AABB2D& bb) const
{
	Vector2 vmin,vmax;
	//get the extreme points of the box relative to the plane's normal
	if(normal.x > Zero)
	{
		vmin.x = bb.bmin.x;
		vmax.x = bb.bmax.x;
	}
	else
	{
		vmin.x = bb.bmax.x;
		vmax.x = bb.bmin.x;
	}
	if(normal.y > Zero)
	{
		vmin.y = bb.bmin.y;
		vmax.y = bb.bmax.y;
	}
	else
	{
		vmin.y = bb.bmax.y;
		vmax.y = bb.bmin.y;
	}
	//intersects if the extreme points are on opposite sides of the plane
	return (distance(vmin) < Zero) != (distance(vmax) < Zero);

}

int Plane2D::allIntersections(const Plane2D& p,Vector2& pt) const
{
  Real x = p.normal.y*offset - normal.y*p.offset;
  Real y = normal.x*p.offset - p.normal.x*offset;
  Real w = normal.x*p.normal.y - normal.y*p.normal.x;
  if(w < Epsilon) {
    //linearly dependent
    if(x < Epsilon && y < Epsilon)  //overlap
      return 2;
    //pt at infinity
    pt.x = x;
    pt.y = y;
    return 0;
  }
  else {
    pt.x = x/w;
    pt.y = y/w;
    return 1;
  }
}

/*
void LocalCoordinates2D::getBasis(Matrix3& basis) const
{
	basis.set(xbasis,ybasis,zbasis,origin);
}

void LocalCoordinates2D::getBasisInv(Matrix3& basis) const
{
	//transpose basis
	basis.setIdentity();
	basis.setRow1(xbasis);
	basis.setRow2(ybasis);
	basis.setRow2(zbasis);
	Vector2 v;
	basis.mulVector(origin, v);
	basis.setCol4(-v);
}

void LocalCoordinates2D::toLocalReorient(const Vector2& vin, Vector2& vout) const
{
	vout.x=dot(vin,xbasis);
	vout.y=dot(vin,ybasis);
	vout.z=dot(vin,zbasis);
}

void LocalCoordinates2D::toLocal(const Vector2& vin, Vector2& vout) const
{
	toLocalReorient(vin - origin, vout);	
}

void LocalCoordinates2D::fromLocalReorient(const Vector2& vin, Vector2& vout) const
{
	vout = vin.x*xbasis + vin.y*ybasis + vin.z*zbasis;
}

void LocalCoordinates2D::fromLocal(const Vector2& vin, Vector2& vout) const
{
	fromLocalReorient(vin, vout);	
	vout = vout+origin;
}

void LocalCoordinates2D::toLocal(const Line2D& l, Line2D& out) const
{
	toLocalReorient(l.direction, out.direction);
	toLocal(l.source, out.source);
}

void LocalCoordinates2D::fromLocal(const Line2D& l, Line2D& out) const
{
	fromLocalReorient(l.direction, out.direction);
	fromLocal(l.source, out.source);
}

void LocalCoordinates2D::toLocal(const Segment2D& l, Segment2D& out) const
{
	toLocal(l.A, out.A);
	toLocal(l.B, out.B);
}

void LocalCoordinates2D::fromLocal(const Segment2D& l, Segment2D& out) const
{
	fromLocal(l.A, out.A);
	fromLocal(l.B, out.B);
}

void LocalCoordinates2D::toLocal(const Plane2D& p, Plane2D& out) const
{
	toLocalReorient(p.normal, out.normal);
	Vector2 v = p.normal*p.offset;
	Vector2 v_out;
	toLocal(v, v_out);
	out.offset = dot(v_out, out.normal);
}

void LocalCoordinates2D::fromLocal(const Plane2D& p, Plane2D& out) const
{
	fromLocalReorient(p.normal, out.normal);

	Vector2 v = p.normal*p.offset;
	Vector2 v_out;
	fromLocal(v, v_out);
	out.offset = dot(v_out, out.normal);
}


void ScaleXBasis(Matrix3& xform, Real scale)
{
	xform(0,0) *= scale;
	xform(1,0) *= scale;
	xform(2,0) *= scale;
}

void ScaleYBasis(Matrix3& xform, Real scale)
{
	xform(0,1) *= scale;
	xform(1,1) *= scale;
	xform(2,1) *= scale;
}

void ScaleZBasis(Matrix3& xform, Real scale)
{
	xform(0,2) *= scale;
	xform(1,2) *= scale;
	xform(2,2) *= scale;
}

void ScaledLocalCoordinates2D::getBasisScaled(Matrix3& basis) const
{
	basis.set(xbasis*dims.x, ybasis*dims.y, zbasis*dims.z, origin);
}

void ScaledLocalCoordinates2D::getBasisScaledInv(Matrix3& basis) const
{
	//transpose basis
	basis.setIdentity();
	basis.setRow1(xbasis);
	basis.setRow2(ybasis);
	basis.setRow2(zbasis);

	Vector2 v;
	basis.mulVector(origin, v);
	basis.setCol4(-v);

	ScaleXBasis(basis,Inv(dims.x));
	ScaleYBasis(basis,Inv(dims.y));
	ScaleZBasis(basis,Inv(dims.z));
}

void ScaledLocalCoordinates2D::normalize(const Vector2& vin, Vector2& vout) const
{
	vout.x=vin.x/dims.x;
	vout.y=vin.y/dims.y;
	vout.z=vin.z/dims.z;
}

void ScaledLocalCoordinates2D::denormalize(const Vector2& vin, Vector2& vout) const
{
	vout.x=vin.x*dims.x;
	vout.y=vin.y*dims.y;
	vout.z=vin.z*dims.z;
}

void ScaledLocalCoordinates2D::toLocalNormalized(const Point2D& p, Point2D& out) const
{
	toLocal(p, out);
	normalize(out, out);
}

void ScaledLocalCoordinates2D::fromLocalNormalized(const Point2D& p, Point2D& out) const
{
	denormalize(out, out);
	fromLocal(p, out);
}

void ScaledLocalCoordinates2D::toLocalNormalized(const Line2D& l, Line2D& out) const
{
	toLocal(l, out);
	normalize(out.direction, out.direction);
	normalize(out.source, out.source);
}

void ScaledLocalCoordinates2D::fromLocalNormalized(const Line2D& l, Line2D& out) const
{
	Line2D temp;
	denormalize(l.direction, temp.direction);
	denormalize(l.source, temp.source);
	fromLocal(temp, out);
}

void ScaledLocalCoordinates2D::toLocalNormalized(const Segment2D& l, Segment2D& out) const
{
	toLocal(l, out);
	normalize(out.A, out.A);
	normalize(out.B, out.B);
}

void ScaledLocalCoordinates2D::fromLocalNormalized(const Segment2D& l, Segment2D& out) const
{
	Segment2D temp;
	denormalize(l.A, temp.A);
	denormalize(l.B, temp.B);
	fromLocal(temp, out);
}

void ScaledLocalCoordinates2D::toLocalNormalized(const Plane2D& p, Plane2D& out) const
{
	toLocalReorient(p.normal, out.normal);
	denormalize(out.normal, out.normal);
	out.normal.inplaceNormalize();

	Vector2 v = p.normal*p.offset;
	Vector2 v_out;
	toLocal(v, v_out);
	normalize(v_out, v_out);
	out.offset = dot(v_out, out.normal);
}

void ScaledLocalCoordinates2D::fromLocalNormalized(const Plane2D& p, Plane2D& out) const
{
	Plane2D p_denorm;
	normalize(p.normal, p_denorm.normal);
	p_denorm.normal.inplaceNormalize();

	Vector2 v = p.normal * p.offset;
	Vector2 v_out;
	denormalize(v, v_out);

	p_denorm.offset = dot(v_out, p_denorm.normal);

	fromLocal(p_denorm, out);
}

bool Box2D::contains(const Point2D& pt) const
{
	Point2D out;
	toLocal(pt,out);
	return 0<=out.x&&out.x<=dims.x &&
		0<=out.y&&out.y<=dims.y &&
		0<=out.z&&out.z<=dims.z;
}

bool Box2D::withinDistance(const Point2D& pt, Real dist) const
{
	Circle2D c;
	toLocal(pt, c.center);
	c.radius = dist;

	AABB2D bb_local;
	bb_local.bmin.setZero();
	bb_local.bmax=dims;
	return bb.intersects(c);
}

bool Box2D::intersects(const Box2D& b) const
{
	Box2D temp;
	AABB2D AABB2D_temp, AABB2D_temp2;
	//make temp localized
	temp.dims = b.dims;
	toLocal(b.origin, temp.origin);
	toLocalReorient(b.xbasis, temp.xbasis);
	toLocalReorient(b.ybasis, temp.ybasis);
	toLocalReorient(b.zbasis, temp.zbasis);
	AABB2D_temp.calculate(temp);
	AABB2D_temp2.bmin.setZero();
	AABB2D_temp2.bmax = dims;
	if(!AABB2D_temp2.intersects(AABB2D_temp))
		return false;

	temp.dims = dims;
	b.toLocal(origin, temp.origin);
	b.toLocalReorient(xbasis, temp.xbasis);
	b.toLocalReorient(ybasis, temp.ybasis);
	b.toLocalReorient(zbasis, temp.zbasis);
	AABB2D_temp.calculate(temp);
	AABB2D_temp2.bmax = b.dims;
	if(!AABB2D_temp2.intersects(AABB2D_temp))
		return false;
	return true;
}

*/

bool Circle2D::Read(File& f)
{
	if(!center.Read(f)) return false;
	if(!ReadFile(f,radius)) return false;
	return true;
}

bool Circle2D::Write(File& f) const
{
	if(!center.Write(f)) return false;
	if(!WriteFile(f,radius)) return false;
	return true;
}

Real Circle2D::distance(const Point2D& v) const
{
	return (center-v).norm() - radius;
}

bool Circle2D::contains(const Point2D& v) const
{
	return DistanceLEQ(center,v,radius);
}

bool Circle2D::contains(const Circle2D& s) const
{
	return DistanceLEQ(center,s.center,radius-s.radius);
}

bool Circle2D::withinDistance(const Point2D& v, Real dist) const
{
	return DistanceLEQ(center,v, radius + dist);
}

bool Circle2D::boundaryWithinDistance(const Point2D& v, Real dist) const
{
    return Abs((center-v).norm()-radius) <= dist;
}

bool Circle2D::intersects(const Line2D& l, Real* t1, Real* t2) const
{
	Vector2 offset=center-l.source;
	Real o_o=dot(offset,offset), o_b=dot(offset,l.direction), b_b=dot(l.direction,l.direction);
	//so we know there's a root to |offset-t*b|==r
	//o.o-2t*b.o+t^2*b.b=r^2
	Real a,b,c;
	a=b_b;
	b=-Two*o_b;
	c=o_o-radius*radius;
	Real x1,x2;
	int res=quadratic(a,b,c,x1,x2);
	if(res<=0) return false;
	if(t1 && t2) {
		*t1=x1;
		*t2=x2;
	}
	return true;
}

bool Circle2D::intersects(const Plane2D& p, Segment2D& s) const
{
  Real d=p.distance(center);
  if(Abs(d) <= radius) {
    //project c on plane
    Vector2 c=center,basis;
    c.madd(p.normal,-d);
    p.getBasis(basis);
    Real h=pythag_leg(d,radius);
    s.A = c + h*basis;
    s.B = c - h*basis;
    return true;
  }
  return false;
}

bool Circle2D::intersects(const Circle2D& s) const
{
  return disksIntersect(center,radius,s.center,s.radius);
}

bool Circle2D::boundaryIntersects(const Circle2D& s) const
{
  return diskCircleIntersect(s.center,s.radius,center,radius);
}

bool Circle2D::boundaryIntersectsBoundary(const Circle2D& s) const
{
  return circlesIntersect(center,radius,s.center,s.radius);
}

bool Circle2D::disksIntersect(const Point2D& ca,Real ra,const Point2D& cb,Real rb)
{
	return DistanceLEQ(ca,cb,ra+rb);
}

bool Circle2D::diskCircleIntersect(const Point2D& ca,Real ra,const Point2D& cb,Real rb)
{
  Real r2 = (ca-cb).normSquared();
  if(r2 <= Sqr(ra+rb)) {
    Real r = Sqrt(r2);
    if(r + ra < rb) return false;
    return true;
  }
  return false;
}

bool Circle2D::circlesIntersect(const Point2D& ca,Real ra,const Point2D& cb,Real rb)
{
  Real r2 = (ca-cb).normSquared();
  if(r2 <= Sqr(ra+rb)) {
    Real r = Sqrt(r2);
    if(r + ra < rb) return false;
    if(r + rb < ra) return false;
    return true;
  }
  return false;
}

void Circle2D::getAABB(AABB2D& bb) const
{
	bb.setPoint(center);
	bb.bmin.x-=radius; bb.bmin.y-=radius;
	bb.bmax.x+=radius; bb.bmax.y+=radius;
}

bool Circle2D::intersectsApprox(const AABB2D& bb) const
{
  // approximate 
  return (bb.bmin.x <= center.x + radius && center.x - radius <= bb.bmax.x) &&
    (bb.bmin.y <= center.y + radius && center.y - radius <= bb.bmax.y);
}



/*

bool Ellipse2D::contains(const Point2D& pt) const
{
	Point2D out;
	toLocalNormalized(pt,out);
	return NormLEQ(out,One);
}

bool Ellipse2D::intersects(const Line2D& l, Real* t1, Real* t2) const
{
	Line2D llocal;
	toLocalNormalized(l,llocal);
	Circle2D s;
	s.center.setZero();
	s.radius = One;
	return s.intersects(llocal,t1,t2);
}
*/



} // namespace Math2D
