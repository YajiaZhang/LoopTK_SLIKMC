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

#include "geometry3d.h"
#include "interpolate.h"
#include "misc.h"
#include <math/misc.h>
using namespace std;

namespace Math3D {


void LocalCoordinates3D::getBasis(Matrix4& basis) const
{
	basis.set(xbasis,ybasis,zbasis,origin);
}

void LocalCoordinates3D::getBasisInv(Matrix4& basis) const
{
	//transpose basis
	basis.setIdentity();
	basis.setRow1(xbasis);
	basis.setRow2(ybasis);
	basis.setRow3(zbasis);
	Vector3 v;
	basis.mulVector(origin, v);
	basis.setCol4(-v);
}

void LocalCoordinates3D::toLocalReorient(const Vector3& vin, Vector3& vout) const
{
	vout.x=dot(vin,xbasis);
	vout.y=dot(vin,ybasis);
	vout.z=dot(vin,zbasis);
}

void LocalCoordinates3D::toLocal(const Vector3& vin, Vector3& vout) const
{
	toLocalReorient(vin - origin, vout);	
}

void LocalCoordinates3D::fromLocalReorient(const Vector3& vin, Vector3& vout) const
{
	vout = vin.x*xbasis + vin.y*ybasis + vin.z*zbasis;
}

void LocalCoordinates3D::fromLocal(const Vector3& vin, Vector3& vout) const
{
	fromLocalReorient(vin, vout);	
	vout = vout+origin;
}

void LocalCoordinates3D::toLocal(const Line3D& l, Line3D& out) const
{
	toLocalReorient(l.direction, out.direction);
	toLocal(l.source, out.source);
}

void LocalCoordinates3D::fromLocal(const Line3D& l, Line3D& out) const
{
	fromLocalReorient(l.direction, out.direction);
	fromLocal(l.source, out.source);
}

void LocalCoordinates3D::toLocal(const Segment3D& l, Segment3D& out) const
{
	toLocal(l.A, out.A);
	toLocal(l.B, out.B);
}

void LocalCoordinates3D::fromLocal(const Segment3D& l, Segment3D& out) const
{
	fromLocal(l.A, out.A);
	fromLocal(l.B, out.B);
}

void LocalCoordinates3D::toLocal(const Plane3D& p, Plane3D& out) const
{
	toLocalReorient(p.normal, out.normal);
	Vector3 v = p.normal*p.offset;
	Vector3 v_out;
	toLocal(v, v_out);
	out.offset = dot(v_out, out.normal);
}

void LocalCoordinates3D::fromLocal(const Plane3D& p, Plane3D& out) const
{
	fromLocalReorient(p.normal, out.normal);

	Vector3 v = p.normal*p.offset;
	Vector3 v_out;
	fromLocal(v, v_out);
	out.offset = dot(v_out, out.normal);
}


void ScaleXBasis(Matrix4& xform, Real scale)
{
	xform(0,0) *= scale;
	xform(1,0) *= scale;
	xform(2,0) *= scale;
}

void ScaleYBasis(Matrix4& xform, Real scale)
{
	xform(0,1) *= scale;
	xform(1,1) *= scale;
	xform(2,1) *= scale;
}

void ScaleZBasis(Matrix4& xform, Real scale)
{
	xform(0,2) *= scale;
	xform(1,2) *= scale;
	xform(2,2) *= scale;
}

void ScaledLocalCoordinates3D::getBasisScaled(Matrix4& basis) const
{
	basis.set(xbasis*dims.x, ybasis*dims.y, zbasis*dims.z, origin);
}

void ScaledLocalCoordinates3D::getBasisScaledInv(Matrix4& basis) const
{
	//transpose basis
	basis.setIdentity();
	basis.setRow1(xbasis);
	basis.setRow2(ybasis);
	basis.setRow3(zbasis);

	Vector3 v;
	basis.mulVector(origin, v);
	basis.setCol4(-v);

	ScaleXBasis(basis,Inv(dims.x));
	ScaleYBasis(basis,Inv(dims.y));
	ScaleZBasis(basis,Inv(dims.z));
}

void ScaledLocalCoordinates3D::normalize(const Vector3& vin, Vector3& vout) const
{
	vout.x=vin.x/dims.x;
	vout.y=vin.y/dims.y;
	vout.z=vin.z/dims.z;
}

void ScaledLocalCoordinates3D::denormalize(const Vector3& vin, Vector3& vout) const
{
	vout.x=vin.x*dims.x;
	vout.y=vin.y*dims.y;
	vout.z=vin.z*dims.z;
}

void ScaledLocalCoordinates3D::toLocalNormalized(const Point3D& p, Point3D& out) const
{
	toLocal(p, out);
	normalize(out, out);
}

void ScaledLocalCoordinates3D::fromLocalNormalized(const Point3D& p, Point3D& out) const
{
	denormalize(out, out);
	fromLocal(p, out);
}

void ScaledLocalCoordinates3D::toLocalNormalized(const Line3D& l, Line3D& out) const
{
	toLocal(l, out);
	normalize(out.direction, out.direction);
	normalize(out.source, out.source);
}

void ScaledLocalCoordinates3D::fromLocalNormalized(const Line3D& l, Line3D& out) const
{
	Line3D temp;
	denormalize(l.direction, temp.direction);
	denormalize(l.source, temp.source);
	fromLocal(temp, out);
}

void ScaledLocalCoordinates3D::toLocalNormalized(const Segment3D& l, Segment3D& out) const
{
	toLocal(l, out);
	normalize(out.A, out.A);
	normalize(out.B, out.B);
}

void ScaledLocalCoordinates3D::fromLocalNormalized(const Segment3D& l, Segment3D& out) const
{
	Segment3D temp;
	denormalize(l.A, temp.A);
	denormalize(l.B, temp.B);
	fromLocal(temp, out);
}

void ScaledLocalCoordinates3D::toLocalNormalized(const Plane3D& p, Plane3D& out) const
{
	toLocalReorient(p.normal, out.normal);
	denormalize(out.normal, out.normal);
	out.normal.inplaceNormalize();

	Vector3 v = p.normal*p.offset;
	Vector3 v_out;
	toLocal(v, v_out);
	normalize(v_out, v_out);
	out.offset = dot(v_out, out.normal);
}

void ScaledLocalCoordinates3D::fromLocalNormalized(const Plane3D& p, Plane3D& out) const
{
	Plane3D p_denorm;
	normalize(p.normal, p_denorm.normal);
	p_denorm.normal.inplaceNormalize();

	Vector3 v = p.normal * p.offset;
	Vector3 v_out;
	denormalize(v, v_out);

	p_denorm.offset = dot(v_out, p_denorm.normal);

	fromLocal(p_denorm, out);
}


void Box3D::set(const AABB3D& bb)
{
  origin = bb.bmin;
  xbasis.set(1,0,0);
  ybasis.set(0,1,0);
  zbasis.set(0,0,1);
  dims = bb.bmax-bb.bmin;
}

bool Box3D::contains(const Point3D& pt) const
{
	Point3D out;
	toLocal(pt,out);
	return 0<=out.x&&out.x<=dims.x &&
		0<=out.y&&out.y<=dims.y &&
		0<=out.z&&out.z<=dims.z;
}

bool Box3D::withinDistance(const Point3D& pt, Real dist) const
{
	Sphere3D sphereLocal;
	toLocal(pt, sphereLocal.center);
	sphereLocal.radius=dist;
	AABB3D bb_local;
	bb_local.bmin.setZero();
	bb_local.bmax = dims;
	bb_local.justify();
	return sphereLocal.intersects(bb_local);
}

bool Box3D::intersects(const Box3D& b) const
{
  cout<<"Not quite done... check split planes a's faces, b's faces, and a's edges x b's edges"<<endl;
  abort();
  return false;
}

bool Box3D::intersectsApprox(const Box3D& b) const
{
	Box3D temp;
	AABB3D aabb_temp, aabb_temp2;
	//make temp localized
	temp.dims = b.dims;
	toLocal(b.origin, temp.origin);
	toLocalReorient(b.xbasis, temp.xbasis);
	toLocalReorient(b.ybasis, temp.ybasis);
	toLocalReorient(b.zbasis, temp.zbasis);
	temp.getAABB(aabb_temp);
	aabb_temp2.bmin.setZero();
	aabb_temp2.bmax = dims;
	if(!aabb_temp2.intersects(aabb_temp))
		return false;

	temp.dims = dims;
	b.toLocal(origin, temp.origin);
	b.toLocalReorient(xbasis, temp.xbasis);
	b.toLocalReorient(ybasis, temp.ybasis);
	b.toLocalReorient(zbasis, temp.zbasis);
	temp.getAABB(aabb_temp);
	aabb_temp2.bmax = b.dims;
	if(!aabb_temp2.intersects(aabb_temp))
		return false;
	return true;
}

void Box3D::getAABB(AABB3D& bb) const
{
	Vector3 x(dims.x*xbasis),y(dims.y*ybasis),z(dims.z*zbasis);

	Vector3 tmin, tmax;
	tmin.setZero();
	tmax.setZero();
	for(int i=0; i<3; i++)
	{
		if(x[i] > Zero) tmax[i] = x[i];
		else tmin[i] = x[i];
		if(y[i] > Zero) tmax[i] += y[i];
		else tmin[i] += y[i];
		if(z[i] > Zero) tmax[i] += z[i];
		else tmin[i] += z[i];
	}

	bb.bmin.add(tmin, origin);
	bb.bmax.add(tmax, origin);
}


bool Sphere3D::Read(File& f)
{
	if(!center.Read(f)) return false;
	if(!ReadFile(f,radius)) return false;
	return true;
}

bool Sphere3D::Write(File& f) const
{
	if(!center.Write(f)) return false;
	if(!WriteFile(f,radius)) return false;
	return true;
}

Real Sphere3D::distance(const Point3D& v) const
{
	return (center-v).norm() - radius;
}

bool Sphere3D::contains(const Point3D& v) const
{
	return DistanceLEQ(center,v,radius);
}

bool Sphere3D::contains(const Sphere3D& s) const
{
	return DistanceLEQ(center,s.center,radius-s.radius);
}

bool Sphere3D::withinDistance(const Point3D& v, Real dist) const
{
	return DistanceLEQ(center,v, radius + dist);
}

bool Sphere3D::boundaryWithinDistance(const Point3D& v, Real dist) const
{
    return Abs((center-v).norm()-radius) <= dist;
}

bool Sphere3D::intersects(const Line3D& l, Real* t1, Real* t2) const
{
	Vector3 offset=center-l.source;
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
	if(res==1) {
	  cout<<"Whoa, line just intersects at one point on the sphere"<<endl;
	  cout<<"l= "<<l.source<<"->"<<l.direction<<endl;
	  cout<<"c= "<<center<<", r="<<radius<<endl;
	  cout<<"t="<<x1<<endl;
	  getchar();
	  x2=x1;
	}
	if(x1 > x2) Swap(x1,x2);
	if(t1 && t2) {
		*t1=x1;
		*t2=x2;
	}
	return true;
}

bool Sphere3D::intersects(const Plane3D& p) const
{
  return Abs(p.distance(center)) <= radius;
}

bool Sphere3D::intersects(const Sphere3D& s) const
{
  return ballsIntersect(center,radius,s.center,s.radius);
}

bool Sphere3D::boundaryIntersects(const Sphere3D& s) const
{
  return ballSphereIntersect(s.center,s.radius,center,radius);
}

bool Sphere3D::boundaryIntersectsBoundary(const Sphere3D& s) const
{
  return spheresIntersect(center,radius,s.center,s.radius);
}

bool Sphere3D::ballsIntersect(const Point3D& ca,Real ra,const Point3D& cb,Real rb)
{
	return DistanceLEQ(ca,cb,ra+rb);
}

bool Sphere3D::ballSphereIntersect(const Point3D& ca,Real ra,const Point3D& cb,Real rb)
{
  Real r2 = (ca-cb).normSquared();
  if(r2 <= Sqr(ra+rb)) {
    Real r = Sqrt(r2);
    if(r + ra < rb) return false;
    return true;
  }
  return false;
}

bool Sphere3D::spheresIntersect(const Point3D& ca,Real ra,const Point3D& cb,Real rb)
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

void Sphere3D::getAABB(AABB3D& bb) const
{
  bb.setPoint(center);
  bb.bmin.x-=radius; bb.bmin.y-=radius; bb.bmin.z-=radius;
  bb.bmax.x+=radius; bb.bmax.y+=radius; bb.bmax.z+=radius;
}

bool Sphere3D::intersects(const AABB3D& bb) const
{
//this one's interesting
//quick reject- center outside of a box expanded by r
//check the faces, edges, and vertices for intersection
//faces: check if the sphere center lies within a (2*r) height box  (actually, group opposite faces together)
//edges: check if the center lies within an r radius cylinder (like segment)
//vertices: within a the sphere

  const static Real InvSqrt3 = One/Sqrt(3.0);
  //trivial reject
  if(!bb.withinDistance(center,radius)) return false;
  //trivial accept
  if(bb.withinDistance(center,radius*InvSqrt3)) return true;

  int rgnx,rgny,rgnz;
  if(center.x>=bb.bmin.x) {
    if(center.x<=bb.bmax.x) rgnx=0;
    else rgnx=1;
  }
  else rgnx=-1;
  if(center.y>=bb.bmin.y) {
    if(center.y<=bb.bmax.y) rgny=0;
    else rgny=1;
  }
  else rgny=-1;
  if(center.z>=bb.bmin.z) {
    if(center.z<=bb.bmax.z) rgnz=0;
    else rgnz=1;
  }
  else rgnz=-1;

  //closest point hits edges/planes
  if(rgnx==0) {
    //check y-z plane
    if(rgny==0 || rgnz==0) return true;
    else {
      //check to hit edges
      return (DistanceSquared2D(center.y,center.z,
	     (rgny<0?bb.bmin.y:bb.bmax.y),(rgnz<0?bb.bmin.z:bb.bmax.z))
	      <=Sqr(radius));
    }
  }
  else if(rgny==0) {
    //check x-z plane
    if(rgnz==0) return true;
    else {
      //check to hit edges
      return (DistanceSquared2D(center.x,center.z,
	     (rgnx<0?bb.bmin.x:bb.bmax.x),(rgnz<0?bb.bmin.z:bb.bmax.z))
	      <=Sqr(radius));
    }
  }
  else if(rgnz==0) {
    //check to hit edges on x,y, plane
    return (DistanceSquared2D(center.x,center.y,
	   (rgnx<0?bb.bmin.x:bb.bmax.x),(rgny<0?bb.bmin.y:bb.bmax.y))
	    <=Sqr(radius));

  }
  else { //check to hit corners
    return (center.distanceSquared(Vector3((rgnx<0?bb.bmin.x:bb.bmax.x),
					     (rgny<0?bb.bmin.y:bb.bmax.y),
					     (rgnz<0?bb.bmin.z:bb.bmax.z)))
	    <=Sqr(radius));
  }
}


bool Sphere3D::intersectsApprox(const AABB3D& bb) const
{
  return bb.withinDistance(center,radius);
}



bool Ellipsoid3D::contains(const Point3D& pt) const
{
	Point3D out;
	toLocalNormalized(pt,out);
	return NormLEQ(out,One);
}

bool Ellipsoid3D::intersects(const Line3D& l, Real* t1, Real* t2) const
{
	Line3D llocal;
	toLocalNormalized(l,llocal);
	Sphere3D s;
	s.center.setZero();
	s.radius = One;
	return s.intersects(llocal,t1,t2);
}

void Ellipsoid3D::getAABB(AABB3D& bb) const
{
	//get the bases of world space in ellipsoid space
	Vector3 xb,yb,zb;
	xb.x = xbasis.x;
	xb.y = ybasis.x;
	xb.z = zbasis.x;
	yb.x = xbasis.y;
	yb.y = ybasis.y;
	yb.z = zbasis.y;
	zb.x = xbasis.z;
	zb.y = ybasis.z;
	zb.z = zbasis.z;

	normalize(xb,xb);
	normalize(yb,yb);
	normalize(zb,zb);

	//now find the points on the sphere with the correct tangent planes
	Vector3 xt,yt,zt;
	xt = cross(yb,zb);
	yt = cross(zb,xb);
	zt = cross(xb,yb);

	//these are the normals, just normalize them
	xt.inplaceNormalize();
	yt.inplaceNormalize();
	zt.inplaceNormalize();

	xb = xbasis * dims.x;
	yb = ybasis * dims.y;
	zb = zbasis * dims.z;

	//aliases
	Vector3& bmin=bb.bmin, &bmax = bb.bmax;

	//take these points back to world coordinates- these will be the min and max points
	bmax.x = bmin.x = xt.x * xb.x + xt.y * yb.x + xt.z * zb.x;
	bmax.y = bmin.y = yt.x * xb.y + yt.y * yb.y + yt.z * zb.y;
	bmax.z = bmin.z = zt.x * xb.z + zt.y * yb.z + zt.z * zb.z;

	if(bmax.x < 0)
		bmax.x = -bmax.x;
	else
		bmin.x = -bmin.x;

	if(bmax.y < 0)
		bmax.y = -bmax.y;
	else
		bmin.y = -bmin.y;

	if(bmax.z < 0)
		bmax.z = -bmax.z;
	else
		bmin.z = -bmin.z;

	bmax += origin;
	bmin += origin;
}



} // namespace Math3D
