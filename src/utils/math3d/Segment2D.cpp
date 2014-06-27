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

#include "Segment2D.h"
#include "clip.h"
#include "misc.h"
#include <math/misc.h>
#include <math/Interval.h>
#include "interpolate.h"
using namespace Math3D;
using namespace std;


bool Segment2D::Read(File& f)
{
	if(!A.Read(f)) return false;
	if(!B.Read(f)) return false;
	return true;
}

bool Segment2D::Write(File& f) const
{
	if(!A.Write(f)) return false;
	if(!B.Write(f)) return false;
	return true;
}

void Segment2D::setTransformed(const Segment2D& s, const Matrix3& xform)
{
	xform.mulPoint(s.A,A);
	xform.mulPoint(s.B,B);
}

void Segment2D::getLine(Line2D& l) const
{
	l.source=A;
	l.direction=B-A;
}

Real Segment2D::closestPointParameter(const Point2D& in) const
{
	Vector2 dir=B-A;
	Real numer = dot(in-A,dir);
	Real denom = dot(dir,dir);
	//t = numer/denom, denom always >= 0
	if(numer <= Zero) return Zero;
	if(numer >= denom) return One;
	return numer/denom;
}

Real Segment2D::closestPoint(const Point2D& in,Point2D& out) const
{
	Real t = closestPointParameter(in);
	eval(t,out);
	return t;
}

Real Segment2D::distance(const Point2D& pt) const
{
  Point2D closest;
  closestPoint(pt,closest);
  return (pt-closest).norm();
}

void Segment2D::eval(Real t, Point2D& out) const
{
	out = A;
	out.madd(B-A,t);
}

Real Segment2D::orientation(const Point2D& a,const Point2D& b,const Point2D& x)
{
  return Orient2D(a,b,x);
}

Real Segment2D::orientation(const Vector2& x) const
{
  return Orient2D(A,B,x);
}

bool Segment2D::intersects(const Segment2D& S) const
{
  return intersects(S.A,S.B);
}

bool Segment2D::intersects(const Vector2& a,const Vector2& b) const
{
  Real u,v;
  u = orientation(A,B,a);
  v = orientation(A,B,b);
  if((u < Zero && v < Zero) || (u > Zero && v > Zero)) return false;
  u = orientation(a,b,A);
  v = orientation(a,b,B);
  if((u < Zero && v < Zero) || (u > Zero && v > Zero)) return false;
  return true;
}

bool Segment2D::intersects(const Segment2D& S,Vector2& p) const
{
  return intersects(S.A,S.B,p);
}

bool Segment2D::intersects(const Vector2& a,const Vector2& b,Vector2& p) const
{
  Matrix2 M;
  Vector2 res,uv;
  M.setCol1(A-B);
  M.setCol2(b-a);
  res = A-a;
  if(FuzzyZero(M.determinant())) {
    //they're parallel
    Vector2 t = A-B;
    Vector2 n; n.setPerpendicular(t);
    Real D = dot(n,A);
    Real d = dot(n,a);
    if(d == D) {
      ClosedInterval U,u;
      U.a = 0;
      U.b = t.normSquared();
      u.a = dot(t,a-B);
      u.b = dot(t,b-B);
      if(U.intersects(u)) {
	if(U.contains(u.a)) p=B+t*u.a/U.b;
	else if(U.contains(u.b)) p=B+t*u.b/U.b;
	else   //u contains U
	  p=B;
	return true;
      }
    }
    return false;
  }
  M.inplaceInverse();
  M.mul(res,uv);
  if(uv.x>=Zero && uv.x<=One &&
     uv.y>=Zero && uv.y<=One) {
    interpolate(A,B,uv.x,p);
    Vector2 temp;
    interpolate(a,b,uv.y,temp);
    if(temp.distance(p) > 1e-3) {
      cout<<"Error: intersection points are too far away "<<endl;
      cout<<A<<" -> "<<B<<endl;
      cout<<a<<" -> "<<b<<endl;
      cout<<"u,v "<<uv<<endl;
      cout<<"inverse basis "<<endl<<M<<endl;
      cout<<"p1,p2"<<p<<", "<<temp<<endl;
      abort();
    }
    return true;
  }
  /*
  if(intersects(a,b)) {
    if(FuzzyZero(uv.x)) { p=A; return true; }
    if(FuzzyEquals(uv.x,One)) { p=B; return true; }
    if(FuzzyZero(uv.y)) { p=a; return true; }
    if(FuzzyEquals(uv.y,One)) { p=b; return true; }
    cout<<"Error! segment is supposed to intersect, but we don't have that in the basis!"<<endl;
    cout<<A<<" -> "<<B<<endl;
    cout<<a<<" -> "<<b<<endl;
    cout<<"u,v "<<uv<<endl;
    cout<<"inverse basis "<<endl<<M<<endl;
    abort();
  }
  */
  return false;
}

void Segment2D::getAABB(AABB2D& bb) const
{
  bb.setPoint(A);
  bb.expand(B);
}

bool Segment2D::intersects(const AABB2D& bb) const
{
  Real u1=0,u2=1;
  return intersects(bb,u1,u2);
}

bool Segment2D::intersects(const AABB2D& bb, Real& u1, Real& u2) const
{
  return ClipLine(A, B-A, bb, u1,u2);
}
