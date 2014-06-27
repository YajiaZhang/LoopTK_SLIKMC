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

#include "Polygon2D.h"
#include <math/Interval.h>
#include "misc.h"
#include "clip.h"
#include "interpolate.h"
#include <errors.h>
using namespace Math3D;
using namespace std;

bool Polygon2D::ccw() const
{
  size_t i=0,j=1,k=2;
  if(vertices.size() < 3) return true;
  while(i < vertices.size()) {
    if(Orient2D(vertices[i],vertices[j],vertices[k]) < Zero)
      return false;
    i=j; j=k; k++;
    if(k >= vertices.size()) k=0;
  }
  return true;
}

bool Polygon2D::nonCrossing() const
{
  Segment2D ei,ej;
  //faster method would compare monotonic pieces
  for(size_t i=0;i<vertices.size();i++) {
    getEdge(i,ei);
    for(size_t j=0;j+1<i;j++) {  //don't compare adjacent segments
      getEdge(i,ej);
      if(ei.intersects(ej)) return false;
    }
  }
  return true;
}

void Polygon2D::getEdge(int i,Segment2D& ei) const
{
  size_t j=next(i);
  ei.A=vertices[i];
  ei.B=vertices[j];
}

void Polygon2D::getPlane(int i,Plane2D& pi) const
{
  size_t j=next(i);
  //must make planes point outside, normal order points inward
  pi.setPoints(vertices[j],vertices[i]);
}

void Polygon2D::setTransformed(const Polygon2D& in, const Matrix3& T)
{
	vertices.resize(in.vertices.size());
	for(size_t i=0; i<vertices.size(); i++)
		T.mulPoint(in.vertices[i], vertices[i]);
}

bool Polygon2D::planeSplits(const Plane2D& p) const
{
	ClosedInterval x;	x.setEmpty();
	for(size_t i=0; i<vertices.size(); i++) {
		x.expand(p.distance(vertices[i]));
		if(x.contains(Zero))
			return true;
	}
	return false;
}

bool Polygon2D::planePos(const Plane2D& p) const
{
	for(size_t i=0; i<vertices.size(); i++)	{
		if(p.distance(vertices[i]) < Zero)
			return false;
	}
	return true;
}

bool Polygon2D::planeNeg(const Plane2D& p) const
{
	for(size_t i=0; i<vertices.size(); i++) {
		if(p.distance(vertices[i]) > Zero)
			return false;
	}
	return true;
}

bool Polygon2D::raySplits(const Vector2& a,const Vector2& b) const
{
	ClosedInterval x;	x.setEmpty();
	for(size_t i=0; i<vertices.size(); i++) {
		x.expand(Orient2D(a,b,vertices[i]));
		if(x.contains(Zero))
			return true;
	}
	return false;
}

bool Polygon2D::rayLeft(const Vector2& a,const Vector2& b) const
{
	for(size_t i=0; i<vertices.size(); i++) {
		if(Orient2D(a,b,vertices[i]) < Zero)
			return false;
	}
	return true;
}

bool Polygon2D::rayRight(const Vector2& a,const Vector2& b) const
{
	for(size_t i=0; i<vertices.size(); i++) {
		if(Orient2D(a,b,vertices[i]) > Zero)
			return false;
	}
	return true;
}

int Polygon2D::residue(const Vector2& x) const
{
  //do the residue with a vertical ray
  int res=0;
  for(size_t i=0;i<vertices.size();i++) {
    const Vector2& vi=vertices[i];
    const Vector2& vj=vertices[next(i)];
    if(vi.x < x.x) {
      if(vj.x > x.x) { //crosses x=x axis left to right
	//does it cross the axis above x?
	Real x1=vi.x-x.x;
	Real x2=vj.x-x.x;
	Real y1=vi.y-x.y;
	Real y2=vj.y-x.y;
	//Real y=(-x2*y1+x1*y2)/(vi.x-vj.x);
	//if(y >= Zero) res--;
	//if(Sign(-x2*y1+x1*y2) == Sign(vi.x-vj.x) or 0) res--;
	if(Sign(-x2*y1+x1*y2) != 1) res--;
      }
      else if(vj.x == x.x) { //ends up on x.x
	//count as if it's x+eps, which means it doesn't count at all
      }
    }
    else if(vi.x > x.x) {
      if(vj.x < x.x) { //crosses x=x axis right to left
	//does it cross the axis above x?
	Real x1=vi.x-x.x;
	Real x2=vj.x-x.x;
	Real y1=vi.y-x.y;
	Real y2=vj.y-x.y;
	if(Sign(-x2*y1+x1*y2) != -1) res++;
      }
      else if(vj.x == x.x) {  //ends up on x.x
	//count it if it crosses the axis x=x
	if(vj.y >= x.y) res++;
      }
    }
    else if(vi.x == x.x) { //starts out on x.x
      //count only if it goes to the right
      if(vj.x > x.x) {
	if(vi.y > x.y) res--;
      }
    }
  }
  return res;
}

//bool Polygon2D::intersectsBoundary(const Polygon2D& other) const;

Real Polygon2D::boundaryDistance(const Point2D& v) const
{
  if(vertices.size() == 0) return 0;
  else if(vertices.size() == 1) return (v-vertices[0]).norm();
  else {
    Segment2D S; S.A=vertices[0]; S.B=vertices[1];
    return S.distance(v);
  }
  Segment2D s;
  getEdge(0,s);
  Real dmax=s.distance(v);
  for(size_t i=1; i<vertices.size(); i++) {
    getEdge(i,s);
    dmax = Max(dmax,s.distance(v));
  }
  return dmax;
}

bool Polygon2D::intersects(const Line2D& l, Real& tmin, Real& tmax) const
{
  abort();
  return false;
  Segment2D s;
  for(size_t i=0;i<vertices.size();i++) {
    getEdge(i,s);
  }
}

bool Polygon2D::intersects(const Line2D& l) const
{
  Real tmin=-Inf,tmax=Inf;
  return intersects(l,tmin,tmax);
}

bool Polygon2D::intersects(const Segment2D& l) const
{
  for(size_t i=0;i<vertices.size();i++) {
    if(l.intersects(vertices[i],vertices[next(i)])) return true;
  }
  return false;
}

void Polygon2D::getAABB(AABB2D& bb) const
{
  if(vertices.size() == 0) {
    bb.minimize();
    return;
  }
  bb.setPoint(vertices[0]);
  for(size_t i=1; i<vertices.size(); i++)
    bb.expand(vertices[i]);
}







bool ConvexPolygon2D::isValid() const
{
  return ccw();
}

bool ConvexPolygon2D::intersects(const ConvexPolygon2D& other) const
{
  size_t i;
  for(i=0; i<other.vertices.size(); i++) {
    if(rayRight(other.vertices[i],other.vertices[next(i)])) return false;
  }
  for(i=0; i<vertices.size(); i++) {
    if(other.rayRight(vertices[i],vertices[next(i)])) return false;
  }
  return true;
}


bool ConvexPolygon2D::contains(const Point2D& v) const
{
  for(size_t i=0; i<vertices.size(); i++)
    if(Orient2D(vertices[i],vertices[next(i)],v) < Zero) return false;
  return true;
}

bool ConvexPolygon2D::withinEdgeDistance(const Point2D& v, Real dist) const
{
  Plane2D p;
  for(size_t i=0; i<vertices.size(); i++) {
    getPlane((int)i,p);
    if(p.distance(v) > dist) return false;
  }
  return true;
}

Real ConvexPolygon2D::edgeDistance(const Point2D& v) const
{
  if(vertices.size() == 0) return 0;
  else if(vertices.size() == 1) return (v-vertices[0]).norm();
  else {
    Segment2D S; S.A=vertices[0]; S.B=vertices[1];
    return S.distance(v);
  }
  Plane2D p;
  getPlane(0,p);
  Real dmax=p.distance(v);
  for(size_t i=1; i<vertices.size(); i++) {
    getPlane(i,p);
    dmax = Max(dmax,p.distance(v));
  }
  return dmax;
}

bool ConvexPolygon2D::intersects(const Line2D& l, Real& tmin, Real& tmax) const
{
	return ClipLine(l.source,l.direction,*this,tmin,tmax);
}

bool ConvexPolygon2D::intersects(const Line2D& l) const
{
	Real tmin=-Inf, tmax=Inf;
	return intersects(l,tmin,tmax);
}

int ConvexPolygon2D::planeIntersections(Plane2D& p,int& e1,int& e2,Real& u1,Real& u2) const
{
  int num=0;
  Real di,dj=p.distance(vertices[0]);
  for(size_t i=0;i<vertices.size();i++) {
    di=dj;
    dj=p.distance(vertices[next(i)]);
    if((di<Zero && dj>Zero) ||
       (di>Zero && dj<Zero) ||
       (di==Zero && dj!=Zero)) {
      if(num==0) {
	u1=di/(di-dj);
	e1=i;
      }
      else if(num==1) {
	u2=di/(di-dj);
	e2=i;
      }
      else {
	cout<<"More than 1 intersection???"<<endl;
	abort();
      }
      num++;
    }
  }
  return num;
}

//if v represents a circular list, erase the elements [a..b)
template <class T>
typename std::vector<T>::iterator CircularDelete(std::vector<T>& v,int a,int b)
{
  if(a>b) { //wraps around
    v.erase(v.begin()+a,v.end());
    return v.erase(v.begin(),v.begin()+b);
  }
  else {
    return v.erase(v.begin()+a,v.begin()+b);
  }
}

void ConvexPolygon2D::setIntersection(const ConvexPolygon2D& a,const ConvexPolygon2D& b)
{
  *this = a;
  //progressively do halfspace intersection for each plane of b
  Plane2D p;
  int e1,e2;
  Real u1,u2;
  for(size_t i=0;i<b.vertices.size();i++) {
    b.getPlane(i,p);
    //get the negative side of the halfspace
    int n=planeIntersections(p,e1,e2,u1,u2);
    if(n==0) {
      if(p.distance(vertices[0]) > Zero) { //all on outside of b
	vertices.clear();
	return;
      }
    }
    else if(n==1) {
      Assert(u1==Zero);
      //split at e1
      if(p.distance(vertices[next(e1)]) > Zero) { //all on outside of b
	//clear to a point
	Vector2 v=vertices[e1];
	vertices.clear();
	vertices.push_back(v);
      }
    }
    else if(n==2) {
      //split between e1 and e2
      //order them so e1 is on inside, e2 on outside
      if(p.distance(vertices[e1]) > Zero) {
	Assert(p.distance(vertices[e2]) <= Zero);
	swap(e1,e2);
	swap(u1,u2);
      }
      Vector2 v1,v2;
      interpolate(vertices[e1],vertices[next(e1)],u1,v1);
      interpolate(vertices[e2],vertices[next(e2)],u2,v2);
      //old order is ... e1 f1 ... e2 f2 ...
      //new order is ... e1 v1 v2 f2 ...
      vector<Vector2>::iterator f2=CircularDelete(vertices,e1+1,e2+1);
      vertices.insert(vertices.insert(f2,v2),v1);
    }
  }
}
