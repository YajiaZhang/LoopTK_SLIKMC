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

#include "Ray3D.h"
using namespace Math3D;

Real Ray3D::closestPoint(const Point3D& in, Point3D& out) const
{
  return Line3D::closestPoint(in,out,0,Inf);
}

Real Ray3D::distance(const Point3D& pt) const
{
  Point3D closest;
  closestPoint(pt,closest);
  return (pt-closest).norm();
}

bool Ray3D::intersects(const Line3D& l, Real* t, Real* u, Real epsilon) const
{
  if(Line3D::intersects(l,t,u,epsilon)) {
    if(*t < -epsilon) return false;
    return true;
  }
  return false;
}

bool Ray3D::intersects(const Ray3D& r, Real* t, Real* u, Real epsilon) const
{
  if(Line3D::intersects(r,t,u,epsilon)) {
    if(*t < -epsilon || *u < -epsilon) return false;
    return true;
  }
  return false;
}

void Ray3D::closestPoint(const Line3D& l,Real& t,Real& u) const
{
  Line3D::closestPoint(l,t,u);
  if(t < 0) {
    t=0;
    u=l.closestPointParameter(source);
  }
}

void Ray3D::closestPoint(const Ray3D& r,Real& t,Real& u) const
{
  Line3D::closestPoint(r,t,u);
  if(t < 0) {
    t=0;
    u=r.closestPointParameter(source);
  }
  if(u < 0) {
    u=0;
    t=r.closestPointParameter(r.source);
  }
  if(t < 0) t=0;
}
