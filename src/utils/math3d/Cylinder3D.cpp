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

#include "Cylinder3D.h"
#include <math/misc.h>
#include <assert.h>
using namespace Math3D;
using namespace std;

void Cylinder3D::getBase(Circle3D& c) const
{
  c.center=center;
  c.axis=axis;
  c.radius=radius;
}

void Cylinder3D::getCap(Circle3D& c) const
{
  c.center=center; c.center.madd(axis,height);
  c.axis=axis;
  c.radius=radius;
}

void Cylinder3D::getAABB(AABB3D& aabb) const
{
  Circle3D c;
  getBase(c);
  c.getAABB(aabb);
  if(axis.x > 0) aabb.bmax.x += axis.x*height;
  else aabb.bmin.x -= axis.x*height;
  if(axis.y > 0) aabb.bmax.y += axis.y*height;
  else aabb.bmin.y -= axis.y*height;
  if(axis.z > 0) aabb.bmax.z += axis.z*height;
  else aabb.bmin.z -= axis.z*height;
}

bool Cylinder3D::intersects(const Line3D& line,Real* tmin,Real* tmax) const
{
  Real axistmin,axistmax;

  //quick reject - infinite cylinder
  Vector3 src=line.source-center;
  const Vector3& dir=line.direction;
  assert(FuzzyEquals(axis.normSquared(),One));
  //quadratic equation
  Real a,b,c;
  Real dv,sv;
  dv = dir.dot(axis);
  sv = src.dot(axis);
  a=dir.normSquared()-Sqr(dv);
  b=Two*(dir.dot(src)-sv*dv);
  c=src.normSquared()-Sqr(sv)-Sqr(radius);
  int nroots=quadratic(a,b,c,axistmin,axistmax);
  //TODO: if the line is contained within the cylinder, ignore this
  if(nroots == 0) return false;
  else if(nroots == 1) axistmax=axistmin;
  else if(nroots == 2) {
    if(axistmin > axistmax) std::swap(axistmin,axistmax);
  }
  else if(nroots == -1) return false;
  else return false;  //what case is this?

  //projection of intersection pts on the cyl axis
  Real axisumin,axisumax;
  Vector3 temp;
  line.eval(axistmin,temp);  axisumin = axis.dot(temp-center);
  line.eval(axistmax,temp);  axisumax = axis.dot(temp-center);

  //now check the caps
  Real tc;
  Circle3D cir;
  if(axisumin < 0) {   //hits a cap first
    if(dv > 0) //line points along axis
      getBase(cir);
    else  //line points against axis
      getCap(cir);
    if(!cir.intersects(line,&tc)) return false;
    axistmin = tc;
  }
  if(axisumin > height) {  //hits a cap last
    if(dv > 0) //line points along axis
      getCap(cir);
    else  //line points against axis
      getBase(cir);
    if(!cir.intersects(line,&tc)) return false;
    axistmax = tc;
  }
  //if(axistmin > axistmax) return false;
  if(tmin) *tmin=axistmin;
  if(tmax) *tmax=axistmax;
  return true;
}
