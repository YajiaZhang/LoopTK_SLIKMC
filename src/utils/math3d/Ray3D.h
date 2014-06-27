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

#ifndef MATH3D_RAY3D_H
#define MATH3D_RAY3D_H

#include "Line3D.h"

namespace Math3D {

struct Ray3D : public Line3D
{
  Real closestPoint(const Point3D& in, Point3D& out) const;
  Real distance(const Point3D& pt) const;
  bool intersects(const Line3D&, Real* t=NULL, Real* u=NULL, Real epsilon=0) const;	//t is the parameter of this ray, u is the line 
  bool intersects(const Ray3D&, Real* t=NULL, Real* u=NULL, Real epsilon=0) const;	//t is the parameter of this ray, u is the other ray
  void closestPoint(const Line3D&,Real& t,Real& u) const;  //same comment as above
  void closestPoint(const Ray3D&,Real& t,Real& u) const;  //same comment as above
};


} //namespace Math3D

#endif
