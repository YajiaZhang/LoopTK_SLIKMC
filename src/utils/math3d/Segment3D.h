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

#ifndef MATH3D_SEGMENT3D_H
#define MATH3D_SEGMENT3D_H

#include "Line3D.h"

namespace Math3D {

struct Segment3D
{
  void setTransformed(const Segment3D&, const Matrix4& xform);
  void getLine(Line3D& l) const;
  void getAABB(AABB3D& bb) const;
  Real closestPointParameter(const Point3D& in) const;
  Real closestPoint(const Point3D& in, Point3D& out) const;  //returns the parameter value of the point
  void closestPoint(const Segment3D&,Real& t,Real& u) const;  //same comment as above
  void eval(Real t, Point3D& out) const;
  bool intersects(const AABB3D&) const;
  bool intersects(const AABB3D&, Real& tmin, Real& tmax) const;
  bool Read(File& f);
  bool Write(File& f) const;
  
  Point3D A;
  Point3D B;
};

} //namespace Math3D

#endif
