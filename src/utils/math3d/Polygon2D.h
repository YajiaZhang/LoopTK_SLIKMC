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

#ifndef MATH3D_POLYGON2D_H
#define MATH3D_POLYGON2D_H

#include "Plane2D.h"
#include "Segment2D.h"
#include <vector>

namespace Math3D {

/** @ingroup Math3D
 * @brief An arbitrary connected polygon given by a vertex list
 *
 * Vertices are numbered 0...n-1. 
 * Edge i connects verts i and i+1 (looping around to 0).
 */
struct Polygon2D
{
  inline size_t next(size_t i) const { return (i+1>=vertices.size()? 0: i+1); }
  inline size_t prev(size_t i) const { return (i==0? vertices.size()-1: i-1); }
  bool ccw() const;
  bool nonCrossing() const;
  void setTransformed(const Polygon2D& in, const Matrix3& T);
  
  void getEdge(int i,Segment2D& ei) const;
  ///The plane normals point to the right of the edge
  void getPlane(int i,Plane2D& pi) const;
  bool planeSplits(const Plane2D& p) const;				///<returns true if the plane intersects this polyhedron
  bool planePos(const Plane2D& p) const;				///<returns true if this is entirely on the positive side of this plane
  bool planeNeg(const Plane2D& p) const;				///<returns true if this is entirely on the negative side of this plane
  bool raySplits(const Vector2& a,const Vector2& b) const;	///<returns true if the ray a->b intersects this polyhedron
  bool rayLeft(const Vector2& a,const Vector2& b) const;			///<returns true if this poly is entirely on the left side of a->b
  bool rayRight(const Vector2& a,const Vector2& b) const;			///<returns true if this poly is entirely on the right side of a->b
  int residue(const Vector2& x) const;     ///< 1 if inside, 0 if outside

  bool intersectsBoundary(const Polygon2D& other) const;		///<returns true if the boundaries intersect
  Real boundaryDistance(const Point2D& v) const; 
  bool intersects(const Line2D& l, Real& tmin, Real& tmax) const;
  bool intersects(const Line2D& l) const;
  bool intersects(const Segment2D& l) const;

  void getAABB(AABB2D&) const;
  
  bool Read(File& f);
  bool Write(File& f) const;

  ///a list of points around the boundary
  std::vector<Point2D> vertices;
};

/** @ingroup Math3D
 * @brief A convex connected polygon such that the vertices are in ccw order.
 */
struct ConvexPolygon2D : public Polygon2D
{
  bool isValid() const;

  bool intersects(const ConvexPolygon2D& other) const; 
  bool contains(const Point2D& v) const;
  bool withinDistance(const Point2D& v, Real dist) const;	
  bool withinEdgeDistance(const Point2D& v, Real dist) const;	
  Real distance(const Point2D& v) const;
  Real edgeDistance(const Point2D& v) const;
  bool intersects(const Line2D& l, Real& tmin, Real& tmax) const;
  bool intersects(const Line2D& l) const;
  bool intersects(const Segment2D& l) const; 

  ///returns either 0,1, or 2, and the edge indices and parameters
  int planeIntersections(Plane2D& p,int& e1,int& e2,Real& u1,Real& u2) const;
  void setIntersection(const ConvexPolygon2D& a,const ConvexPolygon2D& b);
};

} //namespace Math3D

#endif
