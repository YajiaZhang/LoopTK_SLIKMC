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

#ifndef MATH_LINEAR_PATH_H
#define MATH_LINEAR_PATH_H

#include "function.h"
#include "Sequence.h"
#include <vector>

namespace Math {

class LinearPath : public VectorFunction
{
public:
  struct ControlPoint {
    Real t;
    Vector x;
  };

  virtual void PreEval(Real t);
  virtual void Eval(Real t,Vector& x);
  virtual void Deriv(Real t,Vector& dx);

  //overridable functions
  //Interpolate x from a->b
  //Set to be the difference x=a-b
  virtual void Interpolate(const Vector& a,const Vector& b,Real u,Vector& x) const;
  virtual void Difference(const Vector& a,const Vector& b,Vector& x) const;
  virtual Real Distance(const Vector& a,const Vector& b) const;

  inline const Vector& operator[](int i) const { return points[i].x; }
  inline Vector& operator[](int i) { return points[i].x; }
  inline size_t size() const { return points.size(); }
  std::vector<ControlPoint>::iterator GetSegment(Real t);
  inline Real BeginTime() const { return points.front().t; }
  inline Real EndTime() const { return points.back().t; }
  Real Length() const;
  void SetSequence(const VectorSequence& s);
  void ArcLengthParameterize();
  void ScaleTime(Real s);
  void OffsetTime(Real off);
  void Concat(const LinearPath& p);  //adds p onto the end of this
  void Append(const Vector& x,Real dt=Zero);

  bool Read(File& f);
  bool Write(File& f) const;

  std::vector<ControlPoint> points;
};

} //namespace Math

#endif
