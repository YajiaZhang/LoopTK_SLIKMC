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

#ifndef MATH_INTERVAL_SET_H
#define MATH_INTERVAL_SET_H

#include "Interval.h"
#include <vector>

namespace Math {

struct OpenIntervalSet : public std::vector<OpenInterval>
{
  typedef std::vector<OpenInterval> BaseT;
  typedef std::vector<ClosedInterval> ClosedBaseT;

  inline void SetFull() { resize(1); (*this)[0].setFull(); }
  inline void SetEmpty() { clear(); }
  void Union(const BaseT&);
  void Intersect(const BaseT&);
  void Subtract(const ClosedBaseT&);
  void Union(const OpenInterval&);
  void Intersect(const OpenInterval&);
  void Subtract(const ClosedInterval&);
};

struct ClosedIntervalSet : public std::vector<ClosedInterval>
{
  typedef std::vector<ClosedInterval> BaseT;
  typedef std::vector<OpenInterval> OpenBaseT;

  inline void SetFull() { resize(1); (*this)[0].setFull(); }
  inline void SetEmpty() { clear(); }
  void Union(const BaseT&);
  void Intersect(const BaseT&);
  void Subtract(const OpenBaseT&);
  void Union(const ClosedInterval&);
  void Intersect(const ClosedInterval&);
  void Subtract(const OpenInterval&);
};

} //namespace Math

#endif
