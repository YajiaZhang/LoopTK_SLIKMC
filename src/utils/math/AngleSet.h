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

#ifndef MATH_ANGLE_SET_H
#define MATH_ANGLE_SET_H

#include "angle.h"
#include <vector>

namespace Math {

/** @ingroup Math
 * @brief A set of AngleIntervals.
 */
struct AngleSet : public std::vector<AngleInterval>
{
  typedef std::vector<AngleInterval> BaseT;

  inline void SetCircle() { resize(1); (*this)[0].setCircle(); }
  inline void SetEmpty() { clear(); }
  void Union(const BaseT&);
  void Intersect(const BaseT&);
  void Subtract(const BaseT&);
  void Union(const AngleInterval&);
  void Intersect(const AngleInterval&);
  void Subtract(const AngleInterval&);
};

} //namespace Math

#endif
