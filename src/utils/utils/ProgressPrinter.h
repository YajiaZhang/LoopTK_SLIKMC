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

#ifndef UTILS_PROGRESS_PRINTER_H
#define UTILS_PROGRESS_PRINTER_H

#include <iostream>

/** @ingroup Utils
 * @brief Prints the progress of an iterative, long computation.
 *
 * The user calls Update() each iteration from 0...max.
 * Print() is called "increments" times at evenly-spaced
 * intervals.  By default, prints the percentage complete to
 * "out"
 */
class ProgressPrinter
{
public:
  ProgressPrinter(std::ostream& out,int max,int increments=100);
  ProgressPrinter(int max,int increments=100);
  void Update();
  void Update(int iter);
  virtual void Print(float fraction);
  virtual void Done();

  std::ostream& out;
  int max;
  int increments;
  int iter;
};

#endif
