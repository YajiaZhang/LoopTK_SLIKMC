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

#ifndef MATH_MATRIX_EQUATION_PRINTER_H
#define MATH_MATRIX_EQUATION_PRINTER_H

#include "matrix.h"
#include <list>

namespace Math
{

struct EquationTerm
{
  EquationTerm();
  int Height() const;
  int Width() const;
  void PrintLine(int line,std::ostream& out) const;

  const Matrix* matrix;
  const Vector* vector;
  const char* text;
  Real scalar;
  bool transpose;
  bool ASCIIshade;
};

class MatrixEquationPrinter
{
public:
  void PushMatrix(const Matrix&);
  void PushMatrixTranspose(const Matrix&);
  void PushVector(const Vector&);
  void PushText(const char*);
  void PushScalar(Real);
  inline void PushAdd() { PushText("+"); }
  inline void PushSub() { PushText("-"); }
  inline void PushTimes() { PushText("*"); }
  inline void PushDot()  { PushText("."); }
  inline void PushEquals()  { PushText("="); }
  void Print(std::ostream& out) const;
  void Clear();

  std::list<EquationTerm> terms;
};


} //namespace Math

#endif
