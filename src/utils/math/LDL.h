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

#ifndef MATH_LDL_H
#define MATH_LDL_H

#include "matrix.h"

namespace Math {

/** @ingroup Math
 * @brief Performs the LDL^t decompositoin of a symmetric matrix A.
 *
 * L is stored in the lower diagonal of LDL, D is stored in the diagonal.
 * If any of the elements of D's diagonal have abs value less than 
 * zeroTolerance (i.e. A is singular), then they are ignored.
 */
template <class T>
struct LDLDecomposition
{
  typedef MatrixTemplate<T> MatrixT;
  typedef VectorTemplate<T> VectorT;

  LDLDecomposition();
  LDLDecomposition(const MatrixT& A);

  void set(const MatrixT& A);
  void backSub(const VectorT& b, VectorT& x) const;
  void LBackSub(const VectorT& b, VectorT& x) const;
  void LTBackSub(const VectorT& b, VectorT& x) const;
  void DBackSub(const VectorT& b, VectorT& x) const;
  void getInverse(MatrixT& Ainv) const;
  void getL(MatrixT& L) const;
  void getD(VectorT& d) const;
  void getA(MatrixT& A) const;
  void mulL(const Vector& x,Vector& y) const;
  void mulLT(const Vector& x,Vector& y) const;
  void mulD(const Vector& x,Vector& y) const;

  /// Update the LDL decomposition for A + xx^t
  void update(const VectorT& x);
  /// "Downdate" the LDL decomposition for A - xx^t (on failure, LDL is undefined)
  bool downdate(const VectorT& x);

  MatrixT LDL;
  T zeroTolerance;
};

} // namespace Math

#endif
