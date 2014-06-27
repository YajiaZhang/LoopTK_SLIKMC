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

#ifndef MATH_CHOLESKY_DECOMPOSITION_H
#define MATH_CHOLESKY_DECOMPOSITION_H

#include "matrix.h"

namespace Math {

/** @ingroup Math
 * @brief Performs the Cholesky decomposition.
 *
 * Decomposes the positive semi-definite matrix A to LL^t.
 */
template <class T>
class CholeskyDecomposition
{
public:
  typedef MatrixTemplate<T> MatrixT;
  typedef VectorTemplate<T> VectorT;

  CholeskyDecomposition();
  CholeskyDecomposition(const MatrixT& A);
  CholeskyDecomposition(const MatrixT& A,MatrixT& L);

  void setDestination(MatrixT& L);
  bool set(const MatrixT& A);
  bool setPerturbed(const MatrixT& A,T& lambda);
  void backSub(const VectorT& b, VectorT& x) const;
  void LBackSub(const VectorT& b, VectorT& x) const;
  void LTBackSub(const VectorT& b, VectorT& x) const;
  void getInverse(MatrixT& Ainv) const;

  /// Update the cholesky decomposition for A + xx^t
  void update(const VectorT& x);
  /// "Downdate" the cholesky decomposition for A - xx^t (on failure, L is undefined)
  bool downdate(const VectorT& x);

  MatrixT L;
  T zeroEpsilon;
};

}
#endif
