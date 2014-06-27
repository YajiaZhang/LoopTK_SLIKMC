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

#ifndef MATH_ROW_ECHELON_H
#define MATH_ROW_ECHELON_H

#include "MatrixTemplate.h"
#include <vector>

namespace Math {

/** @ingroup Math
 * @brief Compute reduced row-eschelon form for a matrix A.
 *
 * Back-substitution of Ax=b is done by providing b at initialization.
 */
template <class T>
class RowEchelon
{
public:
  typedef MatrixTemplate<T> MatrixT;
  typedef VectorTemplate<T> VectorT;

  RowEchelon();
  RowEchelon(const MatrixT& A);
  RowEchelon(const MatrixT& A,const VectorT& b);
  RowEchelon(const MatrixT& A,const MatrixT& B);
  
  void set(const MatrixT& A);
  void set(const MatrixT& A,const VectorT& b);
  void set(const MatrixT& A,const MatrixT& B);
  void backSub(VectorT& x) const;
  int getRank() const;
  int getNull() const;
  /// Returns a orthonormal basis for the nullspace in the columns of N
  void getNullspace(MatrixT& N) const;
  /// Calculates a pseudoinverse x0 as well as a nullspace matrix N
  void getAllSolutions(VectorT& x0,MatrixT& N) const;

  //helpers
  void calcFirstEntries();
  
  MatrixT R;
  MatrixT EB;
  std::vector<int> firstEntry;  ///<indexes the first entry in each row of r
};

template <class T>
bool IsRowEchelon(const MatrixTemplate<T>& A);
template <class T>
bool IsReducedRowEchelon(const MatrixTemplate<T>& A);

} //namespace Math


#endif
