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

#ifndef MATH_BACK_SUBSTITUTE_H
#define MATH_BACK_SUBSTITUTE_H

#include "MatrixTemplate.h"

namespace Math {

// If A is upper triangular nxn, solves Ax=b
template <class T>
bool UBackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x);

// If A is lower triangular nxn, solves Ax=b
template <class T>
bool LBackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x);

// If A is lower triangular nxn, solves A^t*x=b
template <class T>
bool LtBackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x);

// The following are identical to the above, but assume the diagonal of a is all ones
template <class T>
void U1BackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x);
template <class T>
void L1BackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x);
template <class T>
void Lt1BackSubstitute(const MatrixTemplate<T>& a, const VectorTemplate<T>& b, VectorTemplate<T>& x);


} // namespace Math

#endif
