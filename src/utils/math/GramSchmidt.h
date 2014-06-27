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

#ifndef MATH_GRAM_SCHMIDT_H
#define MATH_GRAM_SCHMIDT_H

#include "VectorTemplate.h"

namespace Math {

  /** @addtogroup Math */
  /*@{*/

///Performs the Gram-Schmidt process to get an orthonormal basis
///for the span of the input vectors X.
///Returns the number of nonzero vectors in the output basis
template <class T>
int OrthonormalBasis(const VectorTemplate<T>* x, VectorTemplate<T>* basis, int n);

///Same as above, but does not normalize
template <class T>
int OrthogonalBasis(const VectorTemplate<T>* x, VectorTemplate<T>* basis, int n);

///orthogonalizes a vector w.r.t the orthogonal basis of n vectors
template <class T>
void Orthogonalize(VectorTemplate<T>& x,const VectorTemplate<T>* basis, int n);

/*@}*/
} //namespace Math

#endif
