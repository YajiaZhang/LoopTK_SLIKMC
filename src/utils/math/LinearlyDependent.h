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

#ifndef MATH_LINEARLY_DEPENDENT_H
#define MATH_LINEARLY_DEPENDENT_H

#include "VectorTemplate.h"
#include "MatrixTemplate.h"

namespace Math {

/** @ingroup Math
 * @brief Robust determination of linear dependence.
 *
 * Calculates a constant c s.t. they're linearly dependent within rel error 
 * eps and |c| <= 1.
 *
 *  2 cases, (1) a*c = b with maxabs(a*c-b)/|a| <= eps.
 *  (2) a = c*b with maxabs(a-c*b)/|b| <= eps.
 *  in case 1, cb = false, in case 2, cb=true.
 *  c is chosen by pseudoinverse
 * @return True if the vectors are dependent.
 */
template <class T>
bool LinearlyDependent_Robust(const VectorTemplate<T>& a, const VectorTemplate<T>& b, T& c, bool& cb, T eps = Epsilon);

/** @ingroup Math
 * @brief A matrix version of the function above.
 *
 * @return A vector c s.t. A*c = 0.  All coefficients |ci| <= 1
 *  At least one ci = 1.
 */
template <class T>
bool LinearlyDependent_Robust(const MatrixTemplate<T>& A, VectorTemplate<T>& c, T eps = Epsilon);


} //namespace Math

#endif 
