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

#ifndef MATH_METRIC_H
#define MATH_METRIC_H

#include "VectorTemplate.h"
#include "MatrixTemplate.h"

/** @file metric.h 
 * @ingroup Math
 * @brief Standard vector/matrix metrics
 */

namespace Math {

/** @addtogroup Math */
/*@{*/

template<class T>
T Norm_L1(const VectorTemplate<T>& x);
template<class T>
T Norm_L2(const VectorTemplate<T>& x);
///Same as above, but robust to over/underflow
template<class T>
T Norm_L2_Safe(const VectorTemplate<T>& x);
template<class T>
T Norm_LInf(const VectorTemplate<T>& x);
template<class T>
T Norm_Mahalanobis(const VectorTemplate<T>& x,T k);

template<class T> 
T Distance_L1(const VectorTemplate<T>& x,const VectorTemplate<T>& y);
template<class T>
T Distance_L2(const VectorTemplate<T>& x,const VectorTemplate<T>& y);
///Same as above, but robust to over/underflow
template<class T>
T Distance_L2_Safe(const VectorTemplate<T>& x,const VectorTemplate<T>& y);
template<class T>
T Distance_LInf(const VectorTemplate<T>& x,const VectorTemplate<T>& y);
template<class T>
T Distance_Mahalanobis(const VectorTemplate<T>& x,const VectorTemplate<T>& y,T k);

template<class T>
T Norm_L1(const MatrixTemplate<T>& A);
template<class T>
T Norm_LInf(const MatrixTemplate<T>& A);
template<class T>
T Norm_Frobenius(const MatrixTemplate<T>& A);
///Same as above, but robust to over/underflow
template<class T>
T Norm_Frobenius_Safe(const MatrixTemplate<T>& A);

template<class T>
T Distance_L1(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B);
template<class T>
T Distance_LInf(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B);
template<class T>
T Distance_Frobenius(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B);
///Same as above, but robust to over/underflow
template<class T>
T Distance_Frobenius_Safe(const MatrixTemplate<T>& A,const MatrixTemplate<T>& B);

/*@}*/
} //namespace Math

#endif
