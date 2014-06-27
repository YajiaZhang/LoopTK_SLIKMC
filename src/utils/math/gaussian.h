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

#ifndef MATH_GAUSSIAN_H
#define MATH_GAUSSIAN_H

#include "MatrixTemplate.h"

namespace Math {

/** @ingroup Math
 * @brief Multivariate gaussian N(mu,K) of d dimensions
 *
 * \f$ P(x~N(mu,K)) = c e^{-1/2 (x-mu)^t K^{-1} (x-mu))} \f$
 * where c = 1/sqrt((2pi)^d |K|).
 *
 * Represented by the cholesky decomposition of K = LL^t.
 *
 * \f$ x^t K^{-1} x = 
 *     x^t (L L^t)^{-1} x = 
 *     x^t L^{-t} L^{-1} x = |L^{-1} x|^2 \f$
 *
 * so sqrt(|K|) = |L| = product of diag(L).
 * Equivalent to the change of variable y = L^-1(x-mu), such that
 * P(x~N(mu,K)) = P(y~N(0,I))
 */
template <class T>
class Gaussian
{
public:
  typedef Gaussian<T> MyT;
  typedef VectorTemplate<T> VectorT;
  typedef MatrixTemplate<T> MatrixT;

  Gaussian();
  Gaussian(int d);
  Gaussian(const MatrixT& sigma, const VectorT& mu);
  Gaussian(const MyT& g);
  
  void resize(int d);
  int numDims() const;
  
  ///Covariance matrix must be symmetric, strictly positive definite
  bool setCovariance(const MatrixT& sigma);
  void getCovariance(MatrixT& sigma) const;
  void setMean(const VectorT& mu);
  const VectorT& getMean() const { return mu; }
  
  ///Evaluates the probability of x
  T probability(const VectorT& x) const;
  ///Generates a point x according to this distribution
  void generate(VectorT& x) const;
  
  MatrixT L;		///< Cholesky decomposition of the covariance matrix
  VectorT mu;		///< mean
};

}

#endif
