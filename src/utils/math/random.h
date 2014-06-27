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

#ifndef MATH_RANDOM_H
#define MATH_RANDOM_H

#include <utils/random.h>
#include "math.h"

/** @file math/random.h
 * @ingroup Math
 * @brief Defines a standard method for random floating-point number
 * generation.
 *
 * The random number generator is defined as Math::rng.
 */

/** @addtogroup Math
 * @brief Generates random gaussian distributed floating point values.
 *
 * Uses a given uniform random number generator, generates values
 * in a gaussian with the given mean and standard deviation.
 */
template <class RNG>
float frand_gaussian(RNG& rng,float mean, float stddev);
template <class RNG>
double drand_gaussian(RNG& rng,double mean, double stddev);

namespace Math {
  /** @addtogroup Math */
  /*@{*/

extern StandardRNG rng;
inline void Srand(unsigned long seed) { rng.seed(seed); }
inline long int RandInt() { return rng.randInt(); }
inline long int RandInt(long int n) { return rng.randInt(n); }

#ifdef MATH_DOUBLE
inline Real Rand() { return rng.randDouble(); }
inline Real Rand(Real a,Real b) { return rng.randDouble(a,b); }
inline Real RandGaussian(Real mean, Real stddev) { return drand_gaussian(rng,mean,stddev); }
#else
inline Real Rand() { return rng.randFloat(); }
inline Real Rand(Real a,Real b) { return rng.randFloat(a,b); }
inline Real RandGaussian(Real mean, Real stddev) { return frand_gaussian(rng,mean,stddev); }
#endif //MATH_DOUBLE

/** @fn Srand()
 * \brief Seeds the Math random number generator.
 */

/** @fn Rand() 
 * @brief Generates a random Real uniformly in [0,1].
 */

/** @fn Rand(Real,Real)
 * @brief Generates a random Real uniformly in [a,b].
 */

/** @fn RandInt() 
 * @brief Generates a random int.
 */

/** @fn RandInt(int) 
 * @brief Generates a random int in [0,n).
 */

/** @fn RandGaussian(Real,Real) 
 * \brief Generates a random Real in a Gaussian distribution.
 */

/*@}*/
} //namespace Math



//definition of rand_gaussian functions

template <class RNG>
float frand_gaussian(RNG& rng, float mean, float stddev)
{
  static float t = 0.0f;
  float x,v1,v2,r;
  if (t == 0) {
    do {
      v1 = 2.0f * rng.randFloat() - 1.0f;
      v2 = 2.0f * rng.randFloat() - 1.0f;
      r = v1 * v1 + v2 * v2;
    } while (r>=1.0f);
    r = sqrtf((-2.0f*logf(r))/r);
    t = v2*r;
    return(mean+v1*r*stddev);
  }
  else {
    x = t;
    t = 0.0f;
    return(mean+x*stddev);
  }
}

template <class RNG>
double drand_gaussian(RNG& rng, double mean, double stddev)
{
  static double t = 0.0;
  double x,v1,v2,r;
  if (t == 0) {
    do {
      v1 = 2.0 * rng.randDouble() - 1.0;
      v2 = 2.0 * rng.randDouble() - 1.0;
      r = v1 * v1 + v2 * v2;
    } while (r>=1.0);
    r = sqrt((-2.0*log(r))/r);
    t = v2*r;
    return(mean+v1*r*stddev);
  }
  else {
    x = t;
    t = 0.0;
    return(mean+x*stddev);
  }
}

#endif
