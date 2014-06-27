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

#ifndef MATH_SAMPLE_H
#define MATH_SAMPLE_H

#include <vector>
#include "IntervalSet.h"

/** @file math/sample.h
 * @brief Functions for random sampling of various sets.
 */

namespace Math {
/** @addtogroup Math */
/*@{*/

/** @brief Samples an integer with weighted probability.
 *
 * The probability of sampling integer i in [0,n) is wi / W
 * where W = sum of all wi.
 */
int WeightedSample(const std::vector<Real>& weights);
/// Same as above, but the sum of all weights, W, is provided.
int WeightedSample(const std::vector<Real>& weights,Real totWeight);

/// Uniformly samples the given intervals
Real Sample(const Interval& s);
/// Uniformly samples the interval set
Real Sample(const ClosedIntervalSet& s);

/// Uniform distribution on boundary of a circle with radius r
void SampleCircle(Real r, Real& x, Real& y);
/// Uniform distribution inside a circle with radius r
void SampleDisk(Real r, Real& x, Real& y);

/// Uniform distribution in the triangle whose vertices are (0,0),(1,0),(0,1)
void SampleTriangle(Real& x, Real& y);

/// Uniform distribution on boundary of a sphere with radius r
void SampleBall(Real r, Real& x, Real& y, Real& z);
/// Uniform distribution inside a sphere with radius r
void SampleSphere(Real r, Real& x, Real& y, Real& z);

} //namespace Math
/*@}*/

#endif
