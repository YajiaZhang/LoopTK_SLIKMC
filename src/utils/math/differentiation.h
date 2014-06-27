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

#ifndef MATH_DIFFERENTIATION_H
#define MATH_DIFFERENTIATION_H

#include "function.h"

/** @file math/differentiation.h"
 * @ingroup Math
 * @brief Numerical differentiation routines.
 */

namespace Math {
/** @addtogroup Math */
/*@{*/

/// centered differences
/// f'(x) ~= 1/2h*(f(x+h) - f(x-h)) + O(h^2)
Real dfCenteredDifference(RealFunction& f, Real x, Real h);
void dfCenteredDifference(VectorFunction& f, Real x, Real h, Vector& df);

/// centered differences for 2nd derivative
/// f''(x) ~= 1/h^2*(f(x+h) - 2f(x) + f(x-h)) + O(h^2)
Real ddfCenteredDifference(RealFunction& f, Real x, Real h);
void ddfCenteredDifference(VectorFunction& f, Real x, Real h, Vector& ddf);

/// x is provided as a temporary, it's restored to its original value later
void GradientForwardDifference(ScalarFieldFunction& f,Vector& x,Real h,Vector& g);
void JacobianForwardDifference(VectorFieldFunction& f,Vector& x,Real h,Matrix& J);
void HessianForwardDifference(ScalarFieldFunction& f,Vector& x,Real h,Matrix& H);
void HessianForwardDifference_Grad(ScalarFieldFunction& f,Vector& x,Real h,Matrix& H);

void GradientCenteredDifference(ScalarFieldFunction& f,Vector& x,Real h,Vector& g);
void JacobianCenteredDifference(VectorFieldFunction& f,Vector& x,Real h,Matrix& J);
void HessianCenteredDifference(ScalarFieldFunction& f,Vector& x,Real h,Matrix& H);
void HessianCenteredDifference_Grad(ScalarFieldFunction& f,Vector& x,Real h,Matrix& H);

/// specified hi in each direction xi
void GradientForwardDifference(ScalarFieldFunction& f,Vector& x,const Vector& h,Vector& g);
void JacobianForwardDifference(VectorFieldFunction& f,Vector& x,const Vector& h,Matrix& J);
void HessianForwardDifference(ScalarFieldFunction& f,Vector& x,const Vector& h,Matrix& H);
void HessianForwardDifference_Grad(ScalarFieldFunction& f,Vector& x,const Vector& h,Matrix& H);

void GradientCenteredDifference(ScalarFieldFunction& f,Vector& x,const Vector& h,Vector& g);
void JacobianCenteredDifference(VectorFieldFunction& f,Vector& x,const Vector& h,Matrix& J);
void HessianCenteredDifference(ScalarFieldFunction& f,Vector& x,const Vector& h,Matrix& H);
void HessianCenteredDifference_Grad(ScalarFieldFunction& f,Vector& x,const Vector& h,Matrix& H);

/*@}*/
} //namespace Math

#endif
