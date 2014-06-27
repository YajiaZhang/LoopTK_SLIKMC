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

#ifndef POPTIMIZE_H
#define POPTIMIZE_H

#include "PIKAlgorithms.h"
#include "PConstants.h"
#include "PExtension.h"
#include "PBasic.h"

#define  pi 3.14159265
#define deg2rad pi/180.0e0
#define rad2deg 180.0e0/pi
/**
 * @package Optimization
 *
 * Contains methods to Optimize specified functions over loops.
 */
/** This structure contains the DOF solution and the non-DOF solution after optimization*/
struct OptimalSol {
  vector<IKSolutions> LoopSol;
  vector<double> NonLoopSol;
};

class POptimize{
 public:
   /**
   * Constructs a new <code>POptimize</code> with
   * specified <code>loops</code> and <code>Dofs</code> of the loops. The optimization
   * will be carried over the <code>Dofs</code> of the <code>loops</code>.
   */

	POptimize(vector<PProtein*> loops, vector< vector<CDof> > Dofs);
   
   /**
   * Returns <code>OptimalSol</code> corresponding to <code>loops</code> after
   * optimizing specified functor <code>FunctToOptimize</code>. Numerical derivatives
   * are calculated and projected on to the null space.
   */
	OptimalSol OptimizeNullSpace(FunctFunctor *FunctToOptimize);

   /**
   * Returns <code>OptimalSol</code> corresponding to <code>loops</code> after
   * optimizing specified functor <code>FunctToOptimize</code> with specified
   * derivative functor <code>DerviOfFunct</code>.
   */
	OptimalSol Optimize(FunctFunctor *FunctToOptimize, DerivFunctor *DerivOfFunct);
   /**
   * Returns <code>OptimalSol</code> corresponding to <code>loops</code> after
   * optimizing specified functor <code>FunctToOptimize</code>. Numerical derivatives
   * are calculated and projected on to the null space.In this method, one can specify the initial values of the non-DOF variables in <code>InitC</code>.*/

	OptimalSol OptimizeNullSpace(FunctFunctor *FunctToOptimize, vector<double> InitC);

   /**
   * Returns <code>OptimalSol</code> corresponding to <code>loops</code> after
   * optimizing specified functor <code>FunctToOptimize</code> with specified
   * derivative functor <code>DerviOfFunct</code>.In this method, one can specify the initial values of the non-DOF variables in <code>InitC</code>.*/
	OptimalSol Optimize(FunctFunctor *FunctToOptimize, DerivFunctor *DerivOfFunct, vector<double> InitC);

  private:
	IKSolution IKCloseChain(PProtein *lp, Vector3 endPriorGoal, Vector3 endGoal, Vector3 endNextGoal);
	vector<PProtein*> lps;
	vector<int> nRsds, nbb, ndofs;
	vector<Vector3> endPG,endG,endNG;
	vector<ChainMove> AllMove;
	vector< vector<CDof> > DofsToUse;
};

#endif
