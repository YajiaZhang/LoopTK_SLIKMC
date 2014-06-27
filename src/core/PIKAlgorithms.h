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

/*
 * Contains inverse kinematic
 * algorithms to be called on
 * chains.
*/

#ifndef PIK_ALGORITHMS_H
#define PIK_ALGORITHMS_H

#include "PBasic.h"

typedef vector<ChainMove> IKSolution;
typedef vector<IKSolution> IKSolutions;

/**
 * Contains methods to perform cyclic coordinate
 * descent (CCD) for loop closure.
 */

class PCCDSolver {
  public:

    /**
     * Constructs a new PCCDSolver.
     */

    PCCDSolver(PChain *loop, PAtom *effectorPrior, PAtom *effectorEnd, Vector3 goalPrior, Vector3 goalEnd);

    /**
     * Performs one iteration of the CCD algorithm. Returns the sum of
     * squares of the distances between the end effectors and their goals.
     */

    Real DoDescent();

    void GetGoals(Vector3 &goalPrior, Vector3 &goalEnd);

  private:
    PChain *m_chain;
    Vector3 m_goalPrior, m_goalEnd;
    PAtom *m_effectorPrior, *m_effectorEnd;

    /* Degrees of freedom in the chain. */
    vector<PBond *> m_DOFs;

    /* Current iteration of CCD algorithm. */
    int m_ccdIter;

    void UpdateSides(const Vector3 &axisOfRotation, const Vector3 &basePrior, const Vector3 &atom, const Vector3 &goal, Real &k1, Real &k2);

};

/**
 * Simplifies the construction of a CCD solver
 * for a <code>PProtein</code>, allowing the user to pass
 * only a pointer to the loop and, if desired,
 * the direction as a <code>ProteinSide</code> parameter.
 * The standard <code>PCCDSolver</code> functions can then
 * be used to perform cyclic coordinate descent.
 */

class PProteinCCDSolver {
  public:

    /**
     * Constructs a new <code>PProteinCCDSolver</code> for
     * the specified <code>loop</code>.  The direction is
     * selected automatically based on distances
     * between anchor and end effector atoms:
     * whichever <code>ProteinSide</code> (start/end) has the
     * greater C-N anchor-effector distance is used.
     */

    PProteinCCDSolver(PProtein *loop);

    /**
     * Constructs a new PProteinCCDSolver for
     * the specified protein, using the specified
     * ProteinSide to select anchors/end effectors.
     */

    PProteinCCDSolver(PProtein *loop, ProteinSide side);
    ~PProteinCCDSolver() { delete m_solver; }

    /**
     * Performs the same function as <code>PCCDSolver::DoDescent()</code>.
     *
     * @see PCCDSolver::DoDescent
     */

    Real DoDescent() { return m_solver->DoDescent(); }

  private:
    PCCDSolver *m_solver;

    void ConstructorHelper(PProtein *loop, ProteinSide side, Vector3 goalPrior, Vector3 goalEnd);
};


/**
 * Finds exact inverse kinematic solutions for protein loops.
 */

class PExactIKSolver {
  public:

/**
 Given a closed loop and three residue indices, it finds Exact IK solutions using Coutsias method. Note for a certain pose Exact IK can give maximum 16 solutions. The loop in its current state corresponds to one of the solutions.
 */

    static IKSolutions FindSolutions(PProtein *loop,int Res_indices_to_use[3]);
/**
 Given a loop, three residue indices, and a goal pose defined by three atoms, it finds Exact IK solutions using Coutsias method. The loop conformation can be open or closed.
 */

    static IKSolutions FindSolutions(PProtein *loop, int Res_indices_to_use[3], Vector3 *endPriorG, Vector3 *endG, Vector3 *endNextG);
    //Yajia added this function.
    static int PExactIKSolver::FindSolutions(PProtein *loop, int DOF_indices_to_use[3], Vector3 *endPriorG, Vector3 *endG, Vector3 *endNextG, IKSolutions& solutions);
  private:
    static double AngleBetweenVectors(Vector3 n1,Vector3 n2);

};


#endif
