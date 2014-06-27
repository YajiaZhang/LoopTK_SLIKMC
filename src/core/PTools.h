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

#ifndef PTOOLS_H
#define PTOOLS_H


#include "PBasic.h"
#include "PExtension.h"
#include "PNumRoutines.h"
#include <gsl/gsl_sf_bessel.h>


typedef Real(*ConformationDistFn)(const PLightChain*, const PLightChain*);

struct NullSpaceRet{
	double *Sval;
	double **Svec;
	int *ns;
	int n_ns;
};

class PCluster;
//@package Inverse Kinematics
/**
 * PTools defines a variety of functions related
 * to loop closure, including functions to calculate goal positions,  computation of Jacobian, computation of Null Space, and computation of RMSD distance.
 */
class PTools {
  public:

    /**
     * Computes goal positions for chain loop closure. This version of
     * the overloaded function <code>GetGoals</code> does not make the
     * assumption that the chain in question represents a protein, and
     * so requires information about the location of both anchors and
     * end effectors.
     * 
     *  @param	anchorPrior	The anchor closest to the disconencted chain loop.
     *  @param	anchorEnd	The anchor farthest from the disconnected chain loop.
     *  @param	effectorPrior	The end effector whose goal is farthest from <code>anchorPrior</code>.
     *  @param	effectorEnd	The end effector whose goal is closest to <code>anchorPrior</code>.
     *  @param	len		Distance between <code>goalEnd</code> and <code>anchorPrior</code>.
     *  @param	angle1		Angle formed by <code>goalPrior</code>-<code>goalEnd</code>-<code>anchorPrior</code>.
     *  @param	angle2		Angle formed by <code>goalEnd</code>-<code>anchorPrior</code>-<code>anchorEnd</code>.
     *  @param	goalPrior	Reference to the <code>Vector3</code> in which to store the goal for <code>effectorPrior</code>.
     *  @param	goalEnd		Reference to the <code>Vector3</code> in which to store the goal for <code>effectorEnd</code>.
     */

    static void GetGoals(PAtom *anchorPrior, PAtom *anchorEnd, PAtom *effectorPrior, PAtom *effectorEnd,
                          Real len, Real angle1, Real angle2, Vector3 &goalPrior, Vector3 &goalEnd);

    /**
     * Computes goal positions for the specified <code>protein</code> on
     * the specified <code>ProteinSide</code>.
     */

    static void GetGoals(PProtein *protein, ProteinSide side, Vector3 &goalPrior, Vector3 &goalEnd);

    /**
     * Computes goal positions for the specified <code>protein</code>.
     * The <code>ProteinSide</code> direction is selected automatically 
     * based on distances between anchor and end effector atoms: whichever
     * side (<code>start</code>/<code>end</code>) has the greater C-N
     * anchor-effector distance is used.
     */

    static ProteinSide GetGoals(PProtein *protein, Vector3 &goalPrior, Vector3 &goalEnd);

    /**
     * Caller is responsible for deallocating return.
     */
    static PProtein *CreateSlimProtein(PProtein *protein);
    static PProtein *CreateSlimProtein(PProtein *protein, int startRes, int endRes);

    /**
     * Returns the distance between the <i>closest</i> anchor and end
     * effector atoms on the specified <code>side</code> of <code>loop</code>.
     */

    static Real EffectorAnchorDist(PProtein *loop, ProteinSide side);

    /**
     * Returns smallest shared chain or NULL if different top-level chains.
    */
    static PChain *LowestCommonChain(PChain *c1, PChain *c2);

    /**
     * Performs conformation space clustering on the specified <code>space</code>.
     * The user-provided function <code>compFn</code> is used to determine the
     * distance between each pair of conformations.
     */
    static list<PCluster> GetClusters(const PConformationSpace &space,
                                      ConformationDistFn distFn,
                                      Real threshold);

    /**
     * Returns the index (of atoms in the backbone) for atom
     * <code>atomID</code> in residue number <code>resIndex</code>.
     */
    static int getBackboneAtomIndex(int resIndex, const string &atomID);

    /**
     * Attempts to find a collision-free rotamer for the sidechain of residue
     * <code>res</code>.  If successful, <code>res</code> is updated to reflect
     * the new conformation and <code>ApplyRotamer</code> returns true.  If no
     * collision-free conformation can be found, <code>res</code> is unchanged
     * and <code>ApplyRotamer</code> returns false.
     */

    //static bool ApplyRotamer(PResidue *res);

    /** Computes Jacoboian of specified loop, using specified DOFs (ind), tool frame position(p) and
     * direction (forward/backward).
     */

    static void ComputeJacobian(PProtein *loop, vector<int>& ind, double **Jac, Vector3& p,bool forward);

     /**
     * Computes Jacoboian of specified loop, using specified DOFs and direction
     */

    static void ComputeJacobian(PProtein *loop, vector<int>& ind, double **Jac, bool forward);
    
     /**
     * Computes Jacoboian of specified loop, using all backbone DOFs, C-anchor's position, and direction
     */
 
    static void ComputeJacobian(PProtein *loop, double **Jac);

//    static void ComputeJacobian(PProtein* loop, int atom_index, double **Jac, bool forward = true);

    static void ComputeJacobian(PProtein* loop, double** Jac, Vector3& atom_pos);

    /**
     * Computes Singularity Robust Inverse of m x n matrix A, The result is contained in Ainv.
     */

    static void ComputePseudoInverse(double **A, int m, int n, double **Ainv);
    /**
     * Computes Null Space of the Jacobian <code>Jac</code>. The dimensions of Jac are 6 x dim or 3 x dim depending on <code>SixDimensional</code>
     */
    static void ComputeNullSpace(double **Jac, int dim, bool SixDimensional, NullSpaceRet* Ret);
    /**
     * Projects Vector <code>ToProject</code> on the null space of Jacobian of the <code>loop</code>, given the DOFs of loop(<code>ind</code>) to consider, direction of computation of Jacobian (<code>forward</code>). The Projected vector is returned in <code>AfterProject</code>.
     */
    static void ProjectOnNullSpace(PProtein *loop, vector<int> ind, bool forward, double ToProject[], double AfterProject[]);

    static void ProjectOnNullSpace(PProtein *loop, vector<int> ind, bool forward, double ToProject[], double AfterProject[], bool sixDimensional);
    /**
     * Returns all backbone DOFs of all <code>loops</code>
     */
    static vector<vector<CDof> > GetBBDofs(vector<PProtein *> loops);

    /**
     * Returns RMSD between C_alpha taoms of two proteins <code>protein0</code> and
     * <code>protein1</code> starting from residue <code>loopstart</code> and ending
     * at residue <code>loopend</code>.
     */
    static double RMSDCalpha (PProtein *protein0, PProtein *protein1, int loopstart, int loopend);


    /**
     * Returns RMSD between backbones of two proteins <code>protein0</code>
     * and <code>protein1</code> starting from residue <code>loopstart</code>
     * and ending at residue <code>loopend</code>.
     */
    static double RMSDBackbone(PProtein *protein0, PProtein *protein1, int loopstart, int loopend);

    /**
     * Returns RMSD between all the atoms of two proteins <code>protein0</code>
     * and <code>protein1</code> starting from residue <code>loopstart</code>
     * and ending at residue <code>loopend</code>.
     */
    static double RMSDAllAtom(PProtein *protein0, PProtein *protein1, int loopstart, int loopend);

    /**
     * Modifies a protein loop/ fragment specified by <code>lpD</code> to overlap it with a specified 
     *loop/fragment specified by <code>lpS</code>. The loops start respectively at <code>startS</code> and
     *<code>startD</code> and have <code>numRes</code> number of residues. The N-C-alpha atoms of the starting residues of the loops should have same coordinates.
     */
    static bool CopyBackbone(PProtein *lpS, PProtein *lpD, int startS, int startD, int numRes);
    /**
     *Modifies a protein loop/ fragment specified by <code>lpD</code> to overlap it with a specified 
     *loop/fragment specified by <code>lpS</code>.
     */
    static bool CopyBackbone(PProtein *lpS, PProtein *lpD);
    /**
     *Finds a new conformation of loop specified by <code>lp</code> after a random perturbation in its
     *null space. The magnitude of perturbation is specified by <code>pert_mag</code>. All Backbone DOFs of
     *the loop are perturbed.
     */
    static void RandomNullSpacePerturb(PProtein *lp, double pert_mag);
    /**
     *Finds a new conformation of loop specified by <code>lp</code> after a random perturbation in its null space. The magnitude of perturbation is specified by <code>pert_mag</code>. Only the DOFs specified by <code> Dofs</code> are perturbed.
     */
    static void RandomNullSpacePerturb(PProtein *lp, vector<vector<CDof> > Dofs, double pert_mag);
    /**
     *In many applications one might want to relax the closure constraint by a little bit and in the end the loop specified by <code>lp</code> has to be closed. This method finds a closed conformation closest in terms of RMSD distance to the open conformation. The pose of the closed conformation is defined by the three atoms <code>endPriorGoal</code>, <code>endGoal</code>, and <code>endNextGoal</code>.
     */
    static IKSolution CloseAlmostClosedLoop(PProtein *lp, Vector3 endPriorGoal,Vector3 endGoal, Vector3 endNextGoal);
    static int gsl_test();
  private:

};
//@package Conformation Analysis
/**
 * 
 *
 * Represents a cluster of "similar" loop conformations.
 */

class PCluster {
 public:

  /**
   * Creates a new <code>PCluster</code>.
   */

  PCluster();

  /**
   * Creates a new <code>PCluster</code> as a copy of
   * <code>other</code>.
   */

  PCluster(const PCluster &other);

  ~PCluster();

  /**
   * Adds the specified <code>PConformation</code>
   * to this cluster.
   */

  void AddConformation(PLightChain *toAdd);

  /**
   * Merges the conformation cluster <code>other</code>
   * with this <code>PCluster</code>, clearing the
   * contents of the <code>other</code> cluster.
   */

  void MergeWithCluster(PCluster &other);

  /**
   * Returns a "representative conformation" for this
   * cluster.
   */

  PLightChain *GetRepresentativeConformation() const;

  /**
   * Returns the number of conformations in this cluster.
   */

  int size() const { return m_confs->size(); }

  /**
   * Returns <code>true</code> if this cluster is empty,
   * <code>false</code> otherwise.
   */

  bool empty() const { return m_confs->empty(); }

  /**
   * Returns a new <code>list</code> of the conformations
   * in this cluster.
   */

  list<PLightChain *> getConformations() const {
    return *m_confs;
  }


 private:
  list<PLightChain *> *m_confs;
  
};

#endif
