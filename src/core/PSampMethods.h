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

#include <iostream>
#include <time.h>
#include "PBasic.h"
#include "PChainNavigator.h"
#include "PConstants.h"
#include "PDataCollector.h"
#include "PExtension.h"
#include "PTools.h"
#include "POptimize.h"
#include "PPhiPsiDistribution.h"
#include <string>
#include <list>
#include <vector>
#include <stack>
#include <stdlib.h>
#include <math.h>
#include <map>

class SpaceRelationship {

  private:
	friend class PSampMethods;

        Real angle_c_n_ca;
        Real angle_ca_c_n;
        Real angle_o_c_n;
        Real length_c_n;
        Real length_ca_c;
        Real length_o_c;

        SpaceRelationship () {
                angle_c_n_ca = ANGLE_C_N_CA;
                angle_ca_c_n = ANGLE_CA_C_N;
                angle_o_c_n = ANGLE_O_C_N;
                length_c_n = LENGTH_C_N;
                length_ca_c = LENGTH_CA_C;
                length_o_c = LENGTH_O_C;
        }

        SpaceRelationship (PResidue *res1, PResidue *res2) {
                Vector3 ca1 = res1->getAtomPosition("CA");
                Vector3 c1 = res1->getAtomPosition("C");
                Vector3 o1 = res1->getAtomPosition("O");
                Vector3 n2 = res2->getAtomPosition("N");
                Vector3 ca2 = res2->getAtomPosition("CA");
                angle_c_n_ca = rad2deg*PMath::getAngle(c1-n2, ca2-n2);
                length_c_n = (c1-n2).norm();
                angle_ca_c_n = rad2deg*PMath::getAngle(ca1-c1,n2-c1);
                length_ca_c = (ca1-c1).norm();
                angle_o_c_n = rad2deg*PMath::getAngle(o1-c1,n2-c1);
                length_o_c = (o1-c1).norm();
        }

        void print() {
                cout << "angle_c_n_ca = " << angle_c_n_ca << endl;
                cout << "angle_ca_c_n = " << angle_ca_c_n << endl;
                cout << "angle_o_c_n = " << angle_o_c_n << endl;
                cout << "length_c_n = " << length_c_n << endl;
                cout << "length_ca_c = " << length_ca_c << endl;
                cout << "length_o_c = " << length_o_c << endl;
        }
}; // end class SapceRelationship


/**
 * This contains methods for generating samples using various strategies
 */
class PSampMethods {
  public:
	/**
	* Randomizes the <code>loop</code> and then closes it using Exact IK method
	* Coutsias et al. Sampled conformation is collision free or not
	* depending on <code>clash_free</code>.
	*/
	static IKSolution RandAndIKClose(PProtein *loop, bool clash_free);
	/**
	* Randomizes the <code>loop</code> and then closes it using Exact IK method
	* Coutsias et al. Sampled conformation is collision free or not
	* depending on <code>clash_free</code>. The conformation is output to
	* <code>pdbFileName</code>.
	*/
	
	static IKSolution RandAndIKClose(PProtein *loop,const string &pdbFileName, bool clash_free);
	/**
	* Generates <code>num_wanted</code> (number of) closed conformations starting from a closed conformation of a loop <code>loop</code>, by using various permutations of 6DOFs as needed by ExactIK solution method developed by Coutsias et al. 
	*/
	static IKSolutions PermuteIK(PProtein *loop, int num_wanted);
	
	/**
	* Generates <code>num_wanted</code> (number of) closed clash-free,backbone conformations by deforming a loop specified by residue indices <code>loopSid</code> and <code>loopEid</code> of protein chain specified by <code>orig_protein</code>.
	* The loop is deformed in a random direction in its null/tangent space and the magnitude of deformation is specified by <code>deform_mag</code>.
	* <code>orig_protein</code> stays unchanged. Side chains act like rigid bodies attached to the backbone.
	* New side chains can be placed using addSideChain method defined in this class.
	*/
	static vector<PProtein*> DeformSampleBackbone(PProtein *orig_protein, int loopSid, int loopEid, int num_wanted, double deform_mag);
	static PProtein* PSampMethods::DeformSampleBackbone(PProtein *protein, int loopSid, int loopEid, double deform_mag);

        /**
         * Generates <code>num_wanted</code> (number of)closed, clash-free backbone conformations of a loop specified by residue indices <code>loopSid</code> and
         * <code>loopEid</code> of protein chain specified by <code>protein</code>.
         * Note that, the residue index starts from 0. If <code>num_wanted</code> is not specified, the default is 1.
         * The <code>protein</code> is not changed at all in the function. The output is a vector of pointers to proteins containing the desired loop conformation.
         * The resulting closed, collsion-free loop backbones will only consist of N, C_alpha, C_beta, C and O atoms.
         * There will be no collision between any atom on the loop backbone with any other atom on the loop backbone or with any atom on the rest of the protein.
         * The phi and psi angles in the loop are sampled uniformly between 0 and 2 pi.
	 *
         */
	static vector<PProtein*> SeedSampleBackbone (PProtein *protein, int loopSid, int loopEid, int num_wanted=1);
//	static PProtein* SeedSampleBackbone (PProtein *protein, int loopSid, int loopEid);

	/**
	 * Similar to the above function. In addition, the loops are with side chains, which are also not in collision. Side-chains are added using SCWRL3. User needs to provide the path of the executable of SCWRL3 as <code>scwrl3_path</code>.
	**/ 
	static vector<PProtein*> SeedSampleBackboneWithSidechain (PProtein *protein, int loopSid, int loopEid, string scwrl3_path, int num_wanted=1);

	/**
	 * Similar to the above function. However, rather than returning the entire protein, only the loop portion is returned.
	**/
	static vector<PProtein*> SeedSampleBackboneWithSidechainLoopOnly (PProtein *original_protein, int loopSid, int loopEid, string scwrl3_path, int num_wanted=1);


        /**
         * Same as the other SeedSampleBackbone method, but the phi and psi angles are sampled according a distribution. If the distribution map <code>distri_map</code> only has one element, then this distribution will be applied on all amino acids. Otherwise, there should be 20 distributions in the map. Each distribution corresponds to one amino acid, and the corresponding name of a distribution should be the 3-letter amino acid name in all capital letters.
	 */
         /* Usage example:
         *      map<string,PPhiPsiDistribution> Map = PPhiPsiDistribution::generateRamachandran();
         *      PProtein *protein = PDBIO::readFromFile("pdbfiles/135L.pdb");
         *      int loopSid=64, loopEid=72;
         *      vector<PProtein*> newps = PSampMethods::SeedSampleBackbone(protein,loopSid,loopEid,Map);
         *      PProtein *loop = new PProtein(newps[0],loopSid,loopEid);
         *      // Write the loop to a pdb file
         *      PDBIO::writeToFile(loop,"pdbfiles/135L_seed.pdb");
         *      // Add side chain
         *      PSampMethods::addSidechain("pdbfiles/135L_seed.pdb","pdbfiles/135L_NoLoop.pdb","/home/scwrl3_lin","pdbfiles/135L_seed_withSC.pdb");
         *      // Merge the loop with side chain and the protein
         *      PProtein *new_protein = PSampMethods::MergeProtein(loop,protein,loopSid);
         *      PDBIO::writeToFile(new_protein,"pdbfiles/135_new.pdb");
         */
	static vector<PProtein*> SeedSampleBackbone (PProtein *protein, int loopSid, int loopEid, map<string,PPhiPsiDistribution> &distri_map, int num_wanted=1);

        /**
         * Similar to the above function. In addition, the loops are with side chains, which are also not 
	 * in collision. Side-chains are added using SCWRL3. User needs to provide the path of the 
	 * executable of SCWRL3 as <code>scwrl3_path</code>.
        **/
	static vector<PProtein*> SeedSampleBackboneWithSidechain (PProtein *protein, int loopSid, int loopEid, map<string,PPhiPsiDistribution> &distri_map, string scwrl3_path, int num_wanted);

        /**
         * Similar to the above function. However, rather than returning the entire protein, only the 
	 * loop portion is returned.
        **/
	static vector<PProtein*> SeedSampleBackboneWithSidechainLoopOnly (PProtein* original_protein, int loopSid, int loopEid, map<string,PPhiPsiDistribution> &distri_map, string scwrl3_path, int num_wanted);

	/**
	 * Merges the <code>loop</code> with the <code>original_protein</code>. The 0-th residue in the loop becomes
	 * the <code>startRid</code>-th residue in the resulting protein.
	 * Neither <code>loop</code> or <code>original_protein</code> is changed in the function.
	 * Users are responsible to delete the returned protein if they don't need it any more.
	 */
	static PProtein* MergeProtein (PProtein *loop, PProtein* original_protein, int startRid);

	/**
	 * Similar to the above function. However, the 0-th residue in the loop becomes the residue with PDB ID 
	 * <code>startPdbId</code> in the protein.
	**/
	static PProtein* MergeProteinByPdbId (PProtein *loop, PProtein* original_protein, int startPdbId);

	/**
	 * Add side chains to a loop using SCWRL3 with option -i, -o and -f.
	 * For the details of SCWRL3 and its usage, please go to http://dunbrack.fccc.edu/SCWRL3.php. 
	 * The loop is specified in a pdb file <code>loopFile</code> (option -i).
	 * A boundary is specified in a pdb file <code>boundaryFile</code> (option -f).
	 * Please specify the FULL path of your scwrl3 program at <code>scwrl3_path</code>.
	 * The loop with side chain will be written into a pdb format file <code>outLoopFile</code>.
	 * Note that, in the output loop file <code>outLoopFile</code>, columns after the 3D coordinates
	 * are not meaningful. 
	 */
	static void addSidechain (string loopFile, string boundaryFile, string scwrl3_path, string outLoopFile);

	
	/**
	 * Add side chains to a portion of a protein using SCWRL3 with option -i, -o and -s.
	 * For the details of SCWRL3 and its usage, please go to http://dunbrack.fccc.edu/SCWRL3.php. 
	 * The protein is specified in PDB format in file <code>protein_input</code>.
	 * The protion of the protein from the residue <code>addStart</code> (as in the PDB file)
	 * to the residue <code>addEnd</code> is to be placed side chains, while the rest
	 * is to serve as the boundary.
	 * Please specify the FULL path of your scwrl3 program at <code>scwrl3_path</code>.
	 * The entire protein with side-chain-placement in the portion will be written into the file
	 * <code>protein_output</code>.
	 *  Note that, in the output loop file <code>outLoopFile</code>, columns after the 3D coordinates
         * are not meaningful. 
	 */
	static void AddSidechain (string protein_input, int addStart, int addEnd, string scwrl3_path, string protein_output);

        static vector<PProtein*> SeedSampleBackboneLoopOnly (PProtein* original_protein, int loopSid, int loopEid, int num_wanted=1);

        static vector<PProtein*> SeedSampleBackboneLoopOnly (PProtein* original_protein, int loopSid, int loopEid, map<string,PPhiPsiDistribution> &distri_map, int num_wanted=1);

	/**
	 * Fill in a missing loop in protein <code>original_p</code> from residue ID as in 
	 * the PDB file <code>start_pdb_id</code> to <code>end_pdb_id</code> with a randomly
	 * generated and closed conformation. This conformation is not guaranteed to be 
	 * collision-free. The <code>original_p</code> is not modified. This function is 
	 * useful when you want to do seed sampling on a missing loop. Invoke this function
	 * first, and then invoke one of the seed sampling functions.
	 */
	static PProtein* FillMissingLoop (PProtein *original_p, int start_pdb_id, int end_pdb_id, vector<string> loop_seq);


    private:
        static vector<Real>* generateForwardOpenLoop(PProtein *loopEntire, int collisionFreeResNum);
        static vector<Real>* generateForwardOpenLoop (PProtein* loopEntire, int collisionFreeResNum, map<string,PPhiPsiDistribution> &map_distri, vector<string> &aa_names);
        static vector<Real>* generateBackwardOpenLoop (PProtein* loop);
        static vector<Real>* generateBackwardOpenLoop (PProtein* loop, map<string,PPhiPsiDistribution> &map_distri, vector<string> &aa_names);
        static Vector3 computePreCpos (PResidue *res, Real angle, Real bondLength);
        static void computeGoal (PResidue *res, SpaceRelationship *sr, Vector3 *endPriorG, Vector3 *endG, Vector3 *endNextG);
        static bool IKClose (PProtein* move_loop, Vector3 endPriorG, Vector3 endG, Vector3 endNextG);
        static void UpdateProtein(PProtein *p, vector<Real> *positiveAngles, vector<Real> *negativeAngles, bool forwardDir);
        static void UpdateProtein(PProtein *p, vector<Real> *angles, bool forwardDir);
        static double computeMaxLength (int proteinSize);


};

