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

#ifndef PPHIPSIDISTRIBUTION_H
#define PPHIPSIDISTRIBUTION_H

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <map>

using namespace std;

/**
 * The PPhiPsiDistribution class is for sampling backbone conformations according to a certain distribution.
 * Currenlty, the seed sampling methods can take in only 1 distribution (for all amino acids) or 
 * 20 distributions (one for each amino acid, such as a Ramachandran distribution). In the later case, the name of the 
 * distribution must be 3-letter amino acid name in all capital letters. In the previous case, name doesn't matter, and 
 * set to "anything" by default.
 */
class PPhiPsiDistribution {
  public:

  PPhiPsiDistribution() { }
	/**
	 * Construct a PPhiPsiDistribution object by providing the distribution <code>v</code> and a corresponding <code>name</code>.
         * In v, the 1st dimension is about Phi and the 2nd dimension is about Psi
         * i.e., v[i][j] is the probability of taking the Phi in the ith interval and Psi in the jth interval.
         * Important: 
         * 1. All inner vectors must have the same size.
         * 2. All entries should sum up to 1.
         * 3. v[i][0] to v[i][n-1] store the probability of Psi in 360/n degree intervals from -180 to +179.
         *    v[0][j] to v[n-1][j] store the probability of Phi in 360/n degree intervals from -180 to +179.
	 * For example:
	 *	vector<double> row1, row2;
	 *	row1.push_back(0.1); row1.push_back(0.2); row1.push_back(0.3);
	 *	row2.push_back(0.2); row2.push_back(0.1); row2.push_back(0.1);
	 * 	vector< vector<double> > phipsi;
	 * 	phipsi.push_back(row1); phipsi.push_back(row2);
	 * This means: the phi angle has two intervals: [-180,0) and [0,+179], 
	 *             and psi angle has three intervals: [-180,-60), [-60,+60) and [+60,179].
	 *             The probabilty of a residue taking phi angle in [-180,0) and psi angle in [+60,+179] is 0.3.
	 */
	PPhiPsiDistribution(const vector< vector<double> > &v, const string &name="anything");
	
	/*
	 * Prints a human-readable representation of this
	 * PPhiPsiDistribution to the specified output stream.
   * Useful for debugging.
	 */
	void print(ostream &out) const;

	/*
	 * Given a Phi value <code>phi</code> and Psi value <code>psi</code>, return the probabilty of having them according to 
	 * this distribution.
	 */
	double getProbPhiPsi (double phi, double psi) const;

	/*
	 * Given a Psi value <code>psiValue</code>, return a distribution of Phi.
	 * The ith element in the output is the probability of phi in the ith interval given the psi value.
	 */
	vector<double> getPhiDistribution (double psiValue) const;

	/*
	 * Given a Phi value <code>phiValue</code>, return a distribution of Psi.
	 * The ith element in the output is the probability of psi in the ith interval given the psi value.
	 */
	vector<double> getPsiDistribution (double phiValue) const;

	/*
	 * Unconditional probability of Phi.
	 * The ith element in the output is the probability of phi in the ith interval if we do not have any priori knowledge on psi.
	 */
	vector<double> getPhiDistribution() const { return PhiDistribution; }

	/*
	 * Unconditional probability of psi.
	 * The ith element in the output is the probability of psi in the ith interval if we do not have any priori knowledge on phi.
	 */
	vector<double> getPsiDistribution() const { return PsiDistribution; }

	/*
 	 * Return true if the distribution is empty; false otherwise.
	 */
	bool isEmpty() const { return distribution.empty(); }

  /**
   * Returns the number of phi intervals in this distribution.
   */
  int numPhiIntervals() const { return PhiIntervalNum; }

  /**
   * Returns the number of psi intervals in this distribution.
   */
  int numPsiIntervals() const { return PsiIntervalNum; }

	/*
	 * Create a Ramachandran distribution using the distribution data in our database.
	 */
	static map<string,PPhiPsiDistribution> generateRamachandran();

  private:
    string AA_Name;
    vector< vector<double> > distribution;
    vector<double> PhiDistribution;
    vector<double> PsiDistribution;
    int PhiIntervalNum;
    int PsiIntervalNum;
    double PhiIntervalSize;
    double PsiIntervalSize;

    int getPhiIntervalIdx (double phiValue) const;
    int getPsiIntervalIdx (double psiValue) const;
};

#endif
