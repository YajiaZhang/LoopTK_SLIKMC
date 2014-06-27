/*
 * SLIKMC.h
 *
 *  Created on: Jun 10, 2012
 *      Author: Yajia
 */

/**
\mainpage Overview

Sub-Loop Inverse Kinematics Monte Carlo (SLIKMC) is a novel technique for sampling configurations of a kinematic chain according to a
specified probability density while accounting for loop closure constraints.
This method is an integration of two Markov chain Monte Carlo method, blocked Gibbs sampling and Metropolis-Hastings algorithm and it samples chain configurations in an unbiased manner.
SLIKMC is applicable in sampling close-loop/free-end kinematic chains, protein loops in 3D space.
The samples generated by SLIKMC are verified to have higher quality than ad-hoc fashion loop construction methods and sampling speed is proven to outrun methods based on discrete search.

SLIKMC is developed by Yajia Zhang and Kris Hauser at Intelligent Motion Lab in Indiana University Bloomington.
SLIKMC is implemented using the software package Protein Loop Kinematic Toolkit (LoopTK, https://simtk.org/home/looptk) and released as an extension version of LoopTK.

<H2> Download and Installation </H2>
<B>Download:</B> <A HREF = "../LoopTK_release.zip" > package </A>

<B>System requirement:</B> Cygwin/Linux, Cygwin is recommended.

<B>Required software libraries:</B>

1.GNU Scientific Library (GSL) http://www.gnu.org/software/gsl/

2.Mesa 3D Graphics Library (Mesa) version 8.0.2 is recommended. http://mesa3d.sourceforge.net/

3.GLUI User Interface Library http://glui.sourceforge.net/

<B>Installation</B>

1.Install LoopTK package.

For installing LoopTK, please refer to the instruction here: http://ai.stanford.edu/looptk/download.html

2.Install SLIKMC.

Simply go to LoopTK/slikmc and "make". An executable file will be generated with name slikmc.

<H2>Quick-start tutorial </H2>
Please see

<H2>Contact </H2>

For comments or bug issue, please contact the software maintainer <B>Yajia Zhang</B> (zhangyaj@indiana.edu).
*/

#ifndef SLIKMC_H_
#define SLIKMC_H_

#include "PProtein.h"
#include <vector.h>
#include "RamachandranPlot.h"
#include "BFactor.h"
#include "SideChainRotater.h"
#include "Prior.h"

/**
 * @brief Sub-Loop Inverse Kinematic Markov Chain (SLIKMC) sampler. Support chain/subchain close-loop sampling, free-end sampling. Side-chain sampling is also supported.
 */
class SLIKMCSampler {
public:
	/**
	 * @brief Construct a SLIKMC sampler for specific chain
	 */
	SLIKMCSampler( PProtein* chain);
	/**
	 * @brief Destructor
	 */
	virtual ~SLIKMCSampler();

	/**
	 * @brief Sample conformations of sub-loop from residue s to residue e.
	 * @param time time duration for sampling
	 * @param s index of starting residue
	 * @param e index of ending residue
	 */
	void sample( const double time, const int s, const int e);

	/**
	 * @brief Enable using B-factors as priors.
	 * @param chain a chain conformation with desired atom positions and B-factors
	 */
	void enableBfactors( PProtein* chain = NULL);

	/**
	 * @brief Disable using B-factors as priors.
	 */
	void disableBfactors();

	/**
	 * @brief Enable side-chain conformation sampling.
	 */
	void enableSidechain();

	/**
	 * @brief Disable side-chain in sampling.
	 */
	void disableSidechain();

	/**
	 * @brief Enable free-end sampling (terminal atoms will stay at fixed positions).
	 */
	void enableFreeEnd();

	/**
	 * @brief Disable free-end sampling.
	 */
	void disableFreeEnd();

	/**
	 * @brief Enable saving generated conformations.
	 * @param skip number of conformations to skip after last saving conformation.
	 */
	void enableLog(const int skip = 0);

	/**
	 * @brief Disable saving generated conformations.
	 */
	void disableLog();

	/**
	 * @brief Enable steric clash checking for samples.
	 */
	void enableCollisionChecking();

	/**
	 * @brief Disable steric clash checking for samples.
	 */
	void disableCollisionChecking();

	/**
	 * @brief Enable using Ramachandran plot as prior.
	 */
	void enableRamachandran();

	/**
	 * @brief Disable using Ramachandran plot as prior.
	 */
	void disableRamachandran();

	/**
	 * @brief Print current settings for sampling.
	 */
	void dispSettings();

	/**
	 * @brief Enable custom defined priors.
	 */
	void enableCustomPriors();

	/**
	 * @brief Disable custom defined priors.
	 */
	void disableCustomPriors();

	/**
	 * @brief Add a custom defined prior.
	 */
	void addCustomPrior( Prior& prior);

private:
	/**
	 * @brief Metropolis-Hastings step to decide whether to accept a proposal block.
	 * @param P Probability density of initial block
	 * @param Q Proposal density of initial block to proposal block
	 * @param P_proposal Probability density of proposal block
	 * @param Q_proposal Proposal density of proposal block to initial block
	 * @return true: accept; false: reject
	 */
	bool MHStep( double P, double Q, double P_proposal, double Q_proposal);

	/**
	 * @brief Evaluate probability density of one block or sub-chain.
	 * @return probability density in logarithm
	 */
	double getP_log( PProtein* chain);

	/**
	 * @brief Evaluate proposal density given one block or sub-chain. Initial block and proposal block are assumed to be independent.
	 * @param num_solutions number of IK solution for the block or sub-chain
	 * @param status an indicator if matrix is not invertible in calculating metric tensor
	 * @return proposal density in logarithm
	 */
	double getQ_log( PProtein* chain, const int num_solutions, int& status);

	/**
	 * @brief Calculate the metric tensor for one block.
	 * @param protein a block to be evaluated
	 * @param status an indicator if matrix is not invertible in calculating metric tensor
	 * @return
	 */
	double getMetricTensor_log(PProtein* protein, int& status);

	/**
	 * @brief The protein chain to sample
	 */
	PProtein* protein;

	/**
	 * @brief Collection of 4-residue size blocks
	 */
	vector<PProtein*> subchains;

	/**
	 * @brief Maximum number of dihedral angles we try for the first residue in the subchain in case that IK cannot find a solution.
	 */
	static const int MAX_IK_SAMPLE = 100;

	/**
	 * @brief Maximum number of Metropolis Hasting rejects before we giving up.
	 */
	static const int MAX_METROPOLIS_REJECT = 1;

	/**
	 * @brief Maximum number of collision detected before we giving up.
	 */
	static const int MAX_COLLISION_DETECT = 1;

	BFactor* bfactor;
	RamachandranPlot* rplot;

	SidechainRotater* scRotater;

	bool use_BFactor;
	bool use_Rotamer;
	bool freeEnd;
	bool use_colChecking;
	bool use_RPlot;
	bool use_customPrior;

	bool init_Rotamer;
	bool logFile;
	int skipLength;

	vector<Prior*> priors;
};

const double EPSILON = 0.0000001;
#endif /* SLIKMC_H_ */
