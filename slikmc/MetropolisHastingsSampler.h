/*
 * BasicMetropolisSampling.h
 *
 *  Created on: Jul 24, 2012
 *      Author: Yajia
 */

#ifndef BASICMETROPOLISSAMPLING_H_
#define BASICMETROPOLISSAMPLING_H_
#include "BFactor.h"
#include "RamachandranPlot.h"
#include "PProtein.h"
#include "math/Gaussian.h"
using namespace Math;


/**
 * @brief The class MHSample is a sampler based on standard Metropolis-Hastings algorithm. This class can be used to generate free-end conformations.
 */
class MHSampler {
public:
	/**
	 * @brief Construct a standard Metropolis-Hastings sampler for specific chain.
	 */
	MHSampler(PProtein* chain);
	/**
	 * @brief Destructor
	 */
	virtual ~MHSampler();

	/**
	 * @brief Sample conformations of protein.
	 * @param time time duration for sampling
	 * @param radius perturbation radius in degrees
	 */
	void sample( const double time, const double radius);

	/**
	 * @brief Enable using B-factors as priors.
	 * @param chain a chain conformation with desired atom positions and B-factors.
	 */
	void enableBfactors( PProtein* chain = NULL);

	/**
	 * @brief Disable using B-factors as priors.
	 */
	void disableBfactors();
private:
	BFactor* bfactor;
	RamachandranPlot rplot;
	PProtein* chain;

	/**
	 * @brief Metropolis-Hastings step to decide whether to accept a proposal conformation
	 * @param P Probability density of initial conformation
	 * @param Q Proposal density of initial conformation to proposal conformation
	 * @param P_proposal Probability density of proposal conformation
	 * @param Q_proposal Proposal density of proposal conformation to initial conformation
	 * @return true: accept; false: reject
	 */
	bool MHStep( double P, double Q, double P_proposal, double Q_proposal);

	/**
	 * @brief Evaluate probability density given one conformation
	 * @param conformation to be evaluated
	 * @return probability density in logarithm
	 */
	double getP_log( PProtein* chain);

	bool use_BFactor;
};

#endif /* BASICMETROPOLISSAMPLING_H_ */
