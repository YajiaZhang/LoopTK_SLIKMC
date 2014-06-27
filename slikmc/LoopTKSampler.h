/*
 * LoopTKSampler.h
 *
 *  Created on: Jul 31, 2012
 *      Author: Yajia
 */

#ifndef LOOPTKSAMPLER_H_
#define LOOPTKSAMPLER_H_
#include "PProtein.h"
#include "RamachandranPlot.h"
#include "BFactor.h"

/**
 * @brief An auxiliary class for class LoopTKSampler. This class is a data structure for storing a chain conformation and its score.
 */
class Protein_Score{
public:
	PProtein* protein;
	double score;
};

/**
 * @brief An auxiliary class for class LoopTKSampler. This class is used as a comparator to compare two conformations according to their scores.
 */
class CompareProtein_Score {
public:
    bool operator()(Protein_Score& t1, Protein_Score& t2);
};

/**
 * @brief LoopTK Configuration Sampling technique. This sampler can be used to sample close-loop chain/subchain conformations.
 */
class LoopTKSampler {
public:
	/**
	 * @brief Construct a LoopTK configuration sampler for specific chain.
	 */
	LoopTKSampler(PProtein* chain);
	/**
	 * @brief Destructor
	 */
	virtual ~LoopTKSampler();

	/**
	 * @brief Sample conformations of sub-loop from residue s to residue e.
	 * @param time time duration to sample
	 * @param s index of starting residue
	 * @param e index of ending residue
	 * @param num number of top-scores conformations to be saved, by default, save all conformations.
	 */
	void sample( const double time, const int s, const int e, const int num = -1);

	/**
	 * @brief Evaluate the chain according to priors.
	 * @param chain chain to be evaluated
	 * @return probability density in logarithm
	 */
	double evaluate_log( PProtein* chain);

	/**
	 * @brief Perturb a chain segment for some time.
	 * @param chain the chain to be perturbed
	 * @param time time duration to perturb
	 * @param s index of starting residue
	 * @param e index of ending residue
	 * @return a perturbed chain
	 */
	static PProtein* perturb( PProtein* chain, const double time, const int s, const int e);

	/**
	 * @brief Enable using B-factors as priors.
	 * @param protein a protein conformation with desired atom positions
	 */
	void enableBFactors( PProtein* protein = NULL);

	/**
	 * @brief Disable using B-factors as priors.
	 */
	void disableBFactors();
private:
	BFactor* bfactor;
	RamachandranPlot rplot;
	PProtein* chain;

	bool use_BFactor;
};

#endif /* LOOPTKSAMPLER_H_ */
