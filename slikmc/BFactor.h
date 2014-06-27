/*
 * BFactor.h
 *
 *  Created on: Jun 13, 2012
 *      Author: Yajia
 */

#ifndef BFACTOR_H_
#define BFACTOR_H_
#include "PProtein.h"
#include <vector.h>
using namespace Math3D;
using namespace std;

/**
 * @brief The class BFactor is used to initialize and store B-factors database, evaluate the probability densities of conformations according to B-factors.
 */
class BFactor {
public:
	/**
	 * @brief Construct a B-factors database for a specific chain. The chain must be the top level chain in the sampling problem.
	 */
	BFactor(PProtein* protein);

	/**
	 * @brief Evaluate a chain (or its subchain) according to B-factors.
	 * @param chain the chain to be evaluated
	 * @param s index of starting residue
	 * @param e index of ending residue
	 * @return products of probability densities of all backbone atoms in logarithm
	 */
	double evalAtomPositions( PProtein* chain, const int s, const int e);

	/**
	 * @brief Output atom positions and their B-factors to a file
	 * @param filename output filename
	 */
	void output( char* filename);
private:
	vector< string> residue_name;
	vector< vector<string> > atom_name;
	vector< vector<Vector3> > atom_pos;
	vector< vector<double> > atom_bfactors;
	vector< vector<double> > atom_variance;

	/**
	 * @brief Evaluate probability density of an atom position given desired atom position and B-factor.
	 * @param mean desired atom position
	 * @param variance variance of each dimension
	 * @param position the atom position to be evaluated
	 * @return probability density of the atom position in logarithm
	 */
	double getProbDensity_log( const Vector3& mean, const double variance, const Vector3& position) const;

	/**
	 * @brief Get probability density of one value in a given normal distribution
	 * @param mean mean of the normal distribution
	 * @param variance variance of the normal distribution
	 * @param x value to be evaluated
	 * @return probability density of x in logarithm
	 */
	double getProbDensity_log( const double mean, const double variance, const double x);
};

#endif /* BFACTOR_H_ */
