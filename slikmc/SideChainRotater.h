/*
 * SidechainRotater.h
 *
 *  Created on: Jan 27, 2013
 *      Author: Yajia
 */

#ifndef SIDECHAINROTATER_H_
#define SIDECHAINROTATER_H_

#include <PProtein.h>
#include "Rotamer.h"

/**
 * @brief The class SidechainRotater is used to rotate side-chain conformations for a given chain and evaluate a given chain according to its side-chain structure.
 */
class SidechainRotater {
public:
	/**
	 * @brief Constructor
	 * @param chain the chain to be sampled, used for initialization of side-chain database.
	 */
	SidechainRotater(PProtein* chain);

	/**
	 * @brief Destructor
	 */
	virtual ~SidechainRotater();

	/**
	 * @brief Sample and rotate side-chains.
	 * @param chain the chain to be sampled.
	 * @param aBackbone backbone dihedral angles of the chain.
	 */
	void rotateSidechain( PChain* chain, const vector<DihedralAngle>& aBackbone);

	/**
	 * @brief Evaluate the chain according to side-chain conformations.
	 * @param chain the chain to be evaluated
	 * @return probability in logarithm
	 */
	double evalSidechain( PChain* chain);

	/**
	 * @brief Get all side-chain angles for a given chain.
	 * @param chain the given chain
	 * @param angles stores the side-chain angles
	 */
	static void getSidechainAngles(PProtein* chain, vector< vector<double> >& angles);
private:
	Rotamer* rotamer;
};

#endif /* SIDECHAINROTATER_H_ */
