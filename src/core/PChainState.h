/*
 * PChainState.h
 *
 *  Created on: Jun 4, 2012
 *      Author: Yajia
 */

#ifndef PCHAINSTATE_H_
#define PCHAINSTATE_H_

#include <vector.h>
#include <math3d/primitives.h>
using namespace Math3D;

class PChainState {
	friend class PChain;
public:
	PChainState();
	virtual ~PChainState();
	void clean();
private:
	//The local residue index in a specific chain
	int index_start;
	int index_end;
	vector<vector<Vector3> > atoms_pos;
};

#endif /* PCHAINSTATE_H_ */
