/*
 * BFactor.cc
 *
 *  Created on: Jun 13, 2012
 *      Author: Yajia
 */

#include "BFactor.h"
#include "PConstants.h"
#include <cmath>
#include "Utility.h"

BFactor::BFactor(PProtein* chain) {
	//NOTE: This must be only constructed on top level chain.
	Debug::check( chain->getTopLevelChain() == chain, "The chain in initializing BFactors must be top level chain!");
	int size = chain->size();
	this->residue_name.resize( size);
	this->atom_name.resize( size);
	this->atom_pos.resize( size);
	this->atom_bfactors.resize( size);
	this->atom_variance.resize( size);

	for( int i = 0; i < size; i++) {
		PResidue* residue = chain->getResidue(i);

		this->residue_name[i] = residue->getName();
		PAtom* N = residue->getAtom( PID::N);
		PAtom* Ca = residue->getAtom( PID::C_ALPHA);
		PAtom* C = residue->getAtom( PID::C);

		this->atom_name[i].push_back( PID::N);
		this->atom_name[i].push_back( PID::C_ALPHA);
		this->atom_name[i].push_back( PID::C);
		this->atom_name[i].push_back( PID::O);

		HASH_MAP_STR(PAtom *)* atom_map = residue->getAtomMap();
		for( HASH_MAP_STR(PAtom *)::iterator atom_iter = atom_map->begin(); atom_iter != atom_map->end(); atom_iter++) {
			string name = atom_iter->first;
			if( name == PID::N || name == PID::C_ALPHA || name == PID::C || name == PID::O )
				continue;
			this->atom_name[i].push_back( name);
		}
	}

	for( int i = 0; i < size; i++) {
		PResidue* residue = chain->getResidue(i);
		for( int j = 0; j < this->atom_name[i].size(); j++) {
			PAtom* atom = residue->getAtom( this->atom_name[i][j]);
			this->atom_pos[i].push_back( atom->getPos());
			this->atom_bfactors[i].push_back( atom->getTempFactor());
			this->atom_variance[i].push_back( atom->getTempFactor() / (8 * PI * PI));
		}
	}
	return;
}

double BFactor::evalAtomPositions(PProtein* chain, const int s, const int e) {
	//First, we need to check whether the chain is the top chain or a subchain.
	//Make sure that the chain is a subchain of "this->protein" (only one level down) or is the root chain
//	assert( chain->IsSubChainOf( this->protein) || chain->getTopLevelChain() == NULL);
	Debug::check( e - s + 1 == chain->size(), "indices and subchain size doesn't match");

	int index_start = s;
	//logged probability density;
	double log_pd = 0;
	for( int i = 0; i < chain->size(); i++) {
		PResidue* residue = chain->getResidue(i);
		Vector3 N = residue->getAtomPosition( PID::N);
		Vector3 Ca = residue->getAtomPosition( PID::C_ALPHA);
		Vector3 C = residue->getAtomPosition( PID::C);

		//Get the global index in the vectors
		int index = index_start + i;

		double d_N = (N - this->atom_pos[index][0]).norm();
		double p_d_N = this->getProbDensity_log( 0, this->atom_variance[index][0], d_N);

		double d_Ca = (Ca - this->atom_pos[index][1]).norm();
		double p_d_Ca = this->getProbDensity_log( 0, this->atom_variance[index][1], d_Ca);

		double d_C = (C - this->atom_pos[index][2]).norm();
		double p_d_C = this->getProbDensity_log( 0, this->atom_variance[index][2], d_C);

		log_pd += ( p_d_N + p_d_Ca + p_d_C);
	}
	return log_pd;
}

double BFactor::getProbDensity_log( const Vector3& mean, const double variance, const Vector3& x) const{
	//NOTE: B-Factor is a VARIANCE!
	Vector3 d = x - mean;
	double c = 1 / (pow( 2 * PI, 1.5) * pow( variance * variance * variance, 0.5));
	double exponent = -0.5 * (d.x * d.x + d.y * d.y + d.z * d.z) / variance;

	double log_prob_density = log(c) + exponent;
	return log_prob_density;
}


double BFactor::getProbDensity_log( const double mean, const double variance, const double x)
{
	double c = 1 / sqrt( 2 * PI * variance);
	double exponent = -0.5 * ( mean - x) * (mean - x) / variance;
	return log(c) + exponent;
}

void BFactor::output(char* filename) {
	int size_residue = this->atom_pos.size();
	ofstream out;
	out.open( filename);
	if( !out.is_open()) {
		cerr << "Cannot open file " << filename << endl;
		abort();
	}
	for( int i = 0; i < size_residue; i++) {
		int size_atom = this->atom_pos[i].size();
		for( int j = 0; j < size_atom; j++) {
			out << this->residue_name[i] << "\t" << this->atom_name[i][j] << "\t" << this->atom_pos[i][j] << "\t" << this->atom_bfactors[i][j] << endl;
		}
	}
	out.flush();
	out.close();
	return;
}




