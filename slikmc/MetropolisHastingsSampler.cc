/*
 * BasicMetropolisSampling.cc
 *
 *  Created on: Jul 24, 2012
 *      Author: Yajia
 */

#include "MetropolisHastingsSampler.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <climits>
#include "math/MatrixTemplate.h"
#include "PConstants.h"
#include "PExtension.h"
using namespace std;
using namespace Math;

MHSampler::MHSampler(PProtein* protein) {
	this->chain = protein;
	this->bfactor = NULL;
	this->use_BFactor = false;
	return;
}

MHSampler::~MHSampler() {
	if( this->bfactor != NULL)
		delete this->bfactor;
}

void MHSampler::sample( const double time_duration, const double std_perturb) {
	//In order to adapt to LoopTK tool functions, vector should start from 1.
	int size_dihedral = this->chain->size() * 2 + 1;

	int size = size_dihedral;
	VectorTemplate<double> mean( size);
	mean.setZero();
	MatrixTemplate<double> variance( size, size);
	variance.setZero();
	for( int i = 0; i < size; i++) {
		variance(i, i) = std_perturb * std_perturb;
	}

	Gaussian<double> gaussian;
	gaussian.resize( size);
	gaussian.setMean( mean);
	gaussian.setCovariance( variance);

	int size_residue = this->chain->size();
	int size_rotation = size_residue * 2 + 1;

	clock_t begin = clock();
	int count_success = 0;
	int count_total = 0;

	while( true) {
		PChainState* state = this->chain->saveChainState();
		double P = this->getP_log( this->chain);
		double Q = 1;

		bool success = false;
		cout << "Sampling # " << count_total;
		VectorTemplate<double> perturb( size_rotation);
		gaussian.generate( perturb);

		for( int i = 0; i < size_rotation - 1; i++) {
			this->chain->RotateChain_noGridUpdate("backbone", i, forward, perturb[i]);
		}
//		this->protein->RotateChain_noGridUpdate( "backbone", 3, backward, perturb[ size_rotation - 1]);

		double P_proposal = this->getP_log( this->chain);
		double Q_proposal = 1;
		if( this->MHStep( P, Q, P_proposal, Q_proposal) == true) {
			this->chain->updateAtomsGrid();
			if( this->chain->InAnyCollision() == false) {
				stringstream ss;
				ss << count_success;
				PDBIO::writeToFile(chain, "../pdbfiles_out/mh_" + ss.str() + ".pdb");
				count_success += 1;
				success = true;
			}
		}
		if( success == false) {
			this->chain->restoreChainState_noGridUpdate( state);
			cout << "\tReject." << endl;
		}
		else
			cout << "\tAccept." << endl;

		count_total += 1;
		clock_t curr = clock();
		if( (curr - begin) / 1000 > time_duration ) {
			cout << "Run out of time. Sampling stops" << endl;
			break;
		}
	}

	cout << "Duration: " << time_duration << endl;
	cout << " Total sampling: " << count_total << endl;
	cout << " Accepted: " << count_success << endl;
	cout << " Accept Ratio: " << count_success / (double)count_total << endl;
}

bool MHSampler::MHStep( double P, double Q, double P_proposal, double Q_proposal)
{
	//Be careful that they are the logged probability.
	double ratio_log = P_proposal + Q - P -  Q_proposal;
	if( ratio_log >= 0) {
		return true;
	}
	else {
		double ratio = exp(ratio_log);
		double temp = (double)rand()/(double)RAND_MAX;
		if( temp < ratio)
			return true;
		else
			return false;
	}
}

void MHSampler::enableBfactors(PProtein* protein) {
	if( this->bfactor != NULL) {
		delete this->bfactor;
	}
	if( protein != NULL) {
		assert( this->chain->size() == protein->size());
		this->bfactor = new BFactor( protein);
	}
	else
		this->bfactor = new BFactor( this->chain);
	this->use_BFactor = true;
	return;
}

void MHSampler::disableBfactors() {
	this->use_BFactor = false;
	return;
}

double MHSampler::getP_log( PProtein* chain)
{
	double log_prob_rplot = 0;
	int start = 0;
	int end = chain->size() - 1;
	for( int i = start; i <= end; i++) {
		double temp = this->rplot.getResidueAngleProbability( chain, i);
		log_prob_rplot += log(temp);
	}

	double log_prob_bfactor = 0;
	if( this->use_BFactor)
		log_prob_bfactor = this->bfactor->evalAtomPositions( chain, 0, chain->size() - 1);
	return log_prob_rplot + log_prob_bfactor;
}




