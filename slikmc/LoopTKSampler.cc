/*
 * LoopTKSampler.cc
 *
 *  Created on: Jul 31, 2012
 *      Author: Yajia
 */

#include "LoopTKSampler.h"
#include "PSampMethods.h"
#include <stdio.h>
#include <iostream>
#include <queue>
#include <sstream>
using namespace std;

LoopTKSampler::LoopTKSampler(PProtein* protein) {
	this->chain = protein;
	this->bfactor = NULL;
	this->use_BFactor = false;
}

LoopTKSampler::~LoopTKSampler() {
	if( this->bfactor != NULL)
		delete this->bfactor;
}

void LoopTKSampler::sample( const double time_duration, const int s, const int e, const int num_conformation) {
	clock_t begin = clock();
	int num_generated = 0;
	vector<PProtein*> proteins;
	priority_queue< Protein_Score, vector<Protein_Score>, CompareProtein_Score > pq;

	PProtein* prev = this->chain;
	while( true) {
		cout << "Generating # " << num_generated << endl << flush;
		vector<PProtein*> curr_vector = PSampMethods::SeedSampleBackbone( prev, s, e);
		assert( curr_vector.size() == 1);

		PProtein* curr = curr_vector[0];

		//record every conformation
		if( num_conformation == -1) {
			stringstream ss;
			ss << num_generated;
			PDBIO::writeToFile( curr, "../pdbfiles_out/loopTK_" + ss.str() + ".pdb");
		}
		//record only top-scored conformations
		else {
			double score = this->evaluate_log( curr);
			Protein_Score ps;
			ps.protein = curr;
			ps.score = score;
			pq.push( ps);

			if( pq.size() > num_conformation) {
				Protein_Score low = pq.top();
				if( low.protein == curr) {
					delete curr;
				}
				else {
					prev = curr;
					delete low.protein;
					cout << "current score:" << score << "\t" << "low score:" << low.score << endl;
				}
				pq.pop();
				assert( pq.size() == num_conformation);
			}
		}
		num_generated += 1;
		cout << "done" << endl;

		clock_t current = clock();
		if( (current - begin) / 1000 > time_duration ) {
			cout << "Run out of time. Sampling stops" << endl;
			break;
		}
	}
	clock_t curr = clock();
	cout << "Duration:" << (curr - begin) / 1000.0 << endl;

	if(num_conformation != -1) {
		int i = 0;
		while (!pq.empty()) {
			stringstream ss;
			ss << i;
			PProtein* p = pq.top().protein;
			PDBIO::writeToFile(p, "../pdbfiles_out/loopTK_" + ss.str() + ".pdb");
			delete p;
			pq.pop();
			cout << "***Record # " << i << endl;
			i++;
		}
	}
}

double LoopTKSampler::evaluate_log(PProtein* chain) {
	double log_prob_rplot = 0;
	int start = 0;
	int end = chain->size() - 1;
	for( int i = start; i <= end; i++) {
		double temp = this->rplot.getResidueAngleProbability( chain, i);
		log_prob_rplot += log(temp);
	}

	double log_prob_bfactor = 0;
	if( this->use_BFactor) {
		if( chain != chain->getTopLevelChain())
			log_prob_bfactor = this->bfactor->evalAtomPositions( chain, chain->getTopLevelIndices().first, chain->getTopLevelIndices().second);
		else
			log_prob_bfactor = this->bfactor->evalAtomPositions( chain, 0, chain->size() - 1);
	}
	return log_prob_rplot + log_prob_bfactor;
}

bool CompareProtein_Score::operator ()(Protein_Score& t1, Protein_Score& t2) {
	if (t1.score > t2.score)
		return true;
	else
		return false;
}

PProtein* LoopTKSampler::perturb(PProtein* protein, const double time_duration, const int s, const int e) {
	clock_t begin = clock();

	PProtein* prev = protein;
	int iter = 0;
	while( true) {

		cout << iter++ << endl;

		vector<PProtein*> curr_vector = PSampMethods::SeedSampleBackbone( prev, s, e);
		delete prev;
		prev = curr_vector[0];

		clock_t current = clock();
		if( (current - begin) / 1000 > time_duration )
		{
			cout << "Run out of time. Sampling stops" << endl;
			break;
		}
	}
	return prev;
}

void LoopTKSampler::enableBFactors(PProtein* protein) {
	if( this->bfactor != NULL) {
		delete this->bfactor;
	}
	if( protein != NULL)
		this->bfactor = new BFactor( protein);
	else
		this->bfactor = new BFactor( this->chain);
	this->use_BFactor = true;
}

void LoopTKSampler::disableBFactors() {
	this->use_BFactor = false;
}
