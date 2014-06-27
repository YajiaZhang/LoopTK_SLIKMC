/*
 * SideChainRotater.cc
 *
 *  Created on: Jan 27, 2013
 *      Author: Yajia
 */

#include "SideChainRotater.h"
#include <PExtension.h>

SidechainRotater::SidechainRotater( PProtein* protein) {
	this->rotamer = new Rotamer( protein);
	this->rotamer->initSidechainDatabase( protein);
	return;
}

SidechainRotater::~SidechainRotater() {
	delete this->rotamer;
}

void SidechainRotater::rotateSidechain( PChain* chain, const vector<DihedralAngle>& angles_backbone) {

	int start = 0;
	int end = chain->size() - 1;

	//The first and last residue cannot decide dihedral angle, therefore, cannot decide the rotamer according to dihedral angle.
	int start_top = chain->getTopLevelIndices().first;
	int end_top = chain->getTopLevelIndices().second;
	int size_top = chain->getTopLevelChain()->size();
	if( start_top == 0){
		start = 1;
	}
	if( end_top == size_top - 1) {
		end = end - 1;
	}

	for( int i = start; i <= end; i++) {
		string name = chain->getResidue(i)->getName();
		if( name == "GLY" || name == "ALA") {
			//No rotamer for GLY and ALA
			continue;
		}
		PResidue* residue = chain->getResidue(i);

		double phi = angles_backbone[i].phi;
		double psi = angles_backbone[i].psi;
		vector<double> rotamer;
		int index = 0;
		this->rotamer->sample( name, phi, psi, index, rotamer);

		residue->applyRotamer( index, rotamer);
	}
	return;
}

double SidechainRotater::evalSidechain(PChain* chain) {
	int start = 0;
	int end = chain->size() - 1;

	int start_top = chain->getTopLevelIndices().first;
	int end_top = chain->getTopLevelIndices().second;
	int size_top = chain->getTopLevelChain()->size();
	if( start_top == 0){
		start = 1;
	}
	if( end_top == size_top - 1) {
		end = end - 1;
	}
	return this->rotamer->evalSidechain_log( chain, start, end);
}

void SidechainRotater::getSidechainAngles(PProtein* protein, vector<vector<double> >& angles) {

	for( int i = 0; i < protein->size(); i++) {
		vector<double> chis;
		protein->getResidue(i)->getChiAngles( chis);
		angles.push_back( chis);
	}
	return;
}
