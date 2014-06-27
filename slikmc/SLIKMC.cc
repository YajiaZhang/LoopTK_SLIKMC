/*
 * SLIKMC.cc
 *
 *  Created on: Jun 10, 2012
 *      Author: Yajia
 */

#include "SLIKMC.h"
#include "PBasic.h"
#include "PChain.h"
#include "PChainNavigator.h"
#include "PDataCollector.h"
#include "PExtension.h"
#include "PHydrogenBondTracker.h"
#include "PTools.h"
#include "PConfSpaceNavigator.h"
#include "PIKAlgorithms.h"
#include "PChainState.h"
#include <iostream>
#include <algorithm>
#include <time.h>
#include "Utility.h"
#include "math/MatrixTemplate.h"
#include "math/Gaussian.h"
#include <Math/Matrix.h>
using namespace Math;

SLIKMCSampler::SLIKMCSampler(PProtein* protein) {
	this->protein = protein;
	//every intermediate block overlaps three residues with succeeding block
	int num_subchains = this->protein->size() - 3;
	for (int i = 0; i < num_subchains; i++)
	{
		int start = i;
		int end = i + 3;
		PProtein* chain = new PProtein(protein, start, end);
		this->subchains.push_back(chain);
	}

	this->init_Rotamer = false;
	this->use_BFactor = false;
	this->use_Rotamer = false;
	this->use_colChecking = false;
	this->use_RPlot = false;
	this->use_customPrior = false;

	this->rplot = new RamachandranPlot();
	this->freeEnd = false;

	this->bfactor = NULL;
	this->scRotater = NULL;

	//For logging conformations
	this->logFile = false;
	this->skipLength = 0;
	return;
}

SLIKMCSampler::~SLIKMCSampler() {
	if( this->bfactor != NULL) delete this->bfactor;
	if( this->scRotater != NULL) delete this->scRotater;
	delete this->rplot;
	return;
}

void SLIKMCSampler::sample( const double time, const int s, const int e) {

	Debug::check( s >= 0 && e >= s + 3 && e < this->protein->size(), "Sth wrong with the starting index and ending index");
	clock_t begin = clock();
	int size_protein = this->protein->size();
	int num_subchains = size_protein - 3;

	int stat_distinct = 0;
	int stat_conformation = 0;

	int s_chain = s;
	int e_chain = e - 3;
	int i = 0;
	while( true) {
		cout << "Start sampling conformation # " << i << ":" << endl;
		bool changed = false;
		for( int j = s_chain; j <= e_chain; j++) {
			PProtein* subchain = this->subchains[j];
			subchain->attachResidues(); 								//NOTE:This is necessary!
			PChainState* state_subchain = subchain->saveChainState(); 	//Save the s

			Vector3 endPriorG = subchain->getAtomAtRes(PID::C_ALPHA, 3)->getPos();;
			Vector3 endG = subchain->getAtomAtRes(PID::C, 3)->getPos();
			Vector3 endNextG = subchain->getAtomAtRes(PID::O, 3)->getPos();
			int index_to_use[3] = { 1, 2, 3};
			/*
			 * calculate the importance ratio for the initial block.
			 */
			double P = this->getP_log( subchain);
			IKSolutions iks_initial;
			int n = PExactIKSolver::FindSolutions( subchain, index_to_use, &endPriorG, &endG, &endNextG, iks_initial);
			int status = 0; //Record whether the metric tensor part can be calculated. -1: matrix cannot be inverted; 0: okay
			double Q = this->getQ_log( subchain, n, status);
			if( !(status == 0)) {
				cout << "Calculating matrix failed in the first place" << endl;
				//If calculating metric tensor not successful, then, move on to the next sub-chain.
				continue;
			}

			int	num_MH_reject = 0;
			int num_collision_reject = 0;

			while( (num_MH_reject < MAX_METROPOLIS_REJECT) && (num_collision_reject < MAX_COLLISION_DETECT)) {
				int n_proposal = 0;

				int num_IK_fail = 0;
				bool IK_success = false;
				while( num_IK_fail < this->MAX_IK_SAMPLE) {
					/* The (phi, psi) pair for the first residue is sampled from Ramachandran plot.
					 * If the residue is not the first in the whole protein chain, then, both phi, psi are defined.
					 * Therefore, sample (phi, psi) and rotate the two bonds.
					 * Otherwise, phi will not have definition, then, we just rotate the first bond a little bit.
					 */

					if( this->use_RPlot) {
						PResidue* residue = subchain->getResidue(0);
						PResidue* residue_next = subchain->getResidue(1);
						DihedralAngle* da_goal = this->rplot->getRandomDihedralAngle(residue->getName(), residue_next->getName());
						DihedralAngle* da_curr = subchain->getDihedralAngleAtResidue(0);
						if( j != 0) {
							subchain->RotateChain_noGridUpdate("backbone", 0, forward, da_curr->phi - da_goal->phi);
						}
						else {
							//The first residue in the chain
							double range = 60;
							double angle = (rand() % 100) / 100.0 * range - range / 2;
							subchain->RotateChain_noGridUpdate("backbone", 0, forward, angle);
						}
						subchain->RotateChain_noGridUpdate("backbone", 1, forward, da_curr->psi - da_goal->psi);
						delete da_goal;
						delete da_curr;
					}
					else {
						double phi_change = Random::nextNormal( 0, 10);
						double psi_change = Random::nextNormal( 0, 10);

						subchain->RotateChain_noGridUpdate("backbone", 0, forward, phi_change);
						subchain->RotateChain_noGridUpdate("backbone", 1, forward, psi_change);
					}

					/* Use analytical IK to close the sub-loop using the rest 6 DOFs
					 */
					IKSolutions iks;
					// Here, call a modified IK solver which returns only one solution and the number of possible solutions.
					n_proposal = PExactIKSolver::FindSolutions( subchain, index_to_use, &endPriorG, &endG, &endNextG, iks);
					if( n_proposal > 0) {
						assert( iks.size() == 1);
						subchain->MultiRotate_noGridUpdate(iks[0]);

						/* optional: free end comformation sampling
						 * We perturb the two ends therefore, we get variations on the two ends.
						 */
						if( this->freeEnd)
							if( j == 0) {
								//NOTE: Method 1
								double range = 60;
								double angle = (rand() % 100) / 100.0 * range - range / 2;
								subchain->RotateChain_noGridUpdate( "backbone", 2, backward, angle);
							}
							else if( j == num_subchains - 1) {
								//NOTE: Method 1
								double range = 60;
								double angle = (rand() % 100) / 100.0 * range - range / 2;
								subchain->RotateChain_noGridUpdate( "backbone", 5, forward, angle);
							}
						IK_success = true;
						break;
					}
					else {
						//No IK solution
						subchain->restoreChainState_noGridUpdate( state_subchain);
						num_IK_fail += 1;
					}
				}

				if( IK_success == false) {
//					cout << " No change to subchain # " << j << " due to IK failure." << endl;;
					break;
				}
				/*
				 * optional: add the side chain handling.
				 */
				if( this->use_Rotamer) {
					vector<DihedralAngle> bbangles; subchain->getDihedralAngles( bbangles);
					this->scRotater->rotateSidechain( subchain, bbangles);
				}

				/*
				 * calculate the importance ratio for the proposal block.
				 */
				double P_proposal = this->getP_log( subchain);
				int status = 0;
				double Q_proposal = this->getQ_log( subchain, n_proposal, status);
				bool accept = false;
				if( status != -1)
					accept = this->MHStep( P, Q, P_proposal, Q_proposal);
				else {
					cout << "Reject since matrix is singular; Press any key to continue." << endl;
				}

				if( accept == true) {
					//Must update before calling collision checking!
					subchain->updateAtomsGrid();

					bool collision = false;
					if( this->use_colChecking) {
						collision = subchain->InAnyCollision();
					}

					if( collision == false) {
						changed = true;
						cout << "Succeed: get an accepted conformation" << endl;
						break;
					}
					else {
						cout << "Failed: Collision detected" << endl;
						num_collision_reject += 1;
						subchain->restoreChainState_noGridUpdate( state_subchain);
					}
				}
				else {
					cout << "Failed: Rejected by Metropolis-Hasting step" << endl;
					num_MH_reject += 1;
					subchain->restoreChainState_noGridUpdate( state_subchain);
				}
			}
			//Now we have a closed sub-chain loop, destroy the record for last state.
			delete state_subchain;
		}

		/*
		 * Delete the previous stored protein chain state.
		 * Move on to sample the next conformation based on the current one
		 */
		if( changed == true) {
			stat_distinct += 1;
		}
		stat_conformation += 1;
		if( this->logFile) {
			if( (i + 1) % this->skipLength == 0) {
				int index = (int)(i / this->skipLength);
				stringstream ss;
				ss << index;
				PDBIO::writeToFile( this->protein, "../pdbfiles_out/slikmc_" + ss.str() + ".pdb");
				cout << "***Record # " << index << endl;
			}
		}

		cout << "done." << endl;
		i = i + 1;
		clock_t curr = clock();
		if(( curr - begin) / 1000.0 > time) {
			cout << "Times up!" << endl;
			break;
		}
	}

	if( this->logFile) {
		string filename = "../pdbfiles_out/info.txt";
		ofstream out( filename.c_str());
		out << "stat_distinct:\t" << stat_distinct << endl;
		out << "stat_conformation:\t" << stat_conformation << endl;
		out.flush();
		out.close();
	}
}

bool SLIKMCSampler::MHStep(double P, double Q, double P_proposal, double Q_proposal) {
	//Be careful that they are the logged probability.
	double ratio_log = P_proposal + Q - P -  Q_proposal;
	if( ratio_log >= 0)
		return true;
	else {
		double ratio = exp(ratio_log);
		double temp = (double)rand() / (double)RAND_MAX;
		if( temp < ratio)
			return true;
		else
			return false;
	}
}

double SLIKMCSampler::getP_log(PProtein* chain) {
	double log_prob_rplot = 0;
	if( this->use_RPlot) {
		int start = 0;
		int end = chain->size() - 1;
		for( int i = start; i <= end; i++) {
			double temp = this->rplot->getResidueAngleProbability( chain, i);
			log_prob_rplot += log(temp);
		}
	}

	double log_prob_bfactor = 0;
	if( this->use_BFactor)
		log_prob_bfactor = this->bfactor->evalAtomPositions( chain, chain->getTopLevelIndices().first, chain->getTopLevelIndices().second);

	double log_prob_sidechain = 0;
	if( this->use_Rotamer) {
		log_prob_sidechain = this->scRotater->evalSidechain( chain);
	}

	double log_custom = 0;
	if( this->use_customPrior) {
		for( int i = 0; i < this->priors.size(); i++) {
			log_custom += this->priors[i]->evaluate( chain->getTopLevelChain());
		}
	}
	return log_prob_rplot + log_prob_bfactor + log_prob_sidechain + log_custom;
}

double SLIKMCSampler::getQ_log(PProtein* chain, const int num_solutions, int& status) {

	//NOTE: Currently, just return the probability for the first pair of dihedral angles.
	double log_prob_firstres = 0;
	if( this->use_RPlot)
		log_prob_firstres = log(this->rplot->getResidueAngleProbability(chain, 0));
	//NOTE: The metric tensor part! The result is log.
	double log_Q = log_prob_firstres - 0.5 * this->getMetricTensor_log( chain, status) - log( num_solutions);
	return log_Q;
}

void SLIKMCSampler::enableBfactors( PProtein* p_bfactor) {
	this->use_BFactor = true;
	if (this->bfactor != NULL)
		delete this->bfactor;
	if( p_bfactor == NULL)
		this->bfactor = new BFactor( this->protein);
	else {
		Debug::check(p_bfactor->size() == this->protein->size(), "The chain providing BFactors has different length of the top level chain in sampling!");
		this->bfactor = new BFactor( p_bfactor);
	}
}

void SLIKMCSampler::disableBfactors() {
	this->use_BFactor = false;
}

void SLIKMCSampler::enableSidechain() {
	this->use_Rotamer = true;
	if( this->init_Rotamer == false) {
		assert( this->scRotater == NULL);
		this->scRotater = new SidechainRotater( this->protein);
		this->init_Rotamer = true;
	}
	this->protein->EnableSidechains();
	return;
}

void SLIKMCSampler::disableSidechain() {
	this->use_Rotamer = false;
	this->protein->DisableSidechains();
}

void SLIKMCSampler::enableFreeEnd() {
	this->freeEnd = true;
}

void SLIKMCSampler::disableFreeEnd() {
	this->freeEnd = false;
}

void SLIKMCSampler::enableLog(const int skipLength) {
	this->logFile = true;
	this->skipLength = skipLength + 1;
	assert( skipLength >= 1);
}

void SLIKMCSampler::disableLog() {
	this->logFile = false;
}

void SLIKMCSampler::enableCollisionChecking() {
	this->use_colChecking = true;
}

void SLIKMCSampler::disableCollisionChecking() {
	this->use_colChecking = false;
}

void SLIKMCSampler::enableRamachandran() {
	this->use_RPlot = true;
}

void SLIKMCSampler::disableRamachandran() {
	this->use_RPlot = false;
}

void SLIKMCSampler::dispSettings() {
	string bfactor = this->use_BFactor == true ? "enabled" : "disabled";
	string rplot = this->use_RPlot == true ? "enabled" : "disabled";
	string collision = this->use_colChecking == true ? "enabled" : "disabled";
	string sidechain = this->use_Rotamer == true ? "enabled" : "disabled";
	string freeEnd = this->freeEnd == true ? "enabled" : "disabled";
	string custom = this->use_customPrior == true ? "enabled" : "disabled";

	cout << "Sampling settings:" << endl;
	cout << "  R-Plot:\t" << rplot << endl;
	cout << "  B factors:\t" << bfactor << endl;
	cout << "  Col-Checking:\t" << collision << endl;
	cout << "  Sidechain:\t" << sidechain << endl;
	cout << "  Free end: \t" << freeEnd << endl;
	cout << "  Custom priors: \t" << custom << endl;
}

void SLIKMCSampler::enableCustomPriors() {
	this->use_customPrior = true;
}

void SLIKMCSampler::disableCustomPriors() {
	this->use_customPrior = false;
}

void SLIKMCSampler::addCustomPrior(Prior& prior) {
	this->priors.push_back( &prior);
	return;
}

double SLIKMCSampler::getMetricTensor_log( PProtein* protein, int& status) {
	//NOTE: Their code starts at 1.
	//Therefore, the matrix should be 6 + 1 (x, y, z, rx, ry, rz)by size * 2 + 1with the first row and first column junk values.
	int size = protein->size();
	double** Jac_ca =  Utility::new_Double2D( 6 + 1, size * 2 + 1);
	double** Jac_c = Utility::new_Double2D( 6 + 1, size * 2 + 1);

	//NOTE: By default, ComputeJacobian function calculates the Jacobian for the last backbone atom.
	PTools::ComputeJacobian( protein, Jac_c);
	Vector3 ca_pos = protein->getAtomPos(PID::BACKBONE, size * 3 - 2);
	PTools::ComputeJacobian( protein, Jac_ca, ca_pos);

	assert( (Jac_ca[1][8] < EPSILON) && (Jac_ca[2][8] < EPSILON) && (Jac_ca[3][8] < EPSILON));

	//Next, throw out all garbage lines and all the angular Jacobian entries
	double** dc_dx = Utility::new_Double2D( 6 + 1, 2 + 1);
	double** dc_dy = Utility::new_Double2D( 6 + 1, 6 + 1);

	for( int i = 1; i <= 3; i++)
	{
		for( int j = 1; j <= 2; j++)
		{
			dc_dx[i][j] = Jac_ca[i][j];
		}

		for( int j = 1; j <= 6; j++)
		{
			dc_dy[i][j] = Jac_ca[i][j + 2];
		}
	}

	for( int i = 4; i <= 6; i++)
	{
		for( int j = 1; j <= 2; j++)
		{
			dc_dx[i][j] = Jac_c[i - 3][j];
		}

		for( int j = 1; j <= 6; j++)
		{
			dc_dy[i][j] = Jac_c[i - 3][j + 2];
		}
	}

	double** dc_dy_inverse = Utility::new_Double2D( 6 + 1, 6 + 1);
	PNumRoutines::nr_inverse( dc_dy, 6, dc_dy_inverse, status);
	if( status == -1)
		return 1;

	double** df_dx = Utility::new_Double2D( 6 + 1, 2 + 1);
	Utility::matrix_multiply_trash( 6, 6, 2, dc_dy_inverse, dc_dx, df_dx, -1);

	double** dh_dx = Utility::new_Double2D( 2 + 1, 2 + 1);
	Utility::matrix_square_transposeA( 2, 2, df_dx, dh_dx);

	//Add the identity
	dh_dx[1][1] += 1;
	dh_dx[2][2] += 1;

	//Get the determinant
	double det = dh_dx[1][1] * dh_dx[2][2] - dh_dx[1][2] * dh_dx[2][1];

	double tensor = log(det);

	Utility::delete_Double2D( Jac_c, 6 + 1);
	Utility::delete_Double2D( Jac_ca, 6 + 1);
	Utility::delete_Double2D( dc_dx, 6 + 1);
	Utility::delete_Double2D( dc_dy, 6 + 1);
	Utility::delete_Double2D( dc_dy_inverse, 6 + 1);
	Utility::delete_Double2D( df_dx, 6 + 1);
	Utility::delete_Double2D( dh_dx, 2 + 1);
	return tensor;
}


