#include "PBasic.h"
#include "PChain.h"
#include "PDataCollector.h"
#include "PExtension.h"
#include "PHydrogenBondTracker.h"
#include "PTools.h"
#include "PConfSpaceNavigator.h"
#include <iostream>
#include "Utility.h"
#include "LoopTKSampler.h"
#include "PNumRoutines.h"
#include "Rotamer.h"
#include "SLIKMC.h"
#include "MetropolisHastingsSampler.h"
#include "Prior.h"
using namespace std;

class SamplePrior:public Prior {
public:
	double evaluate(PChain* chain);
};

double SamplePrior::evaluate( PChain* chain) {
	vector<int> index;
	for( int i = 0; i < chain->size(); i++) {
		index.push_back(i);
	}
	Vector3 t1 = chain->getResidue(0)->getAtomPosition(PID::N);
	Vector3 t2 = chain->getResidue(chain->size() - 1)->getAtomPosition( PID::C);
	Vector3 pos = (t1 + t2) / 2;

	double log_prob = 0;
	for (int i = 0; i < index.size(); i++) {
		PResidue* residue = chain->getResidue( index[i]);
		Vector3 Ca = residue->getAtomPosition(PID::C_ALPHA);
		double dCa = pos.distance(Ca);
		log_prob += 2 * dCa;
	}
	return log_prob;
}

int main(int argc, char *argv[]) {
	WarningIndicator wi = SUPPRESS_WARNINGS;

	if (argc > 1) {
		string nowarnings = argv[1];
		if (nowarnings == "-woff") {
			wi = SUPPRESS_WARNINGS;
		}
	}
	LoopTK::Initialize(wi);

	PProtein* chain = PDBIO::readFromFile("../pdbfiles/loop0.pdb");
	SLIKMCSampler sampler( chain);
	sampler.enableCollisionChecking();
//	sampler.enableBfactors();
	sampler.disableBfactors();
	sampler.enableRamachandran();
	sampler.enableLog( 100);
//	sampler.enableSidechain();
	sampler.dispSettings();

	sampler.enableCustomPriors();
	SamplePrior prior;
	sampler.addCustomPrior( prior);
	sampler.sample( 100, 0, chain->size() - 1);

//	LoopTKSampler sampler( chain);
//	sampler.enableBFactors();
//	sampler.sample( 100, 0, chain->size() - 1, 10);
}

