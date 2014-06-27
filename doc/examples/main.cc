#include "PBasic.h"
#include "PChain.h"
#include "PChainNavigator.h"
#include "PDataCollector.h"
#include "PExtension.h"
#include "PHydrogenBondTracker.h"
#include "PTools.h"
#include "PConfSpaceNavigator.h"
#include "PIKAlgorithms.h"
#include "RamachandranPlot.h"
#include "PChainState.h"
#include <iostream>
using namespace std;


int main(int argc, char *argv[]) {
	WarningIndicator wi = SUPPRESS_WARNINGS;

	if (argc > 1) {
		string nowarnings = argv[1];
		if (nowarnings == "-woff") {
			wi = SUPPRESS_WARNINGS;
		}
	}
	LoopTK::Initialize(wi);

	PProtein *protein = PDBIO::readFromFile("loop.pdb");
	PChainNavigator(protein).Run();
	cout << "Before Navigator run!" << endl;
	return 0;
}

