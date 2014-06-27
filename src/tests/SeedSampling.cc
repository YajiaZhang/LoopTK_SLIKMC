#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "fstream"
#include "sstream"
#include "string"
#include <vector>

#include "PBasic.h"
#include "PExtension.h"
#include "PResources.h"
#include "PChain.h"
#include "PTools.h"
#include "POptimize.h"
#include "PSampMethods.h"

using namespace std;

string itostring (int i) {
        stringstream ss;
        string s;
        ss << i;
        ss >> s;
        return s;
}

                                
int main (int argc, char *argv[]) {     
        WarningIndicator wi = SUPPRESS_WARNINGS;
        if (argc>1) {           
                string nowarnings = argv[1];
                if (nowarnings == "-woff") {
                        wi = SUPPRESS_WARNINGS;
                }
        }
        LoopTK::Initialize(wi);

	if (argc!=2) {
		cout << "Need one parameter: full path of SCWRL3" << endl;
		exit(1);
	}
	
	string scwrl3_path = argv[1];
	cout << "This procedure is to generate 2 collision-free loop conformations of a loop in 1MPP." << endl << endl;
	cout << "This procedure should not print out any ``backbone in collision'' message." << endl;
	cout << "However, since SCWRL3 doesn't not guarantee collision-free side chain placement," << endl;
	cout << " it is ok if ``side chain in collision'' is printed." << endl << endl;
	cout << "This procedure will write out 2 complete protein with collision-free loop backbone conformations," << endl;
	cout << " named sample0.pdb and sample1.pdb, respectively." << endl;
	cout << "If collision-free side chains are added successfully, files with _withSC.pdb as part of the name" << endl;
	cout << " will also be generated." << endl << endl;
        
	srand(time(NULL));

        map< string, PPhiPsiDistribution >  distri = PPhiPsiDistribution::generateRamachandran();
                        
        PProtein *protein = PDBIO::readFromFile("pdbfiles/1MPP.pdb");
        int loopSid = 215, loopEid = 223;
	int num_seeds_wanted = 2;
        vector<PProtein*> samples = PSampMethods::SeedSampleBackbone(protein,loopSid,loopEid,distri,num_seeds_wanted);

        for (int i=0; i<samples.size(); ++i) {

		// Check whether the sampled loop is really collision-free
                PProtein *loop = new PProtein(samples[i],loopSid,loopEid);
                loop->DisableSidechains();
                if (loop->InAnyCollision()) {
                        cout << "backbone in collision" << endl;
                }

		// Write the entire protein with sampled loop conformation out
		string file = "sample"+itostring(i);
		PDBIO::writeToFile(samples[i],file+".pdb");

		// Add on side chains and check if the side chains are in any collision
		// Since 
                file = "seed"+itostring(i);
                PDBIO::writeToFile(loop,file+".pdb");
                PSampMethods::addSidechain(file+".pdb","pdbfiles/1MPP_noLoop.pdb",scwrl3_path,file+"_withSC.pdb");
		PProtein *loop_withSC = PDBIO::readFromFile(file+"_withSC.pdb");
		PProtein *sample_withSC = PSampMethods::MergeProtein(loop_withSC,protein,loopSid);
		loop_withSC->Obliterate();
		loop_withSC = new PProtein(sample_withSC,loopSid,loopEid);
		if (loop_withSC->InAnyCollision()) {
			cout << "side chain in collision" << endl;
		}
		else {
			PDBIO::writeToFile(sample_withSC,file+"_withSC.pdb");
		}
	}

	return 0;
}

