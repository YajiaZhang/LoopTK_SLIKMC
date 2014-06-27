#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

void PDBCollisionTest(PChain *protein)
{
  if(protein) {
    auto_ptr<AtomCollisions> allColl(protein->getAllCollidingEither());

    if(allColl->size() > 0) {
      cerr << "Protein contains the following collisions:" << endl;

      for(AtomCollisions::const_iterator it = allColl->begin(); it != allColl->end(); ++it) {
        pair<PAtom *, PAtom *> curPair = *it;
        cerr << "(" << curPair.first->getPos() << ") and (" << curPair.second->getPos() << ")" << endl;
      }

      PUtilities::AbortProgram("Error: Protein contains collisions.");
    }
  }
}

int main() {
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  string fileName = "pdbfiles/2CRO.pdb";
  PProtein *protein = PDBIO::readFromFile(fileName);

  PDBCollisionTest(protein);  /* Run the collision test. */

  delete protein;

  return 0;
}
