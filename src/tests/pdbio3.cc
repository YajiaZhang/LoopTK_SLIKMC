#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

int main (int argc, char *argv[]) {
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  PProtein *protein = PDBIO::readFromFile("pdbfiles/1MPP.pdb");

  // This PDB file has residue IDs that start off as one would
  // expect (1, 2, 3, ...) but eventually get discontinuous.
  for (int i = 1; i < protein->size(); i++) {
    if (i < 109) {
      assert(protein->getResidue(i)->getPdbId() == i);
    } else if (i < 159) {
      assert(protein->getResidue(i)->getPdbId() == i + 1);
    } else if (i < 240) {
      assert(protein->getResidue(i)->getPdbId() == i + 2);
    } else if (i < 289) {
      assert(protein->getResidue(i)->getPdbId() == i + 3);
    } else {
      assert(protein->getResidue(i)->getPdbId() == i + 8);
    }
  }

  delete protein;

  return 0;
}

