#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

int main (int argc, char *argv[]) {
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  PProtein *protein = PDBIO::readFromFile("pdbfiles/1K8U.pdb");

  for (unsigned i = 0; i < protein->size(); i++) {
    assert(protein->getResidue(i)->getPdbId() == i + 2);
  }

  for (int i = 2; i <= protein->size() + 1; i++) {
    PProteinResidue *res = protein->getResidueByPdbIndex(i);

    assert(res != NULL);
    assert(res->getPdbId() == i);
  }

  PProtein *sub_protein = protein->getSubchainByPdbIndices(2, protein->size() + 1);
  assert (sub_protein->size() == protein->size());
  for (unsigned i = 0; i < sub_protein->size(); i++) {
    assert(sub_protein->getResidue(i)->getPdbId() == i + 2);
  }

  delete sub_protein;
  delete protein;

  return 0;
}
