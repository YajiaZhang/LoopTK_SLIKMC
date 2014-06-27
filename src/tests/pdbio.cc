#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

int main() {
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  string fileName = "pdbfiles/1K96.pdb";
  PProtein *protein = PDBIO::readFromFile(fileName);
  PAtom *a;

  for (unsigned i = 0; i < protein->size(); i++) {
    a = protein->getAtomAtRes("N", i);
    assert(a->getOccupancy() == 0.5);

    a = protein->getAtomAtRes("C", i);
    assert(a->getOccupancy() == 0.5);

    a = protein->getAtomAtRes("CA", i);
    assert(a->getOccupancy() == 0.5);

    assert(protein->getResidue(i)->getPdbId() == i + 1);
  }

  for (int i = 1; i <= protein->size(); i++) {
    PProteinResidue *res = protein->getResidueByPdbIndex(i);

    assert(res != NULL);
    assert(res->getPdbId() == i);
  }

  PProtein *sub_protein = protein->getSubchainByPdbIndices(1, protein->size());
  assert (sub_protein->size() == protein->size());
  for (unsigned i = 0; i < sub_protein->size(); i++) {
    a = sub_protein->getAtomAtRes("N", i);
    assert(a->getOccupancy() == 0.5);

    a = sub_protein->getAtomAtRes("C", i);
    assert(a->getOccupancy() == 0.5);

    a = sub_protein->getAtomAtRes("CA", i);
    assert(a->getOccupancy() == 0.5);

    assert(sub_protein->getResidue(i)->getPdbId() == i + 1);
  }

  delete sub_protein;
  delete protein;

  return 0;
}
