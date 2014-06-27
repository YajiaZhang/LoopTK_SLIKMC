#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

void BasicTest(int proteinSize) {
  PProtein *protein = new PProtein(proteinSize);

  // Some basic sanity checks on the protein.
  assert(protein->size() == proteinSize);
  assert(protein->getParent() == NULL);
  assert(protein->getTopLevelChain() == protein);

  // Make sure the backbone atoms are present.
  for (unsigned i = 0; i < protein->size(); i++) {
    assert(protein->getAtomAtRes("C", i) != NULL);
    assert(protein->getAtomAtRes("CA", i) != NULL);
    assert(protein->getAtomAtRes("N", i) != NULL);
  }
  assert(protein->NumAtoms("backbone") == protein->size() * 3);

  // Collision sanity checks.
  if (protein->InAnyCollision()) {
    assert(protein->FindAnyCollision().first != NULL);
    assert(protein->getAllCollidingEither()->size() > 0);
  } else {
    assert(protein->FindAnyCollision().first == NULL);
    assert(protein->getAllCollidingEither()->size() == 0);
  }

  protein->Obliterate();
}

int main() {
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  for(int i = 1; i <= 50; i++) {
    BasicTest(i);
  }

  return 0;
}
