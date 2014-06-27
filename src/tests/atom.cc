#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

bool Implies(bool p, bool q) {
  return (!p && !q) || (p && q);
}

void BasicTest(int proteinSize) {
  PProtein *protein = new PProtein(proteinSize);

  PAtom *a1 = protein->getAtomAtRes("C", 0);
  PAtom *a2 = protein->getAtomAtRes("CA", 0);

  // Test getChain().
  assert(a1->getChain() == protein);
  assert(a2->getChain() == protein);

  // Test WithinActiveBlock.
  assert(a1->WithinActiveBlock());
  assert(a2->WithinActiveBlock());

  // Test isBonded.
  assert(a1->isBonded(a2));
  assert(a2->isBonded(a1));

  // Test isOnBackbone.
  assert(a1->isOnBackbone());
  assert(a2->isOnBackbone());

  // Test various collision methods.
  assert(Implies(a1->InStaticCollision(),
         a1->getAllCollidingStatic()->size() > 0));
  assert(Implies(a1->InSelfCollision(),
         a1->getAllCollidingSelf()->size() > 0));
  assert(Implies(a1->InAnyCollision(),
         a1->getAllCollidingEither()->size() > 0));

  // Test shortestBondPath.
  assert(PAtom::shortestBondPath(a1, a2, 1) == 1);

  protein->Obliterate();
}

int main() {
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  for(int i = 1; i <= 50; i++) {
    BasicTest(i);
  }

  return 0;
}
