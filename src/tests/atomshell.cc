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
  GLColor color = GLColor(1.0, 0.5, 0.0, 1.0);

  PAtomShell shell(1.0, 2.0, "foo", color);

  // Test getCovalentRadius()
  assert(shell.getCovalentRadius() == 1.0);

  // Test getVanDerWaalsRadius()
  assert(shell.getVanDerWaalsRadius() == 2.0);

  // Test getName()
  assert(shell.getName() == "foo");
}

int main() {
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  for(int i = 1; i <= 50; i++) {
    BasicTest(i);
  }

  return 0;
}
