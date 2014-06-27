#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

int main() {
  PBlockBondShell bond(PID::C, PID::N, false);

  assert(bond.getDefinedAtomId() == PID::C);
  assert(bond.getToDefineAtomId() == PID::N);
  assert(bond.isDOF() == false);

  return 0;
}
