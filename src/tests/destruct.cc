#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

void DestructorTest(PChain *protein) {
  PChain *sub1 = new PChain(protein, 0,6);
  PChain *sub2 = new PChain(sub1, 0,5);     
  PChain *sub3 = new PChain(sub2, 0,4);     
  PChain *sub4 = new PChain(sub3, 0,3);

  delete sub4;
  delete sub3;
  delete sub2;
  delete sub1;
}

int main() {
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  string fileName = "pdbfiles/2CRO.pdb";
  PProtein *protein = PDBIO::readFromFile(fileName);

  DestructorTest(protein);  /* Run the multi-destructor test. */

  delete protein;

  return 0;
}
