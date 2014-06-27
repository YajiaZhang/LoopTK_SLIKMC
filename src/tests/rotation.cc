#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

void RotationTest(PChain *protein)
{
  int rotation;

  for(int j = 0; j < 50; j++){
    for(int i = 0; i < protein->size(); i++){
      rotation = rand() % 360;
      protein->RotateChain("backbone", i, forward, rotation);
    }
  }
}

int main() {
  srand(unsigned(time(NULL)));
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  string fileName = "pdbfiles/2CRO.pdb";
  PProtein *protein = PDBIO::readFromFile(fileName);

  RotationTest(protein);  /* Run the rotation test. */

  delete protein;

  return 0;
}
