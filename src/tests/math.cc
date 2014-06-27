#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"

#include <assert.h>
#include <math.h>

#define  THRESHOLD    1e-5

bool closeTo(double x, double y) {
  return (fabs(x-y) <= THRESHOLD);
}

void testGetAngle() {
  assert(closeTo(PMath::getAngle(Vector3(1, 0, 0), Vector3(0, 0, 1)), PI/2));
  assert(closeTo(PMath::getAngle(Vector3(1, 0, 0), Vector3(0, 1, 0)), PI/2));
  assert(closeTo(PMath::getAngle(Vector3(0, 1, 0), Vector3(0, 0, 1)), PI/2));

  assert(closeTo(PMath::getAngle(Vector3(1, 0, 0), Vector3(1, 0, 0)), 0));
  assert(closeTo(PMath::getAngle(Vector3(0, 1, 0), Vector3(0, 1, 0)), 0));
  assert(closeTo(PMath::getAngle(Vector3(0, 0, 1), Vector3(0, 0, 1)), 0));
}

void testSphereIntersectsCube() {
  assert(PMath::sphereIntersectsCube(Vector3(0, 0, 0),
                                     1.0,
                                     Vector3(-1, -1, -1),
                                     Vector3(1, 1, 1)));

  assert(!PMath::sphereIntersectsCube(Vector3(1, 1, 1),
                                      0.5,
                                      Vector3(3, 3, 3),
                                      Vector3(5, 5, 5)));
}

void testSignum() {
  for (int i = -100; i <= 100; i++) {
    if (i < 0) {
      assert(PMath::signum(i) == -1);
    } else {
      assert(PMath::signum(i) == 1);
    }
  }
}

int main() {
  testSphereIntersectsCube();
  testGetAngle();
  testSignum();

  return 0;
}
