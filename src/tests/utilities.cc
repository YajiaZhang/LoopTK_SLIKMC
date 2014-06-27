#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

int main() {
  // Test PUtilities::trimSpaces()
  string s = "                          ";
  PUtilities::trimSpaces(s);
  assert(s.length() == 0);

  // Test PUtilities::toStr()
  assert(PUtilities::toStr(5) == "5");
  assert(PUtilities::toStr(-20) == "-20");

  // Test PUtilities::padLeft() and PUtilities::padRight()
  s = "LoopTK";
  for (int i = 1; i <= s.length(); i++) {
    assert(PUtilities::padLeft(s, i).length() == s.length());
    assert(PUtilities::padRight(s, i).length() == s.length());
  }
  for (int i = s.length() + 1; i <= 100; i++) {
    assert(PUtilities::padLeft(s, i).length() == i);
    assert(PUtilities::padRight(s, i).length() == i);
  }

  // Test PUtilities::getTokens()
  s = "";
  string testTokens[4] = {"this", "is", "a", "test"};
  for (unsigned int i = 0; i < 4; i++) {
    s += testTokens[i] + " ";
  }
  vector<string> tokens = PUtilities::getTokens(s);
  assert(tokens.size() == 4);
  for (unsigned i = 0; i < tokens.size(); i++) {
    assert(tokens[i] == testTokens[i]);
  }

  // Test PUtilities::getLinesStarting()
  tokens = PUtilities::getLinesStarting(tokens, "t", true);
  assert(tokens.size() == 2);
  assert(tokens[0] == testTokens[0].substr(1));
  assert(tokens[1] == testTokens[3].substr(1));

  // Test PUtilities::PointerThatIsNot()
  int a = 0, b = 0;
  int *c = &a, *d = &b, *e = &a;

  assert(PUtilities::PointerThatIsNot(c, d, e) == d);
  assert(PUtilities::PointerThatIsNot(d, c, e) == d);

  // Test isNAN
  assert(isNAN(0) == false);

  return 0;
}
