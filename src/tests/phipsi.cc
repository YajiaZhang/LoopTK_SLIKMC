#include "PExtension.h"
#include "PLibraries.h"
#include "PPhiPsiDistribution.h"
#include "PChain.h"
#include "PBasic.h"

#include "test.h"
#include <assert.h>

#include <map>
using std::map;

static const Real tolerance = 1.0e-6;

bool closeToOne(Real probability) {
  return fabs(1 - probability) <= tolerance;
}

int main() {
  map<string, PPhiPsiDistribution> ramachandran =
      PPhiPsiDistribution::generateRamachandran();
  Real sum;

  assert(ramachandran.size() == kNumStdAminoAcids);

  // Put the 20 standard amino acids defined in "test.h" into an STL set.
  // This allows us to use the "find" method.
  set<string> all_amino_acids(aminoAcids, aminoAcids + kNumStdAminoAcids);
  for (map<string, PPhiPsiDistribution>::const_iterator it = ramachandran.begin();
          it != ramachandran.end(); ++it) {
    assert(all_amino_acids.find(it->first) != all_amino_acids.end());

    PPhiPsiDistribution cur_distribution = it->second;
    assert(!it->second.isEmpty());

    // Verify that the unconditional phi distribution probabilities sum to 1.
    sum = 0;
    vector<double> phi_distribution = cur_distribution.getPhiDistribution();
    assert(!phi_distribution.empty());
    assert(phi_distribution.size() == cur_distribution.numPhiIntervals());
    for (unsigned i = 0; i < phi_distribution.size(); i++) {
      sum += phi_distribution[i];
    }
    assert(closeToOne(sum));

    // Verify that the unconditional psi distribution probabilities sum to 1.
    sum = 0;
    vector<double> psi_distribution = cur_distribution.getPsiDistribution();
    assert(!psi_distribution.empty());
    assert(psi_distribution.size() == cur_distribution.numPsiIntervals());
    for (unsigned i = 0; i < psi_distribution.size(); i++) {
      sum += psi_distribution[i];
    }
  }

  return 0;
}
