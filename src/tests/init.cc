#include "PExtension.h"
#include "PLibraries.h"
#include "PChainNavigator.h"
#include "PChain.h"
#include "PBasic.h"
#include <assert.h>

#include "test.h"

int main() {
  LoopTK::Initialize(SUPPRESS_WARNINGS);

  // Make sure we have at least the standard atoms and amino acids.
  assert(PResources::numAtoms() >= kNumStdAtoms);
  assert(PResources::numBlocks() >= kNumStdAminoAcids);
  assert(PResources::numConnections() >= kNumStdAminoAcids);
  assert(PResources::numResidues() >= kNumStdAminoAcids);

  for (int i = 0; i < kNumStdAminoAcids; i++) {
    string cur_amino_acid = aminoAcids[i];

    // Test that the PResources::Contains* methods work.
    assert(PResources::ContainsBlockShell(cur_amino_acid));
    assert(PResources::ContainsBlockConnection(PID::BACKBONE, cur_amino_acid));
    assert(PResources::ContainsResidueShell(cur_amino_acid));

    // Test PResources::GetBlockShell.
    assert(PResources::GetBlockShell(cur_amino_acid) != NULL);
    assert(PResources::GetBlockShell(cur_amino_acid)->getName() == cur_amino_acid);
    assert(PResources::GetBlockShell(cur_amino_acid)->getType() == PID::SIDECHAIN);

    // Test PResources::GetBlockConnection.
    assert(PResources::GetBlockConnection(PID::BACKBONE, cur_amino_acid) != NULL);

    // Test PResources::GetResidueShell.
    assert(PResources::GetResidueShell(cur_amino_acid) != NULL);
    assert(PResources::GetResidueShell(cur_amino_acid)->getName() == cur_amino_acid);

    // Test that all amino acids have some chi angles (with the exceptions
    // of ALA and GLY), and that each chi index specifies 4 atoms.
    if (cur_amino_acid == PID::ALA || cur_amino_acid == PID::GLY) {
      assert(PResources::numChiIndices(cur_amino_acid) == 0);
      assert(!PResources::ContainsRotamer(cur_amino_acid));
    } else {
      assert(PResources::numChiIndices(cur_amino_acid) > 0);
      assert(PResources::ContainsRotamer(cur_amino_acid));

      vector<vector<Real> > rotamers = PResources::GetRotamer(cur_amino_acid);
      for (unsigned i = 0; i < rotamers.size(); i++) {
        assert(rotamers[i].size() == PResources::numChiIndices(cur_amino_acid));
      }

      for (int j = 0; j < PResources::numChiIndices(cur_amino_acid); j++) {
        assert(PResources::GetChiIndex(cur_amino_acid, j + 1).size() == 4);
      }

      vector<string> cur_chi_indices = PResources::GetChiIndex(cur_amino_acid, 1);
      assert(cur_chi_indices[0] == PID::N);
      assert(cur_chi_indices[1] == PID::C_ALPHA);
      assert(cur_chi_indices[2] == PID::C_BETA);
    }
  }

  // Check that we have atom shells for all the standard atoms.
  for (int i = 0; i < kNumStdAtoms; i++) {
    assert(PResources::ContainsAtomShell(atoms[i]));
    assert(PResources::GetAtomShell(atoms[i]) != NULL);

    for (int j = 0; j < kNumStdAtoms; j++) {
      assert(PResources::GetEpsilonValue(make_pair(atoms[i], atoms[j])) > 0);
      assert(PResources::GetEpsilonValue(make_pair(atoms[i], atoms[j])) ==
             PResources::GetEpsilonValue(make_pair(atoms[j], atoms[i])));
    }
  }

  // Check that FreeResources really does clear out everything
  // in PResources.
  PResources::FreeResources();

  assert(PResources::numAtoms() == 0);
  assert(PResources::numBlocks() == 0);
  assert(PResources::numConnections() == 0);
  assert(PResources::numResidues() == 0);

  for (int i = 0; i < kNumStdAminoAcids; i++) {
    assert(PResources::ContainsBlockShell(aminoAcids[i]) == false);
    assert(PResources::ContainsBlockConnection(PID::BACKBONE, aminoAcids[i]) == false);
    assert(PResources::ContainsResidueShell(aminoAcids[i]) == false);
  }

  for (int i = 0; i < kNumStdAtoms; i++) {
    assert(PResources::ContainsAtomShell(atoms[i]) == false);
  }

  return 0;
}
