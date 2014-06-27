/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "PBasic.h"
#include "PConstants.h"
#include "PExtension.h"
#include "PLibraries.h"
#include "PResources.h"

PProtein::PProtein(): PChain() {}

PProtein::PProtein(const string &firstResidueName): PChain(firstResidueName) {}

PProtein::PProtein(const string &firstResidueName, PResidueSpec &firstResidueSpec): PChain(firstResidueName, firstResidueSpec) {}

PProtein::PProtein(PProtein *protein, int resStartIndex, int resEndIndex): PChain(protein, resStartIndex, resEndIndex) {}

PProtein::PProtein(int numResidues)
{
  string curResidue;
  vector<string> resNames = PResources::getResidueNames();

  for(int i = 0; i < numResidues; i++) {
    curResidue = resNames[rand() % resNames.size()];
    AddResidue(curResidue);
  }

  finalize();
}

PProtein *PProtein::getTopLevelProtein() {
  return (PProtein *)getTopLevelChain();
}

PProteinResidue* PProtein::getResidueByPdbIndex(int pdbIndex) {
  int local_index = pdbIndexToLocalIndex(pdbIndex);
  return (local_index >= 0 ? getResidue(local_index) : NULL);
}

PProtein* PProtein::getSubchainByPdbIndices(int startPdbIndex, int endPdbIndex) {
  return new PProtein(this, pdbIndexToLocalIndex(startPdbIndex), pdbIndexToLocalIndex(endPdbIndex));
}

int PProtein::NumBackboneDOFs() const {
  return NumDOF(PID::BACKBONE);
}

int PProtein::NumSidechainDOFs() const {
  return NumDOF(PID::SIDECHAIN);
}

void PProtein::GetAnchors(ProteinSide side, PAtom *&endAnchor, PAtom *&endPriorAnchor)
{
  if (side == start) {
    PResidue *res = getResidue(0);

    res = res->PreviousResidue();
    if (res == NULL) PUtilities::AbortProgram("Error: residue before head is null!");

    endPriorAnchor = res->getAtom(PID::C);
    endAnchor = res->getAtom(PID::C_ALPHA);
  } else if (side == end) {
    PResidue *res = getResidue(size() - 1);

    res = res->NextResidue();
    if (res == NULL) PUtilities::AbortProgram("Error: residue after tail is null!");

    endPriorAnchor = res->getAtom(PID::N);
    endAnchor = res->getAtom(PID::C_ALPHA);
  }
}

void PProtein::GetEndEffectors(ProteinSide side, PAtom *&endAtom, PAtom *&endPriorAtom)
{
  if (side == start) {
    PResidue *res = getResidue(0);

    endPriorAtom = res->getAtom(PID::C_ALPHA);
    endAtom = res->getAtom(PID::N);
  } else if (side == end) {
    PResidue *res = getResidue(size() - 1);

    endPriorAtom = res->getAtom(PID::C_ALPHA);
    endAtom = res->getAtom(PID::C);
  }
}

void PProtein::RotateBackbone(int DOF_index, BondDirection dir, Real numDegrees)
{
  RotateChain(PID::BACKBONE, DOF_index, dir, numDegrees);
}

void PProtein::RotateSidechain(int DOF_index,float numDegrees) {
  RotateChain(PID::SIDECHAIN,DOF_index,forward,numDegrees);
}

void PProtein::DisableSidechains() {
  DisableSidechains(0,size()-1);
}

void PProtein::DisableSidechains(int resIndex1, int resIndex2) {
  DetachBlocks(PID::SIDECHAIN,PID::BACKBONE,resIndex1,resIndex2);
}

void PProtein::EnableSidechains() {
  EnableSidechains(0,size()-1);
}

void PProtein::EnableSidechains(int resIndex1,int resIndex2) {
  ReattachBlocks(PID::SIDECHAIN,resIndex1,resIndex2);
}

PProtein *PProtein::Clone() {
  PProtein *ret = new PProtein();
  CloneResiduesIntoChain(ret);
  if (getParent()!=NULL) {
    pair<int,int> indices = getTopLevelIndices();
    ret = new PProtein(ret,indices.first,indices.second);
  }
  return ret;
}


void PProtein::RandomizeSidechainAtRes(int resIndex) {
  getResidue(resIndex)->RandomizeDOFs(PID::SIDECHAIN);
}

void PProtein::RandomizeAllSidechains() {
  for(int i=0;i<size();i++) {
    RandomizeSidechainAtRes(i);
  }
}

void PProtein::RandomizeBackbone() {
  for (int i = 0; i < size(); ++i) {
    getResidue(i)->RandomizeDOFs(PID::BACKBONE);
  }
}

int PProtein::pdbIndexToLocalIndex(int pdb_index) {
  for (int i = 0; i < size(); i++) {
    PProteinResidue *cur = getResidue(i);
    if (cur != NULL && cur->getPdbId() == pdb_index) {
      return i;
    }
  }

  return -1;
}

void PProtein::printCollision() {
	for (int i=0; i<size(); ++i) {
		PResidue *res = getResidue(i);
		res->printCollision();
	}
}
