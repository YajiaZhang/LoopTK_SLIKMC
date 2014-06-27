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

#include "PDataCollector.h"
#include "PChainNavigator.h"
#include "PConfSpaceNavigator.h"

#define DISTANCE_FOR_EXACT_IK 4

PDataCollector::PDataCollector(const string &cs2FileName, const string &pdbFileName,int loopStartResIndex, int loopEndResIndex) {
  PProtein *protein = PDBIO::readFromFile(pdbFileName);
  PProtein *loop = new PProtein(protein,loopStartResIndex,loopEndResIndex);
  CS2IO::writeToFile(cs2FileName,loop);
  delete protein; //destructor will destroy all children too (loop)
  InitMembers(cs2FileName);
}


PDataCollector::PDataCollector(const string &cs2FileName, PProtein *loop) {
  CS2IO::writeToFile(cs2FileName,loop);
  InitMembers(cs2FileName);
}


PDataCollector::PDataCollector(const string &cs2FileName) {
  InitMembers(cs2FileName);
}

void PDataCollector::InitMembers(const string &cs2FileName) {
  m_cs2FileName = cs2FileName;
}


void PDataCollector::Clone(const string &newActivecs2Filename) {
  PUtilities::copyFile(newActivecs2Filename,m_cs2FileName);
  m_cs2FileName = newActivecs2Filename;
}

static IKSolutions FindConformation(PProtein *loop, bool shouldBeCollisionFree) {
  PProteinCCDSolver ccd_solver(loop);
  while(ccd_solver.DoDescent()>DISTANCE_FOR_EXACT_IK) { }
  int Res_indices[3];
  Res_indices[0] = 0;
  Res_indices[1] = loop->size()/2;
  Res_indices[2] = loop->size()-1;
  IKSolutions sols = PExactIKSolver::FindSolutions(loop,Res_indices);
  IKSolutions ret;
  if (shouldBeCollisionFree) {
    for(int i=0;i<sols.size();i++) {
      loop->MultiRotate(sols[i]);
      if (!loop->InAnyCollision()) ret.push_back(sols[i]);
      loop->AntiMultiRotate(sols[i]);
    }
  } else ret = sols;
  return ret;
}

static IKSolutions FindConformationIgnoringCollision(PProtein *loop) {
  return FindConformation(loop,false);
}

static IKSolutions FindCollisionFreeConformation(PProtein *loop) {
  return FindConformation(loop,true);
}


void PDataCollector::FindConformations(int amt, bool shouldBeCollisionFree) {
  LoopConfFunction confFinder;
  if (shouldBeCollisionFree) confFinder = &FindCollisionFreeConformation;
  else
    confFinder = &FindConformationIgnoringCollision;
  FindConformations(amt,confFinder);
}

void PDataCollector::FindConformations(int amt, LoopConfFunction FindSomeConformations) {
  for(int i=0;i<amt;) {
    PProtein *loop = CS2IO::readConformation(m_cs2FileName,1);
    cerr<<"read loop"<<endl;
    loop->RandomizeDOFs(forward);
    cerr<<"randomized"<<endl;
    IKSolutions sols = FindSomeConformations(loop);
    cerr<<"got conformations"<<endl;
    for(int j=0;j<sols.size();j++) {
      loop->MultiRotate(sols[j]);
      CS2IO::appendConformation(m_cs2FileName,loop);
      loop->AntiMultiRotate(sols[j]);
    }
    cerr<<"appended conformations"<<endl;
    cerr<<sols.size()<<endl<<endl;

    i+=sols.size();
    loop->Obliterate();
  }
}

void PDataCollector::VisualizeLoop() {
  PProtein *loop = CS2IO::readConformation(m_cs2FileName,1);
  PChainNavigator(loop->getTopLevelProtein(),loop).Run();
}

void PDataCollector::VisualizeConformationSpace() {
  PConformationSpace *space = CS2IO::readConformationSpace(m_cs2FileName);
  PConfSpaceNavigator(space).Run();
}
