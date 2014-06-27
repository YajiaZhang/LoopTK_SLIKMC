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
#include "PResources.h"
#include "POptimize.h"
#include "PTools.h"
class FunctCalculator: public FunctFunctor {
  public:
  FunctCalculator(vector<PProtein*> loop, vector<vector<CDof> > Dofs, vector<int> AtomToMove, vector<Vector3> MovementDirn) {
    funct_value = 0;
    Tfunct_value = 0;
    lps = loop;
    DofsToUse = Dofs;
    Atom = AtomToMove;
    for(int i=0;i<lps.size();i++){
      ndofs.push_back(lps[i]->size()*2);
      Vector3 pos( MovementDirn[i].x+lps[i]->getAtomPos(PID::BACKBONE,Atom[i]).x,
      MovementDirn[i].y+lps[i]->getAtomPos(PID::BACKBONE,Atom[i]).y,MovementDirn[i].z+lps[i]->getAtomPos(PID::BACKBONE,Atom[i]).z);
      FinalPos.push_back(pos);
    }
  }

  double operator()(double p[]) {
    int k = 1;
    Tfunct_value = 0;
    for(int i=0;i<lps.size();i++){
      IKSolutions solutions;
      ChainMove CMove;
      vector<ChainMove> AllMove;
      AllMove.clear();
      for(int j=0;j<DofsToUse[i].size();j++){
        CMove.blockType = (DofsToUse[i])[j].blockType;
        CMove.dir = (DofsToUse[i])[j].dir;
        CMove.DOF_index = (DofsToUse[i])[j].DOF_index;
        CMove.degrees = p[k]*rad2deg;
        AllMove.push_back(CMove);
        k++;
      }
      solutions.push_back(AllMove);
      lps[i]->MultiRotate(solutions[0]);
      double Delta[4];
      Delta[1] = -FinalPos[i].x + lps[i]->getAtomPos(PID::BACKBONE,Atom[i]).x;
      Delta[2] = -FinalPos[i].y + lps[i]->getAtomPos(PID::BACKBONE,Atom[i]).y;
      Delta[3] = -FinalPos[i].z + lps[i]->getAtomPos(PID::BACKBONE,Atom[i]).z;
      Tfunct_value += (Delta[1]*Delta[1]+Delta[2]*Delta[2]+Delta[3]*Delta[3]);
      lps[i]->AntiMultiRotate(solutions[0]);
    }
    return Tfunct_value;
  }
  double operator()(double p[],int i) {
    IKSolutions solutions;
    ChainMove CMove;
    vector<ChainMove> AllMove;
    AllMove.clear();
    for(int j=0;j<DofsToUse[i].size();j++){
      CMove.blockType = (DofsToUse[i])[j].blockType;
      CMove.dir = (DofsToUse[i])[j].dir;
      CMove.DOF_index = (DofsToUse[i])[j].DOF_index;
      CMove.degrees = p[j+1]*rad2deg;
      AllMove.push_back(CMove);
    }
    solutions.push_back(AllMove);
    lps[i]->MultiRotate(solutions[0]);
    double Delta[4];
    Delta[1] = -FinalPos[i].x + lps[i]->getAtomPos(PID::BACKBONE,Atom[i]).x;
    Delta[2] = -FinalPos[i].y + lps[i]->getAtomPos(PID::BACKBONE,Atom[i]).y;
    Delta[3] = -FinalPos[i].z + lps[i]->getAtomPos(PID::BACKBONE,Atom[i]).z;
    funct_value = (Delta[1]*Delta[1]+Delta[2]*Delta[2]+Delta[3]*Delta[3]);
    lps[i]->AntiMultiRotate(solutions[0]);
    return funct_value;
  }
    
  private:
  double funct_value, Tfunct_value;
  vector<PProtein*> lps;
  vector<Vector3> FinalPos;
  vector<int> Atom;
  vector<int> ndofs;
  vector<vector<CDof> > DofsToUse;
};
//Atoms are numbered from
vector<IKSolutions> PMoveAtom::MoveAtom(vector<PProtein*> loops, vector<vector<CDof> > Dofs, vector<int> AtomToMove, vector<Vector3> MovementDirn)
{
  FunctCalculator *FunctToOptimize = new FunctCalculator(loops, Dofs, AtomToMove, MovementDirn);
  POptimize *opt = new POptimize(loops, Dofs);
  vector<IKSolutions> CombinedSoln = opt->OptimizeNullSpace(FunctToOptimize).LoopSol;
  for(int i=0;i<loops.size();i++){
    IKSolutions Soln = CombinedSoln[i];
    if ((Soln[1])[0].DOF_index!=-1){
      loops[i]->MultiRotate(Soln[0]);
      loops[i]->MultiRotate(Soln[1]);
        }
    else cout<<"No solution found"<<endl;
  }
  delete FunctToOptimize;
  delete opt;
  return CombinedSoln;
}

IKSolutions PMoveAtom::MoveAtom(PProtein *loop, int AtomToMove, Vector3 MovementDirn){
  vector<PProtein *> loops;
  loops.push_back(loop);
  vector<vector<CDof> > Dofs;
  Dofs = PTools::GetBBDofs(loops);
  vector<int> Atom;
  Atom.push_back(AtomToMove);
  vector<Vector3> MovementDir;
  MovementDir.push_back(MovementDirn);
  vector<IKSolutions> CombinedSoln = MoveAtom(loops, Dofs, Atom, MovementDir);
  return CombinedSoln[0];
}
