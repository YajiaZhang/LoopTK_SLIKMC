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

#include "POptimize.h"
#include "PNumRoutines.h"
#include "PTools.h"

POptimize::POptimize(vector<PProtein*> loops, vector<vector<CDof> > Dofs){
  lps = loops;
  DofsToUse = Dofs;
  for(int i=0;i<loops.size();i++){
    nRsds.push_back(lps[i]->size());
    nbb.push_back(nRsds[i]*3);
    ndofs.push_back(nRsds[i]*2);
    endPG.push_back(Vector3(lps[i]->getAtomPos(PID::BACKBONE,nbb[i]-2).x,
                            lps[i]->getAtomPos(PID::BACKBONE,nbb[i]-2).y,
                            lps[i]->getAtomPos(PID::BACKBONE,nbb[i]-2).z));
    endG.push_back(Vector3(lps[i]->getAtomPos(PID::BACKBONE,nbb[i]-1).x,
                           lps[i]->getAtomPos(PID::BACKBONE,nbb[i]-1).y,
                           lps[i]->getAtomPos(PID::BACKBONE,nbb[i]-1).z));
    endNG.push_back(lps[i]->getAtomAtRes(PID::O,nRsds[i]-1)->getPos());
  }
}

IKSolution POptimize::IKCloseChain(PProtein *lp, Vector3 endPriorGoal,Vector3 endGoal, Vector3 endNextGoal){
  int num_res = lp->size();
  IKSolutions final_solutions;
  int manip_pos[3][3] = {
    {num_res-3, num_res-2, num_res-1},
    {num_res-4, num_res-2, num_res-1},
    {num_res-4, num_res-3, num_res-1},
  };
  for (int posIndex=0; posIndex<3; posIndex++) {
      IKSolutions pos_solutions = PExactIKSolver::FindSolutions(lp, manip_pos[posIndex], &endPriorGoal, &endGoal, &endNextGoal);
    for (int solIndex=0; solIndex<pos_solutions.size(); solIndex++)
      final_solutions.push_back(pos_solutions[solIndex]);
  }
  // if no solutions are returned
  if (!final_solutions.size()){
    cerr<<"\nNo solutions"<<endl;
    IKSolution solution;
    ChainMove CMove;
    CMove.blockType = PID::BACKBONE;
    CMove.dir = forward;
    CMove.DOF_index = -1;
    CMove.degrees = 0;
    solution.push_back(CMove);
    return solution;
  }
  int best_index = -1;
  double lowest_rmsd = 1e+10;
  
  for (int solIndex=0; solIndex < final_solutions.size(); solIndex++) {
    PProtein* temp_lp = lp->Clone();
    temp_lp->MultiRotate(final_solutions[solIndex]);
    double this_rmsd = PTools::RMSDBackbone(lp, temp_lp, 0, lp->size()-1);
    if (this_rmsd < lowest_rmsd) {
      lowest_rmsd = this_rmsd;
      best_index=solIndex;
    }
    temp_lp->Obliterate();
  }
  if (best_index!=-1)return final_solutions[best_index];
}

class DerivCalculator: public DerivFunctor {
  public:
  DerivCalculator(vector<PProtein*> loops,  vector<vector<CDof> > dofs, FunctFunctor *FunctToOpt, vector<double> InitC) {
    lps = loops;
    DofsToUse = dofs;
    derivv.clear();
    m_FunctToOpt  = FunctToOpt;
    derivv.clear();
    IC= InitC;
  }
  ~DerivCalculator(){
  }
  
  void operator()(double p[], double d[]){
    ChainMove CMove;
    int piterate = 1;//iterator for p[]
    derivv.clear();
    int totloopDOF=0;
    for(int i=0;i<lps.size();i++) totloopDOF += DofsToUse[i].size();
    int totDOF = totloopDOF+IC.size();
    for(int i=0;i<lps.size();i++){
      IKSolutions solutions;
      int loopsize = DofsToUse[i].size();
      vector<ChainMove> AllMove;
      AllMove.clear();
      for(int j=0;j<loopsize;j++){
        CMove.dir = (DofsToUse[i])[j].dir;
        CMove.blockType = (DofsToUse[i])[j].blockType;
        CMove.DOF_index = (DofsToUse[i])[j].DOF_index;
        CMove.degrees = p[piterate++]*rad2deg;
        AllMove.push_back(CMove);
      }
      solutions.push_back(AllMove);
      lps[i]->MultiRotate(solutions[0]);
      vector<double> ybv,ysv;
      double pv[loopsize+IC.size()+1];
      double deltaang = 0.000002;//1E-7
      for(int j=1;j<=loopsize;j++) pv[j]=0.0;
      for(int j=1;j<=IC.size();j++) pv[j+loopsize]=p[totloopDOF+j];
      for(int j=1;j<=loopsize;j++){
        pv[j] = deltaang;
        double temp1 = (*m_FunctToOpt)(pv,i);
        pv[j] = -deltaang;
        double temp2 = (*m_FunctToOpt)(pv,i);
        pv[j] = 0.0;
        if ((DofsToUse[i])[j-1].blockType==PID::BACKBONE) ybv.push_back((temp1 - temp2)/(2.0*deltaang));
        else if ((DofsToUse[i])[j-1].blockType==PID::SIDECHAIN)
        ysv.push_back((temp1 - temp2)/(2.0*deltaang));
      }
      //Enumerate backbone's DOFs for Jacobian computation
      vector<int> JacInd;
      for(int j=0;j<loopsize;j++){
        if ((DofsToUse[i])[j].blockType==PID::BACKBONE){
          JacInd.push_back((DofsToUse[i])[j].DOF_index);
        }
      }
      double yb[ybv.size()+1], ys[ysv.size()+1];
      double derivb[ybv.size()+1];
      for(int j=1;j<=ybv.size();j++) yb[j]=ybv[j-1];
      for(int j=1;j<=ysv.size();j++) ys[j]=ysv[j-1];
      PTools::ProjectOnNullSpace(lps[i],JacInd, true, yb, derivb);
      double deriv[loopsize+1];
      int iterb = 1;
      int iters = 1;
      for(int j=1;j<=loopsize;j++){
        if ((DofsToUse[i])[j-1].blockType==PID::BACKBONE) deriv[j]=derivb[iterb++];
        if ((DofsToUse[i])[j-1].blockType==PID::SIDECHAIN) deriv[j]=ys[iters++];
      }
      for(int j=1;j<=loopsize;j++) derivv.push_back(deriv[j]);
      lps[i]->AntiMultiRotate(solutions[0]);
    }
    
    double pt[totDOF+1];
    double derivIC[IC.size()+1];
    piterate = 1;
    for(int i=1;i<=totloopDOF;i++) pt[i]=p[piterate++];
    for(int j=1;j<=IC.size();j++){ 
      pt[j+totloopDOF]=p[piterate++];
    }
    double delta = 0.0001;
    for(int j=1;j<=IC.size();j++){
      pt[j+totloopDOF] = pt[j+totloopDOF]+delta;
      double temp1 = (*m_FunctToOpt)(pt);
      pt[j+totloopDOF] = pt[j+totloopDOF]-2*delta;
      double temp2 = (*m_FunctToOpt)(pt);
      pt[j+totloopDOF] = pt[j+totloopDOF]+delta;
      derivIC[j] = ((temp1 - temp2)/(2.0*delta));
    }

    for(int i=0;i<derivv.size();i++) d[i+1]=derivv[i];
    for(int j=1;j<=IC.size();j++) d[j+derivv.size()] = derivIC[j]; 
    //double deriv_norm = norm(d, totDOF);
    //if (deriv_norm > 1E-4) for(int j=1;j<=totDOF;j++) d[j] = d[j]/deriv_norm;
  }
  
  // Compute the norm of a vector.
  double  norm(double * a, int len){
    double s = 0;
    for (int i = 1; i <= len; i++) s += a[i]*a[i];
    return sqrt(s);
  }
  
  private:
  Real funct_value;
  vector<PProtein*> lps;
  FunctFunctor *m_FunctToOpt;
  double** Jac;
  int n_ns;
  NullSpaceRet *Ret;
  vector<double> derivv;
  vector<double> IC;
  vector<vector<CDof> > DofsToUse;
};

OptimalSol POptimize::Optimize(FunctFunctor *FunctToOptimize, DerivFunctor *DerivOfFunct){
  vector<double> InitC;
  InitC.clear();
  return Optimize(FunctToOptimize, DerivOfFunct, InitC);
}

OptimalSol POptimize::OptimizeNullSpace(FunctFunctor *FunctToOptimize){
  vector<double> InitC;
  InitC.clear();
  return OptimizeNullSpace(FunctToOptimize, InitC);
}

OptimalSol POptimize::OptimizeNullSpace(FunctFunctor *FunctToOptimize, vector<double> InitC){
  DerivCalculator *DerivOfFunct = new DerivCalculator(lps, DofsToUse, FunctToOptimize, InitC);
  OptimalSol OptimalSolution = Optimize(FunctToOptimize,DerivOfFunct,InitC);
  delete DerivOfFunct;
  return OptimalSolution;
}

OptimalSol POptimize::Optimize(FunctFunctor *FunctToOptimize, DerivFunctor *DerivOfFunct, vector<double> InitC){
  int iter;
  double fret;
  int n_DofsToUse = 0;
  for(int i=0;i<DofsToUse.size();i++) n_DofsToUse+=DofsToUse[i].size(); 
  double sol[n_DofsToUse+InitC.size()+1];
  for(int i=1;i<=n_DofsToUse;i++) sol[i] = 0.0; 
  for(int i=1;i<=InitC.size();i++) sol[n_DofsToUse+i] = InitC[i-1];
  //cout<<"wts:"<<InitC[0]<<","<<InitC[1]<<endl;
  PNumRoutines::nr_multimin(sol, n_DofsToUse,1.0E-2, &iter, &fret, FunctToOptimize,DerivOfFunct);
  ChainMove CMove;
  int k = 1;
  vector<IKSolutions> CombinedSoln;
  for(int i=0;i<lps.size();i++){
    IKSolution solution_jac;
    for(int j=0;j<DofsToUse[i].size();j++){
      CMove.dir = (DofsToUse[i])[j].dir;
      CMove.blockType = (DofsToUse[i])[j].blockType;
      CMove.DOF_index = (DofsToUse[i])[j].DOF_index;
      CMove.degrees = sol[k++]*rad2deg;
      solution_jac.push_back(CMove);
    }
    IKSolutions solutions;
    solutions.push_back(solution_jac);
    lps[i]->MultiRotate(solutions[0]);
    IKSolution solution_ik = IKCloseChain(lps[i],endPG[i],endG[i],endNG[i]);
    lps[i]->AntiMultiRotate(solutions[0]);
    solutions.push_back(solution_ik);
    CombinedSoln.push_back(solutions);
  }
  vector<double> NLSoln;
  for(int i=1;i<=InitC.size();i++) NLSoln.push_back(sol[k++]);
  OptimalSol OptimalSolution;
  OptimalSolution.LoopSol = CombinedSoln;
  OptimalSolution.NonLoopSol = NLSoln;
  
  return OptimalSolution;

}


