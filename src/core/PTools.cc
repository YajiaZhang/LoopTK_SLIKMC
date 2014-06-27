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
#include "PResources.h"
#include "PTools.h"


PChain *PTools::LowestCommonChain(PChain *c1, PChain *c2) {
	for(PChain *c3 = c1;c3!=NULL;c3=c3->getParent()) {
		if(c2->IsSubChainOf(c3)) return c3;
	}
	return NULL;
}

PProtein *PTools::CreateSlimProtein(PProtein *protein) {
	return CreateSlimProtein(protein,0,protein->size()-1);
}

PProtein *PTools::CreateSlimProtein(PProtein *protein, int startRes, int endRes) {
	PProtein *ret = new PProtein();
	assert(startRes>=0);
	assert(endRes<protein->size());
	assert(startRes<=endRes);
	for(int i=0;i<protein->size();i++) {
		PResidue *res = protein->getResidue(i);
		PResidueSpec spec = res->getSpec();
		string type;
		if(startRes<=i&&i<=endRes) {
			type = PID::CB_RES;
			if(res->getAtom(PID::C_BETA)==NULL) type = PID::BACKBONE_RES;
		} else {
			type = res->getResourceName();
		}
		ret->AddResidue(type,spec)->setName(res->getName());
	}
	

	ret->finalize();
	return ret;
}

void PTools::GetGoals(PAtom *anchorPrior, PAtom *anchorEnd, PAtom *effectorPrior, PAtom *effectorEnd,
      Real len, Real angle1, Real angle2, Vector3 &goalPrior, Vector3 &goalEnd)
{
        Vector3 a, b, n, x;
        Matrix3 m;
        AngleAxisRotation r;

        Real goalDist = effectorPrior->getPos().distance(effectorEnd->getPos());
        goalEnd = PMath::ComputePos(anchorEnd->getPos(), anchorPrior->getPos(), angle2, len);

        // Compute normal vector to the plane of m_goalEnd and the two anchors.
        a = anchorEnd->getPos() - goalEnd;
        b = anchorPrior->getPos() - goalEnd;
        n = cross(a, b);

  m = PMath::FindRotationMatrix(n,DtoR(angle1));
  x = m*b;

        // Scale x to have the proper length and add it to goalEnd
        x.inplaceNormalize();
        x.inplaceScale(goalDist);
        goalPrior = goalEnd + x;
}

void PTools::GetGoals(PProtein *protein, ProteinSide side, Vector3 &goalPrior, Vector3 &goalEnd)
{

  PAtom *effectorEnd, *effectorPrior, *anchorEnd, *anchorPrior;
  Real angle1, angle2;

  protein->GetEndEffectors(side, effectorEnd, effectorPrior);
  protein->GetAnchors(side, anchorEnd, anchorPrior);

  if (side == start) {
    angle1 = ANGLE_C_N_CA;
    angle2 = ANGLE_CA_C_N;
  } else {
    angle1 = ANGLE_CA_C_N;
    angle2 = ANGLE_C_N_CA;
  }

  GetGoals(anchorPrior, anchorEnd, effectorPrior, effectorEnd, LENGTH_C_N, angle1, angle2,
     goalPrior, goalEnd);
}

ProteinSide PTools::GetGoals(PProtein *protein, Vector3 &goalPrior, Vector3 &goalEnd)
{
  if (EffectorAnchorDist(protein, start) > EffectorAnchorDist(protein, end)) {
    GetGoals(protein, start, goalPrior, goalEnd);
    return start;
  } else {
    GetGoals(protein, end, goalPrior, goalEnd);
    return end;
  }
}

Real PTools::EffectorAnchorDist(PProtein *loop, ProteinSide side)
{
  PAtom *effectorEnd, *effectorPrior, *anchorEnd, *anchorPrior;

  loop->GetEndEffectors(side, effectorEnd, effectorPrior);
  loop->GetAnchors(side, anchorEnd, anchorPrior);

  return effectorEnd->getPos().distance(anchorPrior->getPos());
}

int PTools::getBackboneAtomIndex(int resIndex, const string &atomID) {
  int incr;

  if (atomID == PID::N) {
    incr = 0;
  } else if (atomID == PID::C_ALPHA) {
    incr = 1;
  } else if (atomID == PID::C) {
    incr = 2;
  } else {
    PUtilities::AbortProgram("Invalid Atom ID to getBackboneAtomIndex");
  }

  return 3*resIndex + incr;
}

list<PCluster> PTools::GetClusters(const PConformationSpace &space, ConformationDistFn distFn, Real threshold)
{
  vector<PLightChain *> confArray(space.begin(), space.end());
  PCluster *tempArray = new PCluster[confArray.size()];

  list<PCluster> result;

  /* Build temporary array of clusters. */
  for(unsigned i = 0; i < confArray.size(); i++) {
    tempArray[i].AddConformation(confArray[i]);
  }

  /* Merge clusters if pointers are within threshold distance. */
  for(unsigned i = 0; i < confArray.size(); i++) {
    for(unsigned j = i + 1; j < confArray.size(); j++) {
      if (distFn(confArray[i], confArray[j]) < threshold) {
        tempArray[i].MergeWithCluster(tempArray[j]);
      }
    }
  }

  /* Prune out empty clusters. */
  for(unsigned i = 0; i < confArray.size(); i++) {
    if (!tempArray[i].empty()) {
      result.push_back(PCluster(tempArray[i]));
    }
  }

  delete[] tempArray;

  return result;
}

// forward direction Jacobian
//NOTE: I guess the Jacobian is for the last atom
void PTools::ComputeJacobian(PProtein *loop, double **Jac){
  int ndofs = loop->size()*2;
  int numRes = loop->size();
  vector<int> ind;
  int k = 0;
  for(int i=1;i<=loop->size();i++){
    ind.push_back(k);
    ind.push_back(k+1);
    k=k+2;
  }
  Vector3 p;
  p = Vector3(loop->getAtomPos(PID::BACKBONE,numRes*3-1).x,loop->getAtomPos(PID::BACKBONE,numRes*3-1).y,loop->getAtomPos(PID::BACKBONE,numRes*3-1).z);  
  ComputeJacobian(loop, ind, Jac, p, true);
}

void PTools::ComputeJacobian(PProtein* loop, double** Jac, Vector3& atom_pos)
{
	  int ndofs = loop->size()*2;
	  int numRes = loop->size();
	  vector<int> ind;
	  int k = 0;
	  for(int i=1;i<=loop->size();i++){
	    ind.push_back(k);
	    ind.push_back(k+1);
	    k=k+2;
	  }
//	  Vector3 p;
//	  p = Vector3(loop->getAtomPos(PID::BACKBONE,numRes*3-1).x,loop->getAtomPos(PID::BACKBONE,numRes*3-1).y,loop->getAtomPos(PID::BACKBONE,numRes*3-1).z);
	  ComputeJacobian(loop, ind, Jac, atom_pos, true);
}

//NOTE: Currently, just consider the forward direction
//void PTools::ComputeJacobian(PProtein* loop, int atom_index, double **Jac, bool forward)
//{
//	vector<int> index_list;
//	for( int i = 0; i < atom_index; i++)
//	{
//		index_list.push_back(i);
//	}
//	Vector3 p = loop->getAtomPos( PID::BACKBONE, atom_index);
//	ComputeJacobian( loop, index_list, Jac, p, forward);
//	return;
//}

void PTools::ComputeJacobian(PProtein *loop, vector<int>& ind, double **Jac, bool forward){
  int numRes = loop->size();
  Vector3 p;
  if (forward == true){
    p = Vector3(loop->getAtomPos(PID::BACKBONE,numRes*3-1).x,loop->getAtomPos(PID::BACKBONE,numRes*3-1).y,loop->getAtomPos(PID::BACKBONE,numRes*3-1).z);  
    ComputeJacobian(loop, ind, Jac, p, true);
  }
  else{
    p = Vector3(loop->getAtomPos(PID::BACKBONE,0).x,loop->getAtomPos(PID::BACKBONE,0).y,loop->getAtomPos(PID::BACKBONE,0).z);  
    ComputeJacobian(loop, ind, Jac, p, false);
  }
}

/*ind contains the atoms corresponding to the center of rotation*/
void PTools::ComputeJacobian(PProtein *loop, vector<int>& ind, double **Jac, Vector3& p,bool forward){
  double len;
  int numRes = loop->size();
  Vector3 q,r,t1,c1;
  int index;
  for(int i=0;i<ind.size();i++){
    index = 3*(ind[i]/2)+(ind[i]%2)+1;
    if (forward == true)
    {
      q = Vector3(loop->getAtomPos(PID::BACKBONE,index-1).x,loop->getAtomPos(PID::BACKBONE,index-1).y,loop->getAtomPos(PID::BACKBONE,index-1).z);
      r = Vector3(loop->getAtomPos(PID::BACKBONE,index).x,loop->getAtomPos(PID::BACKBONE,index).y,loop->getAtomPos(PID::BACKBONE,index).z);
    }
    else{
      q = Vector3(loop->getAtomPos(PID::BACKBONE,index).x,loop->getAtomPos(PID::BACKBONE,index).y,loop->getAtomPos(PID::BACKBONE,index).z);
      r = Vector3(loop->getAtomPos(PID::BACKBONE,index-1).x,loop->getAtomPos(PID::BACKBONE,index-1).y,loop->getAtomPos(PID::BACKBONE,index-1).z);
    }  
    len = (r-q).norm();
    t1=Vector3((r-q).x/len,(r-q).y/len,(r-q).z/len);  
    c1.setCross(t1,(p-r)); 
    Jac[1][i+1]=c1.x;
    Jac[2][i+1]=c1.y;
    Jac[3][i+1]=c1.z;
    Jac[4][i+1]=t1.x;
    Jac[5][i+1]=t1.y;
    Jac[6][i+1]=t1.z;
  }
}

//Calculates Singularity Robust Invers of a m x n matrix A 
void PTools::ComputePseudoInverse(double **A, int m, int n, double **Ainv){
  if (m<n){ //right inverse
    double AT[n+1][m+1];
    double **AF;
    double **AFinv;
    AF = (double **)calloc(m+1,sizeof(double*));
    AFinv = (double **)calloc(m+1,sizeof(double*));
    for(int i=0;i<=m;i++){
      AF[i]=(double*)calloc(m+1,sizeof(double));
      AFinv[i]=(double*)calloc(m+1,sizeof(double));
    }
    for(int i=1;i<=m;i++){
      for(int j=1;j<=n;j++){
        AT[j][i] = A[i][j];
      }
    }
    for(int i=1;i<=m;i++){
      for(int j=1;j<=m;j++){
        AF[i][j]=0.0;
        for(int k=1;k<=n;k++){
          AF[i][j] = AF[i][j]+A[i][k]*AT[k][j];
        }
      }
    }

    for(int i=1;i<=m;i++) AF[i][i]+=0.0004;//Singularity Robust Inverse

    PNumRoutines::nr_inverse(AF,m,AFinv);
    
    for(int i=1;i<=n;i++){
      for(int j=1;j<=m;j++){
        Ainv[i][j]=0.0;
        for(int k=1;k<=m;k++){
          Ainv[i][j] = Ainv[i][j]+AT[i][k]*AFinv[k][j];
        }
      }
    }
    
    free(AF);
    free(AFinv);
  }
  else if (m=n){ //simple inverse
    PNumRoutines::nr_inverse(A,m,Ainv);  
  }
  else if (m>n){ //left inverse
    double AT[n+1][m+1];
    double **AF;
    double **AFinv;
    AF = (double **)calloc(n+1,sizeof(double*));
    AFinv = (double **)calloc(n+1,sizeof(double*));
    for(int i=0;i<=n;i++){
      AF[i]=(double*)calloc(n+1,sizeof(double));
      AFinv[i]=(double*)calloc(n+1,sizeof(double));
    }
    for(int i=1;i<=m;i++){
      for(int j=1;j<=n;j++){
        AT[j][i] = A[i][j];
      }
    }
    for(int i=1;i<=n;i++){
      for(int j=1;j<=n;j++){
        AF[i][j]=0.0;
        for(int k=1;k<=m;k++){
          AF[i][j] = AF[i][j]+AT[i][k]*A[k][j];
        }
      }
    }
    for(int i=1;i<=n;i++) AF[i][i]+=0.0004;//Singularity Robust Inverse

    PNumRoutines::nr_inverse(AF,n,AFinv);
    for(int i=1;i<=n;i++){
      for(int j=1;j<=m;j++){
        Ainv[i][j]=0.0;
        for(int k=1;k<=n;k++){
          Ainv[i][j] = Ainv[i][j]+AFinv[i][k]*AT[k][j];
        }
      }
    }
    free(AF);
    free(AFinv);
  }
}


void PTools::ComputeNullSpace(double **Jac, int dim, bool SixDimensional, NullSpaceRet* Ret){
  if (SixDimensional) PNumRoutines::nr_svd(Jac, 6, dim, Ret->Sval, Ret->Svec);
  else PNumRoutines::nr_svd(Jac, 3, dim, Ret->Sval, Ret->Svec);
  double m = -1.0;
  for (int i = 1; i <= dim; i++)
  if (m < Ret->Sval[i])
  m = Ret->Sval[i];
  
  Ret->n_ns = 0;
  for (int i = 1; i <= dim; i++){
    if (Ret->Sval[i]/m < 0.0001){
      Ret->Sval[i] = 0.0;
      Ret->ns[Ret->n_ns+1] = i;
      Ret->n_ns++;
    }
  }
}
void PTools::ProjectOnNullSpace(PProtein *loop, vector<int> ind, bool forward, double ToProject[], double AfterProject[], bool sixDimensional){
	int loopsize = ind.size();
  double **Jac;
  NullSpaceRet * Ret = new NullSpaceRet;
  Jac  = (double**)calloc(6+1,sizeof(double*));
  for(int i=0;i<6+1;i++) Jac[i] = (double*)calloc(loopsize+1,sizeof(double));

  Ret->ns = (int*) new int[loopsize+1];
  Ret->Sval = (double*)calloc(loopsize+1,sizeof(double));
  Ret->Svec = (double**)calloc(loopsize+1,sizeof(double*));
  for(int i=0;i<loopsize+1;i++) (Ret->Svec)[i] = (double*)calloc(loopsize+1,sizeof(double));
  ComputeJacobian(loop, ind, Jac, forward);
  ComputeNullSpace(Jac, loopsize, sixDimensional, Ret);
  double U[loopsize+1][Ret->n_ns+1];
  double UTranspose[Ret->n_ns+1][loopsize+1];
  double V[loopsize+1][loopsize+1];
  for(int j = 1;j<=loopsize;j++){
    for(int k = 1;k<=Ret->n_ns;k++){
      U[j][k]=Ret->Svec[j][Ret->ns[k]];
      UTranspose[k][j]=U[j][k];
    }
  }
  for(int j = 1;j<=loopsize;j++){
    for(int k = 1;k<=loopsize;k++){
      for(int t = 1;t<=Ret->n_ns;t++) V[j][k] = U[j][t]*UTranspose[t][k];
    }
  }
  int proj_ind[Ret->n_ns+1];
  double proj[Ret->n_ns+1];
  double projf[Ret->n_ns+1];
  int t=1;
  for(int j=1;j<=Ret->n_ns;j++){
    proj[j] = 0.0;
    for(int k = 1;k<=loopsize;k++) proj[j] += UTranspose[j][k]*ToProject[k];
    if (proj[j] > 0.0){
      proj_ind[t] = j;
      projf[t]=proj[j];
      t++;
    }
  }
  int tot_k;
  tot_k = t-1;
  for(int j = 1; j <= loopsize; j++){
    AfterProject[j]=0.0;
    for(int k=1;k<=tot_k;k++) AfterProject[j] = AfterProject[j]+U[j][proj_ind[k]]*projf[k];
  }
free(Jac);
free(Ret->Sval);
free(Ret->Svec);
  delete Ret->ns;
  delete Ret;
	
}
void PTools::ProjectOnNullSpace(PProtein *loop, vector<int> ind, bool forward, double ToProject[], double AfterProject[]){
	PTools::ProjectOnNullSpace(loop, ind,forward, ToProject, AfterProject, true);
}
  
  
  

vector<vector<CDof> > PTools::GetBBDofs(vector<PProtein *> loops){
  vector<vector<CDof> > Dofs;
  for(int i=0;i<loops.size();i++){
    vector<CDof> Dof;
    for(int j=0;j<loops[i]->size()*2;j++){
      CDof temp;
      temp.dir = forward;
      temp.blockType = PID::BACKBONE;
      temp.DOF_index = j; 
      Dof.push_back(temp);
    }
    Dofs.push_back(Dof);
  }
  return Dofs; 
}

double PTools::RMSDCalpha (PProtein *protein0, PProtein *protein1, int loopstart, int loopend){
  double rmsd=0.0;
  Vector3 temp;
  int num_atoms=0;
  for(int j=loopstart; j<=loopend; j++){
    PResidue *res0 = protein0->getResidue(j);
    PResidue *res1 = protein1->getResidue(j);
    PAtom * atom0 = res0->getAtom(PID::C_ALPHA);
    PAtom * atom1 = res1->getAtom(PID::C_ALPHA);
    ++num_atoms;
    temp = (atom0->getPos() - atom1->getPos());
    rmsd += (temp.norm())*(temp.norm());
  }
  return sqrt(rmsd/num_atoms);
}

double PTools::RMSDBackbone(PProtein *protein0, PProtein *protein1, int loopstart, int loopend){
  double rmsd=0.0;
  Vector3 temp;
  int num_atoms=0;
  for(int j=loopstart; j<=loopend; j++){
    PResidue *res0 = protein0->getResidue(j);
    PResidue *res1 = protein1->getResidue(j);
    vector<PAtom *> result0;
    result0.push_back(res0->getAtom(PID::N));
    result0.push_back(res0->getAtom(PID::C_ALPHA));
    result0.push_back(res0->getAtom(PID::C));
    vector<PAtom *> result1;
    result1.push_back(res1->getAtom(PID::N));
    result1.push_back(res1->getAtom(PID::C_ALPHA));
    result1.push_back(res1->getAtom(PID::C));
    num_atoms = num_atoms + result1.size();
    for(int k=0; k<=2; k++){
      PAtom * atom0 = result0[k];
      PAtom * atom1 = result1[k];
      temp = (atom0->getPos() - atom1->getPos());
      rmsd += (temp.norm())*(temp.norm());
    }
  }
  return sqrt(rmsd/num_atoms);
}

double PTools::RMSDAllAtom(PProtein *protein0, PProtein *protein1, int loopstart, int loopend){
  double rmsd=0.0;
  Vector3 temp;
  int num_atoms = 0;
  for(int j=loopstart; j<=loopend; j++){
    PResidue *res0 = protein0->getResidue(j);
    PResidue *res1 = protein1->getResidue(j);
    vector<PAtom *> * result0 = protein0->getResidue(j)->getAtoms();
    vector<PAtom *> * result1 = protein1->getResidue(j)->getAtoms();
    num_atoms = num_atoms + (*result1).size();
    for(int k=0; k<protein0->getResidue(j)->getAtoms()->size(); k++){
      PAtom * atom0 = (*result0)[k];
      PAtom * atom1 = (*result1)[k];
      temp = (atom0->getPos() - atom1->getPos());
      rmsd += (temp.norm())*(temp.norm());
    }
  }
  return sqrt(rmsd/num_atoms);
}

bool PTools::CopyBackbone(PProtein *lpS, PProtein *lpD){
	return CopyBackbone(lpS, lpD, 0, 0, lpS->size());
	
}

bool PTools::CopyBackbone(PProtein *lpS, PProtein *lpD, int startS, int startD, int numRes){
	int endS, endD;
	double dist=lpS->getAtomAtRes(PID::N,startS)->getPos().distanceSquared(lpD->getAtomAtRes(PID::N,startD)->getPos()) + lpS->getAtomAtRes(PID::C_ALPHA,startS)->getPos().distanceSquared(lpD->getAtomAtRes(PID::C_ALPHA,startD)->getPos()) ;
	if(dist > 1E-8){
		cerr<<"BASES OF SOURCE AND DEST DO NOT MATCH"<<endl;
		return false;
	}
	endS = startS+numRes;
	endD = startD+numRes;
	if(endS > lpS->size()){
		cerr<<"SOURCE LOOP AND GIVEN INDICES DO NOT MATCH"<<endl;
	}
	if(endD > lpD->size()){
		cerr<<"DEST LOOP AND GIVEN INDICES DO NOT MATCH"<<endl;
	}
	
	Vector3 u0,u1,u2,u3,u4;
	double tang;
	int DOF_index, kS, kD;
	for(int k=0;k<numRes/3;k++){
		kS = startS+k*3;
		kD = startD+k*3;
		Vector3 endPriorGoal = lpS->getAtomAtRes(PID::C_ALPHA,kS+2)->getPos();
		Vector3 endGoal = lpS->getAtomAtRes(PID::C,kS+2)->getPos();
		Vector3 endNextGoal = lpS->getAtomAtRes(PID::O,kS+2)->getPos();
		int Res_indices[3]={kD,kD+1,kD+2};
		IKSolutions solns = PExactIKSolver::FindSolutions(lpD,Res_indices, &endPriorGoal, &endGoal, &endNextGoal);
		if(solns.size()==0) return false;
		double rmsd_f=1E8;
		int index=0;
		for(int i=0;i<solns.size();i++){
			lpD->MultiRotate(solns[i]);
			double rmsd=0.0;
			for(int j=0;j<(k+1)*3;j++){
				PResidue *res0 = lpS->getResidue(startS+j);
				PResidue *res1 = lpD->getResidue(startD+j);
				vector<PAtom *> result0;
				result0.push_back(res0->getAtom(PID::N));
				result0.push_back(res0->getAtom(PID::C_ALPHA));
				result0.push_back(res0->getAtom(PID::C));
				result0.push_back(res0->getAtom(PID::O));
				vector<PAtom *> result1;
				result1.push_back(res1->getAtom(PID::N));
				result1.push_back(res1->getAtom(PID::C_ALPHA));
				result1.push_back(res1->getAtom(PID::C));
				result1.push_back(res1->getAtom(PID::O));
				for(int t=0; t<=2; t++){
					PAtom * atom0 = result0[t];
					PAtom * atom1 = result1[t];
					rmsd += (atom0->getPos()).distanceSquared(atom1->getPos());
				}
			}
			if(rmsd < rmsd_f){
				rmsd_f = rmsd;
				index = i;
			}
			lpD->AntiMultiRotate(solns[i]);
		}
		lpD->MultiRotate(solns[index]);
	}
	for (int k=(numRes/3)*3;k<numRes;k++){
		kS = startS+k;
		kD = startD+k;
		u0 = lpD->getAtomAtRes(PID::N,kD)->getPos();
		u1 = lpD->getAtomAtRes(PID::C_ALPHA,kD)->getPos();
		u2 = lpD->getAtomAtRes(PID::C,kD)->getPos();
		u3 = lpS->getAtomAtRes(PID::C,kS)->getPos();
		tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
		DOF_index = kD*2;
		lpD->RotateChain(PID::BACKBONE,DOF_index,forward,-tang);
		
		u0 = lpD->getAtomAtRes(PID::C_ALPHA,kD)->getPos();
		u1 = lpD->getAtomAtRes(PID::C,kD)->getPos();
		u2 = lpD->getAtomAtRes(PID::O,kD)->getPos();
		u3 = lpS->getAtomAtRes(PID::O,kS)->getPos();
		tang = PMath::TorsionAngle(u2-u1,u3-u1,u1-u0);
		DOF_index = kD*2+1;
		lpD->RotateChain(PID::BACKBONE,DOF_index,forward,-tang);
	}
	return true;
}

void PTools::RandomNullSpacePerturb(PProtein *lp, double pert_mag){
    vector<PProtein*> lps;
    lps.push_back(lp);
    vector<vector<CDof> > Dofs = PTools::GetBBDofs(lps);
    RandomNullSpacePerturb(lp, Dofs, pert_mag);
}

void PTools::RandomNullSpacePerturb(PProtein *lp, vector<vector<CDof> > Dofs, double pert_mag){
	vector<int> JacInd;
	int i=0;
	int loopsize = Dofs[i].size();
//	srand(time(NULL));
	for(int j=0;j<loopsize;j++){
		if((Dofs[i])[j].blockType==PID::BACKBONE){
		JacInd.push_back((Dofs[i])[j].DOF_index);
		}
	}
	double derivb[loopsize+1];
	for(int j=1;j<=loopsize;j++) derivb[j]=0.0;
	double yb[loopsize+1];
	
	for(int j=1;j<=loopsize;j++)
		yb[j]=((double)rand()/(double)RAND_MAX)*pert_mag-(pert_mag)/2.0;
		
	PTools::ProjectOnNullSpace(lp,JacInd, true, yb, derivb);
	for(int j=0;j<Dofs[i].size();j++)
		lp->RotateBackbone(j,forward,derivb[j+1]);
} 

IKSolution PTools::CloseAlmostClosedLoop(PProtein *lp, Vector3 endPriorGoal,Vector3 endGoal, Vector3 endNextGoal){
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
		cerr<<"No solution found in CloseAlmostClosedLoop"<<endl;
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

int PTools::gsl_test(){
        double x=5.0;
        double y = gsl_sf_bessel_J0 (x);
        printf ("J0(%g) = %.18e\n", x, y);
        return 0;
}

