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

#include "PSampMethods.h"

IKSolution PSampMethods::RandAndIKClose(PProtein *loop, bool clash_free){
  string s;
  return RandAndIKClose(loop,s,clash_free);
}

IKSolution PSampMethods::RandAndIKClose(PProtein *loop,const string &pdbFileName, bool clash_free) {
	Vector3 goalPE = loop->getAtomAtRes(PID::C_ALPHA,loop->size()-1)->getPos();
	Vector3 goalE = loop->getAtomAtRes(PID::C,loop->size()-1)->getPos();
	Vector3 goalEN = loop->getAtomAtRes(PID::O,loop->size()-1)->getPos();
	loop->RandomizeDOFs(forward);
	int counter = 0;
	int Res_indices[3];
	Res_indices[2] = loop->size()-1;
	IKSolutions sols;
	for(int i=1;i<=Res_indices[2]-1;i++){
		Res_indices[1]=i;
		for(int j=0;j<=Res_indices[1]-1;j++){
			Res_indices[0] = j;
			sols = PExactIKSolver::FindSolutions(loop,Res_indices,&goalPE,&goalE,&goalEN);
			for(int k=0;k<sols.size();k++){
				loop->MultiRotate(sols[k]);
				if (clash_free == true){
				if (!loop->InAnyCollision()){
				if (pdbFileName.size()>0)PDBIO::writeToFile(loop->getTopLevelProtein(), pdbFileName);
				return sols[k];
				}
				else loop->AntiMultiRotate(sols[k]);
				}
				else {
				if (pdbFileName.size()>0) PDBIO::writeToFile(loop->getTopLevelProtein(), pdbFileName);
				return sols[k];
				}
			}
		}
	}
	IKSolution solr;
	solr.clear();
	return solr;  
}  

IKSolutions PSampMethods::PermuteIK(PProtein *loop, int num_wanted) {
	int Res_indices[3];
	IKSolutions ret;
	while(1){
		int loopsize = loop->size()-1;
		Res_indices[2] = 2+(int)((loopsize-3)*(double)rand()/(RAND_MAX+1.0));
		Res_indices[1] = 1+(int)((Res_indices[2]-1)*(double)rand()/(RAND_MAX+1.0));
		Res_indices[0] = 0+(int)((Res_indices[1]-1)*(double)rand()/(RAND_MAX+1.0));
		cout<<Res_indices[0]<<","<<Res_indices[1]<<","<<Res_indices[2]<<endl;
		IKSolutions sols = PExactIKSolver::FindSolutions(loop,Res_indices);
		for(int k=0;k<sols.size();k++){
			ret.push_back(sols[k]);
			if (ret.size()>=num_wanted) return ret;
		}
	}
	return ret;
}

vector<PProtein*> PSampMethods::DeformSampleBackbone(PProtein *protein, int loopSid, int loopEid, int num_wanted, double deform_mag) {
	vector<PProtein*> ret;
	for (int i=1; i<=num_wanted; ++i){
		PProtein *temp = protein->Clone();
		PProtein *loop = new PProtein(temp, loopSid, loopEid);
		
		Vector3 AtomCAP,AtomOP,AtomCP;
		AtomCAP=loop->getAtomAtRes(PID::C_ALPHA,loop->size()-1)->getPos();
		AtomCP=loop->getAtomAtRes(PID::C,loop->size()-1)->getPos();
		AtomOP=loop->getAtomAtRes(PID::O,loop->size()-1)->getPos();
		PTools::RandomNullSpacePerturb(loop, deform_mag);
		IKSolutions solns;
		solns.push_back(PTools::CloseAlmostClosedLoop(loop, AtomCAP,AtomCP, AtomOP));
		if(solns.size()>0)loop->MultiRotate(solns[0]);
		Vector3 AtomCAN,AtomON,AtomCN;
		AtomCAN=loop->getAtomAtRes(PID::C_ALPHA,loop->size()-1)->getPos();
		AtomCN=loop->getAtomAtRes(PID::C,loop->size()-1)->getPos();
		AtomON=loop->getAtomAtRes(PID::O,loop->size()-1)->getPos();
		double mag=0.0;
		mag=AtomCAN.distanceSquared(AtomCAP)+AtomCN.distanceSquared(AtomCP)+AtomON.distanceSquared(AtomOP);
		if(mag > 0.0001){
			i--;
			delete loop;
			temp->Obliterate();
			continue;
		}
		
		if(loop->InAnyCollision()){
			i--;
			delete loop;
			temp->Obliterate();
		}
		else{
			ret.push_back(temp);
		//	delete loop;
		}
	}
	return ret;
}

//NOTE: Yajia added this function
PProtein* PSampMethods::DeformSampleBackbone(PProtein *protein, int loopSid, int loopEid, double deform_mag) {
	PProtein* ret;
	for (int i=1; i<=1; ++i){
		PProtein *temp = protein->Clone();
		PProtein *loop = new PProtein(temp, loopSid, loopEid);

		Vector3 AtomCAP,AtomOP,AtomCP;
		AtomCAP=loop->getAtomAtRes(PID::C_ALPHA,loop->size()-1)->getPos();
		AtomCP=loop->getAtomAtRes(PID::C,loop->size()-1)->getPos();
		AtomOP=loop->getAtomAtRes(PID::O,loop->size()-1)->getPos();
		PTools::RandomNullSpacePerturb(loop, deform_mag);
		IKSolutions solns;
		solns.push_back(PTools::CloseAlmostClosedLoop(loop, AtomCAP,AtomCP, AtomOP));
		if(solns.size()>0)loop->MultiRotate(solns[0]);
		Vector3 AtomCAN,AtomON,AtomCN;
		AtomCAN=loop->getAtomAtRes(PID::C_ALPHA,loop->size()-1)->getPos();
		AtomCN=loop->getAtomAtRes(PID::C,loop->size()-1)->getPos();
		AtomON=loop->getAtomAtRes(PID::O,loop->size()-1)->getPos();
		double mag=0.0;
		mag=AtomCAN.distanceSquared(AtomCAP)+AtomCN.distanceSquared(AtomCP)+AtomON.distanceSquared(AtomOP);
		if(mag > 0.0001){
			i--;
			delete loop;
			temp->Obliterate();
			continue;
		}

		if(loop->InAnyCollision()){
			i--;
			delete loop;
			temp->Obliterate();
		}
		else{
//			ret.push_back(temp);
			ret = temp;
		//	delete loop;
		}
	}
	return ret;
}

/* The following is for seed sampling */
void nChoose2 (int n, int id, int *num1, int *num2) {
        int m1=0, m2=0;
        ++id;
        for (int i=1; i<n; ++i) {
                id -= n-i;
                if (id>0) {
                        ++m1;
                }
                else {
                        m2 = m1+id+n-i;
                        break;
                }
        }
        *num1 = m1;
        *num2 = m2;
}

// Sample a value randomly according to uniform distribution in interval A and B
double uniformSample (double A, double B) {
        if (A==B)
                return A;

        double upper = (A>B?A:B);
        double lower = (A>B?B:A);
        double i = rand()%100;
        return i*(upper-lower)/100+lower;
}

// Sample an angle value between [-180,180) according to a distribution distri
double sampleAngle (const vector<double> &distri) {
        double angle=0;
        if (distri.size()<=1)
                return (rand()%360)-180;
        else {
                double num = (rand()%10000)/(double)10000; // num is accurate up to 4 decimal points
                for (int i=0; i<distri.size(); ++i) {
                        if (num<=distri[i]) { // generate a angle value in interval i
                                double lower = i*(360/(double)distri.size())-180;
                                double upper = (i+1)*(360/(double)distri.size())-180;
                                double range = upper-lower;
                                angle = (rand()%10000)/(double)10000*range+lower;
                                break;
                        }
                        else {
                                num -= distri[i];
                        }
                } 
        }
        return angle;
}

class forwardOpenLoopNode {
  public:
        int Rid;
        string Type; // "Psi" and "Phi"
        Real Value; // angle value for Psi and Phi node
        int LastVisitedChildId;
        int LoopSize; // = the number of residues to be collision free
        // Root node is defined as (Rid=-1, Type="Psi", Value=-1, LastVistedChildId=-1, LoopSize=...)
        forwardOpenLoopNode (int i, string s, Real v, int id, int size) {
                Rid = i;
                Type = s;
                Value = v;
                LastVisitedChildId = id;
                LoopSize = size;
        }

        forwardOpenLoopNode (int i, string s, Real v, int size) {
                Rid = i;
                Type = s;
                Value = v;
                LastVisitedChildId = -1;
                LoopSize = size;
        }

        forwardOpenLoopNode (int i, string s, int size) {
                Rid = i;
                Type = s;
                Value = -1;
                LastVisitedChildId = -1;
                LoopSize = size;
        }

        bool equal (forwardOpenLoopNode n) {
                if (Rid==n.Rid && Type==n.Type && Value==n.Value && LastVisitedChildId==n.LastVisitedChildId && LoopSize==n.LoopSize)
                        return true;
                return false;
        }

        int getDof () {
                int dof;
                if (Type=="Phi")
                        dof = Rid*2;
                else // child is Psi
                        dof = Rid*2+1;
                return dof;
        }

        int getChildDof() {
                int dof;
                if (Type=="Phi") // child is Psi
                        dof = Rid*2+1;
                else { // child is Phi
                        if (Rid>=LoopSize-1) // no more child residue
                                dof = -1;
                        else
                                dof = (Rid+1)*2;
                }
                return dof;
        }

        int getChildRid() {
                int rid;
                if (Type=="Phi")
                        rid = Rid;
                else { // child is Phi
                        if (Rid>=LoopSize-1) // no more child residue
                                rid = -1;
                        else
                                rid = Rid+1;
                }
                return rid;
        }

        forwardOpenLoopNode* nextChild (vector<Real> *PhiValues, vector<Real> *PsiValues) {
                forwardOpenLoopNode *nc;
                int childId = LastVisitedChildId+1;
                if (Type=="Psi") { // next child is Phi of the next residue
                        if (Rid>=LoopSize-1) // no more residue
                                nc = NULL;
                        else {
                                if (PhiValues[Rid+1].size()<=childId)
                                        nc = NULL;
                                else {
                                        nc = new forwardOpenLoopNode(Rid+1,"Phi",PhiValues[Rid+1].at(childId),LoopSize);
                                        LastVisitedChildId = childId;
                                }
                        }
                }
                else { // next child is Psi of the same residue
                        if (PsiValues[Rid].size()<=childId)
                                nc = NULL;
                        else {
                                nc = new forwardOpenLoopNode(Rid,"Psi",PsiValues[Rid].at(childId),LoopSize);
                                LastVisitedChildId = childId;
                        }
                }
                return nc;
        } // end method nextChild

        int generateChildren(vector<Real> *PhiValues, vector<Real> *PsiValues, PPhiPsiDistribution &PhiPsiDistri) {
                if (Type=="Psi") {
                        if (Rid>=LoopSize-1)
                                return 0;
                        else
                                PhiValues[Rid+1].clear();
                }
                else
                        PsiValues[Rid].clear();
                int angleNum = 360/20;
                if (PhiPsiDistri.isEmpty()) {
                        // Take samples uniformly distributed over [0,359]
                        Real sample;
                        for (int i=0; i<angleNum; ++i) {
                                sample = uniformSample(0,359);
                                if (Type=="Psi") {
                                        PhiValues[Rid+1].push_back(sample);
                                }
                                else {
                                        PsiValues[Rid].push_back(sample);
                                }
                        } // end for
                }
                else {
                        double sample;
                        vector<double> distri;
                        if (Type=="Psi") {
                                distri = PhiPsiDistri.getPhiDistribution();
                        }
                        else {
                                distri = PhiPsiDistri.getPsiDistribution(Value);
                        }
                        for (int i=0; i<angleNum; ++i) {
				sample = sampleAngle(distri);
                                if (Type=="Psi") { // child type is Phi
                                        PhiValues[Rid+1].push_back(sample);
                                }
                                else {
                                        PsiValues[Rid].push_back(sample);
                                }
                        }
                }

                return 1;
        } // end method generateChildren

        void print () {
                cout << "node(" << Rid << "," << Type << "," << Value << "," << LastVisitedChildId << "," << LoopSize << ")" << endl;
        }

}; // end class forwardOpenLoopNode

class backwardOpenLoopNode {
  public:
        int Rid;
        string Type; // "Psi" and "Phi"
        Real Value; // angle value for Psi and Phi node
        int LastVisitedChildId;
        int LoopSize;
        // Root node is defined as (Rid=LoopSize, Type="Phi", Value=-1, LastVistedChildId=-1, LoopSize=...)

        backwardOpenLoopNode (int i, string s, Real v, int id, int loopsize) {
                Rid = i;
                Type = s;
                Value = v;
                LastVisitedChildId = id;
                LoopSize = loopsize;
        }

        backwardOpenLoopNode (int i, string s, Real v, int loopsize) {
                Rid = i;
                Type = s;
                Value = v;
                LastVisitedChildId = -1;
                LoopSize = loopsize;
        }

        backwardOpenLoopNode (int i, string s, int loopsize) {
                Rid = i;
                Type = s;
                Value = -1;
                LastVisitedChildId = -1;
                LoopSize = loopsize;
        }

        bool equal (backwardOpenLoopNode n) {
                if (Rid==n.Rid && Type==n.Type && Value==n.Value && LastVisitedChildId==n.LastVisitedChildId)
                        return true;
                return false;
        }

        int getDof () { // same as in forwardOpenLoopNode
                int dof;
                if (Type=="Phi")
                        dof = Rid*2;
                else // Type = "Psi"
                        dof = Rid*2+1;
                return dof;
        }

        int getChildDof() {
                int dof;
                if (Type=="Phi") { // child is Psi of the previous residue
                        if (Rid==0) // no child
                                dof = -1;
                        else
                                dof = (Rid-1)*2+1;
                }
                else { // child is Phi of the same residue
                        dof = Rid*2;
                }
                return dof;
        }

        int getChildRid () {
                int rid;
                if (Type=="Phi") { // child is Psi of the previous residue
                        if (Rid==0) // no child
                                rid = -1;
                        else
                                rid = Rid-1;
                }
                else { // child Phi of the same residue
                        rid = Rid;
                }
                return rid;
        }

        backwardOpenLoopNode* nextChild (vector<Real> *PhiValues, vector<Real> *PsiValues) {
        // Child of Rotamer node is Psi, child of Psi node is Phi, and child of Phi node is Rotamer.
                backwardOpenLoopNode *nc;
                int childId = LastVisitedChildId+1;
                if (Type=="Phi") { // child is Psi of the previous residue
                        if (Rid==0) // no more previous child
                                nc = NULL;
                        else {
                                if (PsiValues[Rid-1].size()<=childId)
                                        nc = NULL;
                                else {
                                        nc = new backwardOpenLoopNode(Rid-1,"Psi",PsiValues[Rid-1].at(childId),LoopSize);
                                        LastVisitedChildId = childId;
                                }
                        }
                }
                else { // child is Phi of the same residue
                        if (PhiValues[Rid].size()<=childId)
                                nc = NULL;
                        else {
                                nc = new backwardOpenLoopNode(Rid,"Phi",PhiValues[Rid].at(childId),LoopSize);
                                LastVisitedChildId = childId;
                        }
                }
                return nc;
        } // end method nextChild

        int generateChildren(vector<Real> *PhiValues, vector<Real> *PsiValues, PPhiPsiDistribution &PhiPsiDistri) {
                // Create no-collision intervals, and sample in those intervals
                if (Type=="Phi") { // children are Psi of the previous residue
                        if (Rid==0)
                                return 0;
                        else
                                PsiValues[Rid-1].clear();
                }
                else { // children are Phi of the same residue
                        PhiValues[Rid].clear();
                }
                int angleNum = 360/20;
                if (PhiPsiDistri.isEmpty()) { // Uniform distribution sampling
                        Real sample;
                        for (int i=0; i<angleNum; ++i) {
                                sample = uniformSample(0,359);
                                if (Type=="Psi") {
                                        PhiValues[Rid].push_back(sample);
                                }
                                else {
                                        PsiValues[Rid-1].push_back(sample);
                                }
                        } // end for
                }
                else { // sample angles according to the given distribution
                        vector<double> distri;
                        if (Type=="Psi") {
                                distri = PhiPsiDistri.getPhiDistribution(Value);
                        }
                        else {
                                distri = PhiPsiDistri.getPsiDistribution();
                        }
                        double sample;
                        for (int i=0; i<angleNum; ++i) {
                                sample = sampleAngle(distri);
                                if (Type=="Psi") {
                                        PhiValues[Rid].push_back(sample);
                                }
                                else {
                                        PsiValues[Rid-1].push_back(sample);
                                }
                        }
                }
                 
                return 1;
        } // end method generateChildren

        void print () {
                cout << "node(" << Rid << "," << Type << "," << Value << "," << LastVisitedChildId << "," << LoopSize << ")" << endl;
        }

}; // end class backwardOpenLoopNode 

vector<Real>* PSampMethods::generateForwardOpenLoop (PProtein* loopEntire, int collisionFreeResNum) {
        vector< vector<double> > empty_vector;
        PPhiPsiDistribution empty_distri(empty_vector);
        map<string,PPhiPsiDistribution> empty_map;
        empty_map.insert( make_pair("empty",empty_distri) );
        vector<string> aa_names;
        return generateForwardOpenLoop(loopEntire,collisionFreeResNum,empty_map,aa_names);
}

vector<Real>* PSampMethods::generateForwardOpenLoop (PProtein* loopEntire, int collisionFreeResNum, map<string,PPhiPsiDistribution> &map_distri, vector<string> &aa_names) {
// loopSid and loopEid is the first part of the loop which is supposed to collision-free
// The entire loopEntire ends at loopEid2
// The second part of loopEntire after loopEid also rotate when DoFs in loop change, but they are not guaranteed to be collision-free finally
// ** Finally, the residues 0 to collisionFreeResNum-1 are active for collision, but not for the rest residues. **


        bool single_distri = (map_distri.size()==1);
        bool random_distri = (map_distri.begin()->second).isEmpty();

        bool success=false;
        string child_aa_name;

        if (collisionFreeResNum>loopEntire->size())
                collisionFreeResNum = loopEntire->size();

        int loopSize = collisionFreeResNum;
        vector<Real> PhiValues[loopSize];
        vector<Real> PsiValues[loopSize];

        stack<forwardOpenLoopNode*> stk;
        int localDof;

        vector<Real> *angles = new vector<Real>;
        for (int i=loopSize*2-1; i>=0; --i)
                angles->push_back(0);

        forwardOpenLoopNode *root = new forwardOpenLoopNode(-1,"Psi",collisionFreeResNum); // root node
        if (single_distri)
                root->generateChildren(PhiValues,PsiValues,map_distri.begin()->second);
        else {
//              child_aa_name = loopEntire->getResidue(root->getChildRid())->getName();
                child_aa_name = aa_names[root->getChildRid()];
                root->generateChildren(PhiValues,PsiValues,map_distri.find(child_aa_name)->second);
        }
        stk.push(root);

        forwardOpenLoopNode *node;
        while (!stk.empty()) {
                node = stk.top();

                if (node->Rid==collisionFreeResNum-1 && node->Type=="Psi") {
                        success = true;
                        break;
                }

                else {
                        forwardOpenLoopNode *child = node->nextChild(PhiValues,PsiValues);
                        if (child==NULL) {
                                stk.pop();
                                // Rotate back the backbone by -node.Value 
                                if (node->Rid!=-1) {
                                        localDof = node->getDof();
                                        loopEntire->RotateBackbone(localDof,forward,-node->Value);
                                        (*angles)[localDof] += -node->Value;
                                }
                                // If Phi, detach the residue
                                if (node->Type=="Phi") {
                                        loopEntire->inactivateResidue(node->Rid,node->Rid);
                                }
				delete node;
                        }
                        else if (child->Type=="Phi") {
                                // Attach the new residue
                                loopEntire->activateResidue(child->Rid,child->Rid);
                                // Rotate the backbone if there is no collision
                                Real rotate_angle = child->Value;
                                if ( !random_distri) {
                                        Real curPhi = loopEntire->getResidue(child->Rid)->GetPhi();
                                        rotate_angle = curPhi-rotate_angle;
                                }
                                loopEntire->RotateBackbone(child->getDof(),forward,rotate_angle);
                                PProteinResidue *res = loopEntire->getResidue(child->Rid);
                                res->getAtom("O")->inactivate(); // O is not to be check collision with
                                if ( !res->getAtom("C")->InAnyCollision() &&
                                    ( res->getAtom("CB") == NULL || !res->getAtom("CB")->InAnyCollision() ) ) {
                                        stk.push(child);
                                        (*angles)[child->getDof()] += rotate_angle;
                                        if (single_distri)
                                                child->generateChildren(PhiValues,PsiValues,map_distri.begin()->second);
                                        else {
//                                              child_aa_name = loopEntire->getResidue(child->getChildRid())->getResourceName();
                                                child_aa_name = aa_names[child->getChildRid()];
                                                child->generateChildren(PhiValues,PsiValues,map_distri.find(child_aa_name)->second);
                                        }
                                }
                                else {
                                        loopEntire->RotateBackbone(child->getDof(),forward,-rotate_angle);
                                        loopEntire->inactivateResidue(child->Rid,child->Rid);
                                        delete child;
                                }
                        }
                        else { // child->Type=="Psi"
                                Real rotate_angle = child->Value;
                                if (!random_distri) {
                                        Real curPsi = loopEntire->getResidue(child->Rid)->GetPsi();
                                        rotate_angle = curPsi-rotate_angle;
                                }
                                loopEntire->RotateBackbone(child->getDof(),forward,rotate_angle);
                                PProteinResidue *res = loopEntire->getResidue(child->Rid);
				PAtom *currentO = res->getAtom("O");
				currentO->activate();
                                if ( !currentO->InAnyCollision() ) {
                                        if (child->Rid<loopEntire->size()-1) {
                                                res = loopEntire->getResidue(child->Rid+1);
                                                res->getAtom("N")->activate();
                                                res->getAtom("CA")->activate();
                                                if ( !res->getAtom("N")->InAnyCollision() && !res->getAtom("CA")->InAnyCollision() ) {
                                                        stk.push(child);
                                                        (*angles)[child->getDof()] += rotate_angle;

                                                        if (single_distri)
                                                                child->generateChildren(PhiValues,PsiValues,map_distri.begin()->second);
                                                        else if (child->getChildRid()>0) {
//                                                              child_aa_name = loopEntire->getResidue(child->getChildRid())->getName();
                                                                child_aa_name = aa_names[child->getChildRid()];
                                                                child->generateChildren(PhiValues,PsiValues,map_distri.find(child_aa_name)->second);
                                                        }
                                                }
                                                else {
                                                        loopEntire->RotateBackbone(child->getDof(),forward,-rotate_angle);
                                                        delete child;
                                                }
						res->getAtom("N")->inactivate();
						res->getAtom("CA")->inactivate();	
                                        }
                                        else {
                                                stk.push(child);
                                                (*angles)[child->getDof()] += rotate_angle;

                                                if (single_distri)
                                                        child->generateChildren(PhiValues,PsiValues,map_distri.begin()->second);
                                                else if (child->getChildRid() > 0) {
//                                                      child_aa_name = loopEntire->getResidue(child->getChildRid())->getName();
                                                        child_aa_name = aa_names[child->getChildRid()];
                                                        child->generateChildren(PhiValues,PsiValues,map_distri.find(child_aa_name)->second);
                                                }
                                        }
                                }
                                else {
                                        loopEntire->RotateBackbone(child->getDof(),forward,-rotate_angle);
					currentO->inactivate();
                                        delete child;
                                }
//				currentO->inactivate();
                        }

                }
        } // end while

        while (!stk.empty()) {
                delete stk.top();
                stk.pop();
        }

        if (!success) {
                angles->clear();
		loopEntire->inactivateResidue(0,loopEntire->size()-1); 
        }

        return angles;
}

vector<Real>* PSampMethods::generateBackwardOpenLoop (PProtein* loop) {
        vector< vector<double> > empty_vector;
        PPhiPsiDistribution empty_distri(empty_vector);
        map<string,PPhiPsiDistribution> empty_map;
        empty_map.insert( make_pair("empty",empty_distri) );
        vector<string> aa_names;
        return PSampMethods::generateBackwardOpenLoop(loop,empty_map,aa_names);
}

vector<Real>* PSampMethods::generateBackwardOpenLoop (PProtein* loop, map<string,PPhiPsiDistribution> &map_distri, vector<string> &aa_names) {


        bool single_distri = (map_distri.size()==1);
        bool random_distri = (map_distri.begin()->second).isEmpty();

        bool success=false;
        int loopSize = loop->size();
        string aa_name;

        vector<Real> PhiValues[loopSize];
        vector<Real> PsiValues[loopSize];
        stack<backwardOpenLoopNode*> stk;
        int localDof;

        vector<Real> *angles = new vector<Real>;
        for (int i=0; i<loopSize*2; ++i) {
                angles->push_back(0);
        }

        backwardOpenLoopNode *node = new backwardOpenLoopNode(loopSize,"Phi",-1,loopSize); // Root node
        if (single_distri)
                node->generateChildren(PhiValues,PsiValues,map_distri.begin()->second);
        else if (node->getChildRid()>=0) {
//              aa_name = loop->getResidue(node->getChildRid())->getName();
                aa_name = aa_names[node->getChildRid()];
                node->generateChildren(PhiValues,PsiValues,map_distri.find(aa_name)->second);
        }

        stk.push(node);

        while (!stk.empty()) {
                node = stk.top();
                if (node->Rid==0 && node->Type=="Phi") {
                        success = true;
                        break;
                }
                else {
                        backwardOpenLoopNode *child = node->nextChild(PhiValues,PsiValues);
                        if (child==NULL) {
                                stk.pop();

                                if (node->Rid!=loopSize) {
                                        loop->RotateBackbone(node->getDof(),backward,-node->Value);
                                        (*angles)[node->getDof()] += -node->Value;
                                }
                                if (node->Type=="Psi") { // detach the residue
                                        loop->inactivateResidue(node->Rid,node->Rid);
                                }
				delete node;
                        }
                        else if (child->Type=="Psi") {
                                // Attach a new residue
                                loop->activateResidue(child->Rid,child->Rid);
                                Real rotate_angle = child->Value;
                                if (!random_distri) {
                                        Real curPsi = loop->getResidue(child->Rid)->GetPsi();
                                        rotate_angle = curPsi-rotate_angle;
                                }
                                loop->RotateBackbone(child->getDof(),backward,rotate_angle);
                                PResidue *res = loop->getResidue(child->Rid);
                                if ( !res->getAtom("N")->InAnyCollision() &&
                                    ( res->getAtom("CB")==NULL || !res->getAtom("CB")->InAnyCollision()) ) {
                                        stk.push(child);
                                        (*angles)[child->getDof()] += rotate_angle;
                                        if (single_distri)
                                                child->generateChildren(PhiValues,PsiValues,map_distri.begin()->second);
                                        else if (child->getChildRid()>=0) {
//                                              aa_name = loop->getResidue(child->getChildRid())->getName();
                                                aa_name = aa_names[child->getChildRid()];
                                                child->generateChildren(PhiValues,PsiValues,map_distri.find(aa_name)->second);
                                        }
                                }
                                else {
                                        loop->RotateBackbone(child->getDof(),backward,-rotate_angle);
                                        loop->inactivateResidue(child->Rid,child->Rid);
                                        delete child;
                                }
                        }
                        else {// child->Type == "Phi"
                                Real rotate_angle = child->Value;
                                if (!random_distri) {
                                        Real curPhi = loop->getResidue(child->Rid)->GetPhi();
                                        rotate_angle = curPhi-rotate_angle;
                                }
                                if (child->Rid>0) {
                                        Real rotate_angle = child->Value;
                                        if (!random_distri) {
                                                Real curPhi = loop->getResidue(child->Rid)->GetPhi();
                                                rotate_angle = curPhi-rotate_angle;
                                        }
                                        loop->RotateBackbone(child->getDof(),backward,rotate_angle);
                                        PResidue *res = loop->getResidue(child->Rid-1);
                                        bool ok = true;
                                        PAtom *atom = res->getAtom("C");
                                        atom->activate();
                                        if (atom->InAnyCollision()) {
                                                ok = false;
					}
                                        atom->inactivate();
                                        if (ok) {
                                                atom = res->getAtom("O");
                                                atom->activate();
                                                if (atom->InAnyCollision()) {
                                                        ok = false;
						}
                                                atom->inactivate();
                                        }
                                        if (ok) {
                                                atom = res->getAtom("CA");
                                                if (atom!=NULL) {
                                                        atom->activate();
                                                        if (atom->InAnyCollision()) {
                                                                ok = false;
							}
                                                        atom->inactivate();
                                                }
                                        }
                                        if (ok) {
                                                stk.push(child);
                                                (*angles)[child->getDof()] += rotate_angle;
                                                if (single_distri) 
                                                        child->generateChildren(PhiValues,PsiValues,map_distri.begin()->second);
                                                else if (child->getChildRid()>=0) {
//                                                      aa_name = loop->getResidue(child->getChildRid())->getName();
                                                        aa_name = aa_names[child->getChildRid()];
                                                        child->generateChildren(PhiValues,PsiValues,map_distri.find(aa_name)->second);
                                                }
                                        }
                                        else {
                                                loop->RotateBackbone(child->getDof(),backward,-rotate_angle);
                                                delete child;
                                        }
                                }
                                else {
                                        stk.push(child);
                                        (*angles)[child->getDof()] += rotate_angle;
                                        if (single_distri)
                                                child->generateChildren(PhiValues,PsiValues,map_distri.begin()->second);
                                        else if (child->getChildRid() >= 0) {
//                                              aa_name = loop->getResidue(child->getChildRid())->getName();
                                                aa_name = aa_names[child->getChildRid()];
                                                child->generateChildren(PhiValues,PsiValues,map_distri.find(aa_name)->second);
                                        }
                                }
                        }
                }
        } // end while

        while (!stk.empty()) {
                delete stk.top();
                stk.pop();
        }
        if (!success) {
		loop->inactivateResidue(0,loop->size()-1); 
                angles->clear();
        }
        return angles;
}

Vector3 PSampMethods::computePreCpos (PResidue *res, Real angle, Real bondLength) {
//        Vector3 C_pos = res->getAtomPosition("C");
        Vector3 Ca_pos = res->getAtomPosition("CA");
        Vector3 N_pos = res->getAtomPosition("N");

        return PMath::ComputePos(Ca_pos,N_pos,angle,bondLength);
}

void PSampMethods::computeGoal (PResidue *res, SpaceRelationship *sr, Vector3 *endPriorG, Vector3 *endG, Vector3 *endNextG) {
// The goal residue is the one before res. sr is the space relationship between the goal residue and res in pdb file.
// Output is endPriorG, endG, and endNextG
// endPriorG is its CA atom position, endG is its C atom position, and endNextG is its O atom position.
// endG is computed "randomly" in the sense that the C1-N2-Ca2 plane might not be the same as in origial pdb file. But the length C1-N2 and angle C1-N2-Ca2 are the same.
// Then compute Ca1 and O1 according to the space relationship.
        Vector3 a,b,n,x;
        Matrix3 m;
        Vector3 Ca2_pos = res->getAtomPosition("CA");
        Vector3 N2_pos = res->getAtomPosition("N");
        Vector3 C2_pos = res->getAtomPosition("C");

        // Compute C1_pos
        Vector3 C1_pos = PSampMethods::computePreCpos(res,sr->angle_c_n_ca,sr->length_c_n);

        // Compute normal vector to the plane of N2_pos, Ca2_pos and C1_pos
        a = Ca2_pos - C1_pos;
        b = N2_pos - C1_pos;
        n = cross(a,b);

        // Compute Ca1_pos
        m = PMath::FindRotationMatrix(n,DtoR(-sr->angle_ca_c_n));
        x = m*b;
        x.inplaceNormalize();
        x.inplaceScale(sr->length_ca_c);
        Vector3 Ca1_pos = C1_pos + x;

        // Compute O1_pos
        m = PMath::FindRotationMatrix(n,DtoR(sr->angle_o_c_n));
        x = m*b;
        x.inplaceNormalize();
        x.inplaceScale(sr->length_o_c);
        Vector3 O1_pos = C1_pos + x;

        *endPriorG = Ca1_pos;
        *endG = C1_pos;
        *endNextG = O1_pos;
}

bool PSampMethods::IKClose (PProtein* move_loop, Vector3 endPriorG, Vector3 endG, Vector3 endNextG) {
// IKClose enumerates all possible residue sets, and randomly chooses one to do exact IK close. If it doesn't succeed, will choose another set randomly and unrepeatedly.
// move_loop is the protein loop to be closed.
// endPriorG is the target last CA atom position.
// endG is the target last C atom position.
// endNextG is the target last O atom position.
                /// Exact IK close ///
                int res_indices[3];
                bool no_sol = true;
                int pattern_id;
                int nC2 = (move_loop->size()-1)*(move_loop->size()-2)/2; // n = loop->size()-1
                bool *visited = new bool[nC2];
                for (int m=0; m<nC2; ++m)
                        visited[m] = false;
                for (int m=1; no_sol && m<=nC2; ++m) {
                        pattern_id = rand()%nC2;
                        if (visited[pattern_id]==true) {
                                --m;
                                continue;
                        }
                        else
                                visited[pattern_id] = true;
                        nChoose2(move_loop->size()-1,pattern_id,res_indices,res_indices+1);
                        res_indices[2] = move_loop->size()-1;
                        IKSolutions sols = PExactIKSolver::FindSolutions(move_loop,res_indices,&endPriorG,&endG,&endNextG);

                        if (sols.size()>0) {
                                int sol_id = rand()%sols.size();
                                vector<ChainMove> sol = sols[sol_id];
                                move_loop->MultiRotate(sol);
                                no_sol = false;
                        } // end if
                } // end for
		delete[] visited;
        return no_sol;
}

void PSampMethods::UpdateProtein(PProtein *p, vector<Real> *positiveAngles, vector<Real> *negativeAngles, bool forwardDir) {
        for (int i=0; i<positiveAngles->size(); ++i) {
                Real angle = positiveAngles->at(i)-negativeAngles->at(i);
                if (forwardDir) {
                        p->RotateBackbone(i,forward,angle);
		}
                else
                        p->RotateBackbone(i,backward,angle);
        }
}

void PSampMethods::UpdateProtein(PProtein *p, vector<Real> *angles, bool forwardDir) {
        for (int i=0; i<angles->size(); ++i) {
                if (forwardDir)
                        p->RotateBackbone(i,forward,-angles->at(i));
                else
                        p->RotateBackbone(i,backward,-angles->at(i));
        }
}

double PSampMethods::computeMaxLength (int proteinSize) {
        double ca_n, ca_ca, c_ca;
        double ca_c_n_radian = ANGLE_CA_C_N*deg2rad;
        double c_n_ca_radian = ANGLE_C_N_CA*deg2rad;
        double ca_n_c_radian, ca_n_ca_radian;

        // Compute the max distance from CA_i to CA_(i+1)       
        ca_n = sqrt(pow(LENGTH_CA_C,2)+pow(LENGTH_C_N,2)-2*LENGTH_CA_C*LENGTH_C_N*cos(ca_c_n_radian));
        ca_n_c_radian = asin(LENGTH_CA_C*sin(ca_c_n_radian)/ca_n);
        ca_n_ca_radian = ca_n_c_radian+c_n_ca_radian;
        ca_ca = sqrt(pow(ca_n,2)+pow(LENGTH_N_CA,2)-2*ca_n*LENGTH_N_CA*cos(ca_n_ca_radian));

        // Compute distance from C_i to CA_(i+1)
        c_ca = sqrt(pow(LENGTH_C_N,2)+pow(LENGTH_N_CA,2)-2*LENGTH_C_N*LENGTH_N_CA*cos(c_n_ca_radian));

        return (proteinSize-1)*ca_ca+c_ca+ca_n;
}

/*
        vector<PProtein*> loops = SeedSampleBackboneLoopOnly(protein,loopSid,loopEid,num_wanted);
        vector<PProtein*> merged;
        for (int i=0; i<loops.size(); ++i) {
                PProtein *p = PSampMethods::MergeProtein(loops[i],protein,loopSid);
                merged.push_back(p);
		delete loops[i];
        }
        return merged;
}
*/

vector<PProtein*> PSampMethods::SeedSampleBackbone (PProtein *protein, int loopSid, int loopEid, int num_wanted) {
	vector<PProtein*> loops = SeedSampleBackboneLoopOnly(protein,loopSid,loopEid,num_wanted);
	vector<PProtein*> result;
	int loopsize = loopEid-loopSid+1;
	for (int i=0; i<loops.size(); ++i) {
		PProtein *p = PTools::CreateSlimProtein(protein,loopSid,loopEid);
		for (int j=0; j<loopsize; ++j) {
			p->getResidue(loopSid+j)->setPositions(loops[i]->getResidue(j)->getSpec().Atom_Positions);
		}
		result.push_back(p);
		loops[i]->Obliterate();
	}
	return result;
}

//PProtein* PSampMethods::SeedSampleBackbone( PProtein* protein, int loopSid, int loopEid)
//{
//	vector<PProtein*> loops = SeedSampleBackboneLoopOnly(protein,loopSid,loopEid,num_wanted);
//	vector<PProtein*> result;
//	int loopsize = loopEid-loopSid+1;
//	for (int i=0; i<loops.size(); ++i) {
//		PProtein *p = PTools::CreateSlimProtein(protein,loopSid,loopEid);
//		for (int j=0; j<loopsize; ++j) {
//			p->getResidue(loopSid+j)->setPositions(loops[i]->getResidue(j)->getSpec().Atom_Positions);
//		}
//		result.push_back(p);
//		loops[i]->Obliterate();
//	}
//	return result;
//}
			
/*
vector<PProtein*> PSampMethods::SeedSampleBackbone (PProtein *protein, int loopSid, int loopEid, map<string,PPhiPsiDistribution> &distri_map, int num_wanted) {
        vector<PProtein*> loops = SeedSampleBackboneLoopOnly(protein,loopSid,loopEid,distri_map,num_wanted);
        vector<PProtein*> merged;
        for (int i=0; i<loops.size(); ++i) {
                PProtein *p = PSampMethods::MergeProtein(loops[i],protein,loopSid);
                merged.push_back(p);
		delete loops[i];
        }
        return merged;
}
*/

vector<PProtein*> PSampMethods::SeedSampleBackbone (PProtein *protein, int loopSid, int loopEid, map<string,PPhiPsiDistribution> &distri_map, int num_wanted) {
	vector<PProtein*> loops = SeedSampleBackboneLoopOnly(protein,loopSid,loopEid,distri_map,num_wanted);
	vector<PProtein*> result;
	int loopsize = loopEid-loopSid+1;
	for (int i=0; i<loops.size(); ++i) {
		PProtein *p = PTools::CreateSlimProtein(protein,loopSid,loopEid);
		for (int j=0; j<loopsize; ++j) {
			p->getResidue(loopSid+j)->setPositions(loops[i]->getResidue(j)->getSpec().Atom_Positions);
		}
		result.push_back(p);
		loops[i]->Obliterate();
	}
	return result;
}

PProtein* PSampMethods::MergeProtein (PProtein *loop, PProtein* original_protein, int startRid) {
                int endRid = startRid+loop->size()-1;

                PProtein *temp_protein = new PProtein();
                PResidue *res;
                for (int rid=0; rid<startRid; ++rid) {
                        res = original_protein->getResidue(rid);
			PResidueSpec prs = res->getSpec();
                        temp_protein->AddResidue(res->getResourceName(),prs);
                }
                for (int rid=startRid; rid<=endRid; ++rid) {
                        res = loop->getResidue(rid-startRid);
			PResidueSpec prs = res->getSpec();
                        temp_protein->AddResidue(original_protein->getResidue(rid)->getResourceName(),prs);
//			temp_protein->AddResidue(res->getName(),res->getSpec());
                }
                for (int rid=endRid+1; rid<original_protein->size(); ++rid) {
                        res = original_protein->getResidue(rid);
			PResidueSpec prs = res->getSpec();
                        temp_protein->AddResidue(res->getResourceName(),prs);
                }
                temp_protein->finalize();
                return temp_protein;
}

PProtein* PSampMethods::MergeProteinByPdbId (PProtein *loop, PProtein* original_protein, int startPdbId) {
	int original_start_pdb_id = original_protein->getResidue(0)->getPdbId();
	int original_end_pdb_id = original_protein->getResidue(original_protein->size()-1)->getPdbId();

	int loop_start_pdb_id = loop->getResidue(0)->getPdbId();
	int loop_end_pdb_id = loop->getResidue(loop->size()-1)->getPdbId();
	int loop_pdb_id_diff = startPdbId - loop_start_pdb_id;
	int loop_new_end_pdb_id = loop_end_pdb_id + loop_pdb_id_diff;

	int cur_pdb_id = min(original_start_pdb_id,startPdbId);
	int final_pdb_id = max(original_end_pdb_id,loop_new_end_pdb_id);
	
	PProtein *temp_protein = new PProtein();
	PResidueSpec spec;
	PProteinResidue *res;
	int index;

	while (cur_pdb_id<=final_pdb_id) {
		if (cur_pdb_id>=startPdbId && cur_pdb_id<=loop_new_end_pdb_id) {
			index = loop->pdbIndexToLocalIndex(cur_pdb_id-loop_pdb_id_diff);
			if (index!=-1)
				res = loop->getResidue(index);
			else {
				++cur_pdb_id;
				continue;
			}
		}
		else {
			index = original_protein->pdbIndexToLocalIndex(cur_pdb_id);
			if (index!=-1)
				res = original_protein->getResidue(index);
			else {
				++cur_pdb_id;
				continue;
			}
		}
		spec = res->getSpec();
		temp_protein->AddResidue(res->getResourceName(),spec);
		temp_protein->getResidue(temp_protein->size()-1)->setPdbId(cur_pdb_id);
		++cur_pdb_id;
	}
	temp_protein->finalize();
	return temp_protein;
}	

/*
PProtein* PSampMethods::MergeProteinByPdbId (PProtein *loop, PProtein* original_protein, int startPdbId) {
	int endPdbId = startPdbId+loop->size()-1;
	int original_end = original_protein->getResidue(original_protein->size()-1)->getPdbId();

	PProtein *temp_protein = new PProtein();
	PResidueSpec prs;
	PResidue *original_cur_res, *loop_cur_res;
	original_cur_res = original_protein->getResidue(0);
	loop_cur_res = loop->getResidue(0);
	int original_cur_pdb_id = original_cur_res->getPdbId();
	int loop_cur_pdb_id = startPdbId;

	while (original_cur_pdb_id<=original_end && loop_cur_pdb_id<=endPdbId) {
		if (original_cur_pdb_id<loop_cur_pdb_id) {
			prs = original_cur_res->getSpec();
			temp_protein->AddResidue(original_cur_res->getResourceName(),prs);
			if (original_cur_pdb_id<original_end) {
				original_cur_res = original_cur_res->NextResidue();
				original_cur_pdb_id = original_cur_res->getPdbId();
			}
			else if (original_cur_pdb_id<loop_cur_pdb_id) {
				original_cur_pdb_id = endPdbId;
			}
			else {
				original_cur_pdb_id++;
			}
		}
		else {
			prs = loop_cur_res->getSpec();
			temp_protein->AddResidue(loop_cur_res->getResourceName(),prs);
			temp_protein->getResidue(temp_protein->size()-1)->setPdbId(loop_cur_pdb_id);
			if (loop_cur_pdb_id==original_cur_pdb_id) {
				if (original_cur_pdb_id<original_end) {
					original_cur_res = original_cur_res->NextResidue();
					original_cur_pdb_id = original_cur_res->getPdbId();
				}
			}
			if (loop_cur_pdb_id<endPdbId) {
				loop_cur_res = loop_cur_res->NextResidue();
				loop_cur_pdb_id++;
			}
			else if (loop_cur_pdb_id<original_cur_pdb_id) {
				loop_cur_pdb_id = original_end;
			}
			else {
				loop_cur_pdb_id++;
			}
		}
	}
	temp_protein->finalize();
	return temp_protein;
}
*/	

// Only the residue order, atom name and atom position are correct in the final output. The temperature, residue name and ect. might not be correct.
void PSampMethods::addSidechain (string loopFile, string boundaryFile, string scwrl3_path, string outLoopFile) {
        string command = scwrl3_path+"/scwrl3 -i "+loopFile+" -o temp.out -f "+boundaryFile+" > scwrl3.log";
        system(command.c_str());

        // re-format the output file so that it is in proper PDB format 
        ifstream fin("temp.out");
        ofstream fout(outLoopFile.c_str());
        string s, s1, s2=" 1.00                  ";
        while (!fin.eof()) {
                getline(fin,s);
                if (s.size()<30)
                        break;
                s.erase(30,1);
                s1 = s.substr(0,55);
                fout << s1 << s2 << endl;
        }
        fin.close();
        fout.close();
	
	system("rm temp.out");
	system("rm scwrl3.log");
}

vector<PProtein*> PSampMethods::SeedSampleBackboneLoopOnly (PProtein* original_protein, int loopSid, int loopEid, int num_wanted) {

        static int MIN_MOVE_LOOP_SIZE=4;
        static int MAX_TRIAL_PER_PAIR=50;
        static int C_THRESHOLD=0, MAX_ENDS_NUM=max(num_wanted*2,20);
        static double MIDDLE_SIZE_RATIO=0.5;
        static double MAX_LENGTH_DISCOUNT_RATIO = 1;
        int cfFrontEndLength=0, cfBackEndLength=0;
        PProtein *frontLoop, *backLoop, *loop, *move_loop, *protein;
        Vector3 endG, endPriorG, endNextG, frontLastAtomPos, backFirstAtomPos;
        SpaceRelationship *sr;
        bool workFrontEnd=false, okEndPair=true;
        double endsDistanceThreshold=0;
        vector<PProtein*> result;
        vector<Real> *angles;

        protein = PTools::CreateSlimProtein(original_protein,loopSid,loopEid);

        int loopSize = loopEid-loopSid+1;
        int middleSize = (int)(floor(MIDDLE_SIZE_RATIO*loopSize));
        bool split = middleSize>=MIN_MOVE_LOOP_SIZE;

        loop = new PProtein(protein,loopSid,loopEid);
        if (split) {
                cfFrontEndLength = (loopSize-middleSize)/2;
                cfBackEndLength = loopSize-middleSize-cfFrontEndLength;
                frontLoop = new PProtein(loop,0,loop->size()-cfBackEndLength-1);
                backLoop = new PProtein(loop,loop->size()-cfBackEndLength,loop->size()-1);
                move_loop = new PProtein(frontLoop,cfFrontEndLength,frontLoop->size()-1);
                endsDistanceThreshold = computeMaxLength(middleSize) * MAX_LENGTH_DISCOUNT_RATIO;
                sr = new SpaceRelationship(
                                original_protein->getResidue(loopEid-cfBackEndLength),
                                original_protein->getResidue(loopEid-cfBackEndLength+1));
        }
        else {
                move_loop = loop;
                PResidue *res = protein->getResidue(loopEid);
                endPriorG = res->getAtomPosition("CA");
                endG = res->getAtomPosition("C");
                endNextG = res->getAtomPosition("O");
        }

        while (result.size() < num_wanted) {
                if (split) {
                                // Randomize the loop
                                for (int i=0; i<2*frontLoop->size(); ++i)
                                        frontLoop->RotateBackbone(i,forward,rand()%360);
                                for (int i=0; i<2*backLoop->size(); ++i)
                                        backLoop->RotateBackbone(i,backward,rand()%360);

                                // Generate new initial forward and backward ends
                                loop->inactivateResidue(0,loop->size()-1);
                                bool gotends=false;                                                               
                                while (!gotends) {
                                        while (true) {                                                            
                                                angles = generateForwardOpenLoop(frontLoop,cfFrontEndLength);
                                                if (!angles->empty()) {                                           
							angles->clear();
                                                        delete angles;                                            
                                                        break;                                                    
                                                }                                                                 
                                        }                                                                         
                                        for (int i=0; i<10; ++i) {                                                
                                                angles = generateBackwardOpenLoop(backLoop);
                                                if (!angles->empty()) {                                           
                                                        gotends = true;                                           
							angles->clear();
                                                        break;                                                    
                                                }                                                                 
                                        }                                                                         
                                }

                                frontLastAtomPos = frontLoop->getResidue(cfFrontEndLength-1)->getAtomPosition("C");
                                backFirstAtomPos = backLoop->getResidue(0)->getAtomPosition("N");
                        okEndPair = true;
                        if (frontLastAtomPos.distance(backFirstAtomPos) > endsDistanceThreshold) {
                                okEndPair = false;
                        }
                         // Check if these two ends collide
                        if (okEndPair) {
                                move_loop->inactivateResidue(0,move_loop->size()-1);
                                okEndPair = !(backLoop->InAnyCollision());
                                move_loop->activateResidue(0,move_loop->size()-1);
                        }

                        if (okEndPair) {
                                PResidue *backHeadResInProtein = protein->getResidue(loopEid-cfBackEndLength+1);
                                PSampMethods::computeGoal(backHeadResInProtein,sr,&endPriorG,&endG,&endNextG);
                        }
                } // end if

                for (int tryNum=0; okEndPair && tryNum<MAX_TRIAL_PER_PAIR; ++tryNum) {
                        // Randomize backbone DOFs
                        vector<ChainMove> cms;
                        ChainMove cm;
                        for (int j=0; j<2*move_loop->size(); ++j) {
                                cm.blockType = PID::BACKBONE;
                                cm.dir = forward;
                                cm.DOF_index = j;
                                cm.degrees = rand()%360;
                                cms.push_back(cm);
                        }
                        move_loop->MultiRotate(cms);

                        // Close the loop
                        bool no_sol = PSampMethods::IKClose(move_loop,endPriorG,endG,endNextG);
                        // Check collision and avoid collision if less enough
                        if (no_sol) {
                                continue;
                        }       
                        if (!move_loop->InAnyCollision()) {
                        	result.push_back(loop->Clone());
                                break;
                        }
                } // end for
        } // end while  
                        
        if (split)      
                delete sr;      
        delete protein;         

        return result;
}

vector<PProtein*> PSampMethods::SeedSampleBackboneLoopOnly (PProtein* original_protein, int loopSid, int loopEid, map<string,PPhiPsiDistribution> &distri_map, int num_wanted) {

        if (distri_map.size()==20) {
                if (distri_map.count("ALA")==0 || distri_map.count("ARG")==0 || distri_map.count("ASN")==0
                || distri_map.count("ASP")==0 || distri_map.count("CYS")==0 || distri_map.count("GLN")==0
                || distri_map.count("GLU")==0 || distri_map.count("GLY")==0 || distri_map.count("HIS")==0
                || distri_map.count("ILE")==0 || distri_map.count("LEU")==0 || distri_map.count("LYS")==0
                || distri_map.count("MET")==0 || distri_map.count("PHE")==0 || distri_map.count("PRO")==0
                || distri_map.count("SER")==0 || distri_map.count("THR")==0 || distri_map.count("TRP")==0
                || distri_map.count("TYR")==0 || distri_map.count("VAL")==0 ) {
                        cerr << "Wrong distribution names. Please read the documentation carefull." << endl;
                        exit(1);
                }
        }

        static int MIN_MOVE_LOOP_SIZE=4;
        static int MAX_TRIAL_PER_PAIR=50;
        static int C_THRESHOLD=0, MAX_ENDS_NUM=max(num_wanted*2,20);
        static double MIDDLE_SIZE_RATIO=0.5;
        static double MAX_LENGTH_DISCOUNT_RATIO = 1;
        int cfFrontEndLength=0, cfBackEndLength=0;
        PProtein *frontLoop, *backLoop, *loop, *move_loop, *protein;
        Vector3 endG, endPriorG, endNextG, frontLastAtomPos, backFirstAtomPos;
        SpaceRelationship *sr;
        bool workFrontEnd=false, okEndPair=true;
        double endsDistanceThreshold=0;
        vector<PProtein*> result;
        vector<Real> *angles;

        protein = PTools::CreateSlimProtein(original_protein,loopSid,loopEid);
        
        int loopSize = loopEid-loopSid+1;
        int middleSize = (int)(floor(MIDDLE_SIZE_RATIO*loopSize));
//        bool split = middleSize>=(MIN_MOVE_LOOP_SIZE+1);
        bool split = middleSize>=MIN_MOVE_LOOP_SIZE;

	vector<string> front_aa_names, back_aa_names, move_aa_names;
	map<string,PPhiPsiDistribution>::iterator map_iter;
	string aa_name;
	double phi, psi, curPhi, curPsi;
        
        loop = new PProtein(protein,loopSid,loopEid);
        if (split) {
                cfFrontEndLength = (loopSize-middleSize)/2;
                cfBackEndLength = loopSize-middleSize-cfFrontEndLength;
                frontLoop = new PProtein(loop,0,loop->size()-cfBackEndLength-1);
                backLoop = new PProtein(loop,loop->size()-cfBackEndLength,loop->size()-1);
                move_loop = new PProtein(frontLoop,cfFrontEndLength,frontLoop->size()-1);
                endsDistanceThreshold = computeMaxLength(middleSize) * MAX_LENGTH_DISCOUNT_RATIO;
                sr = new SpaceRelationship(
                                original_protein->getResidue(loopEid-cfBackEndLength),
                                original_protein->getResidue(loopEid-cfBackEndLength+1));

        for (int i=0; i<cfFrontEndLength; ++i) {
                front_aa_names.push_back(original_protein->getResidue(loopSid+i)->getResourceName());
        }
        for (int i=0; i<cfBackEndLength; ++i) {
                back_aa_names.push_back(original_protein->getResidue(loopSid+frontLoop->size()+i)->getResourceName());
        }

        }
        else {
                move_loop = loop;
                PResidue *res = protein->getResidue(loopEid);
                endPriorG = res->getAtomPosition("CA");
                endG = res->getAtomPosition("C");
                endNextG = res->getAtomPosition("O");
        }

        for (int i=0; i<move_loop->size(); ++i) {
                move_aa_names.push_back(original_protein->getResidue(loopSid+cfFrontEndLength+i)->getResourceName());
        }

        while (result.size() < num_wanted) {
                if (split) {
                                // Randomize the loop
                                for (int i=0; i<2*frontLoop->size(); ++i)
                                        frontLoop->RotateBackbone(i,forward,rand()%360);
                                for (int i=0; i<2*backLoop->size(); ++i)
                                        backLoop->RotateBackbone(i,backward,rand()%360);

                                // Generate new initial forward and backward ends
				loop->inactivateResidue(0,loop->size()-1);
				bool gotends=false;
				while (!gotends) {
					while (true) {
						angles = generateForwardOpenLoop(frontLoop,cfFrontEndLength,distri_map,front_aa_names);
						if (!angles->empty()) {
							angles->clear();
							delete angles;
							break;
						}
					}
					for (int i=0; i<10; ++i) {
						angles = generateBackwardOpenLoop(backLoop,distri_map,back_aa_names);
						if (!angles->empty()) {
							gotends = true;
							angles->clear();
							delete angles;
							break;
						}
					}
				}
			
                                frontLastAtomPos = frontLoop->getResidue(cfFrontEndLength-1)->getAtomPosition("C");
                                backFirstAtomPos = backLoop->getResidue(0)->getAtomPosition("N");
                        okEndPair = true;
                        if (frontLastAtomPos.distance(backFirstAtomPos) > endsDistanceThreshold) {
                                okEndPair = false;
                        }
                         // Check if these two ends collide
                        if (okEndPair) {
                                move_loop->inactivateResidue(0,move_loop->size()-1);
                                okEndPair = !(backLoop->InAnyCollision());
                                move_loop->activateResidue(0,move_loop->size()-1);
                        }

                        if (okEndPair) {
                                PResidue *backHeadResInProtein = protein->getResidue(loopEid-cfBackEndLength+1);
                                PSampMethods::computeGoal(backHeadResInProtein,sr,&endPriorG,&endG,&endNextG);
                        }
                } // end if

                for (int tryNum=0; okEndPair && tryNum<MAX_TRIAL_PER_PAIR; ++tryNum) {
                        // Randomize backbone DOFs
                        vector<ChainMove> cms;
                        ChainMove cm;
                        for (int j=0; j<2*move_loop->size(); ++j) {
                                cm.blockType = PID::BACKBONE;
                                cm.dir = forward;
                                cm.DOF_index = j;
//                                cm.degrees = rand()%360;
                                        aa_name = move_aa_names[j/2];
                                        map_iter = distri_map.find(aa_name);
				        if (j%2==0) { // Phi 
                                                phi = sampleAngle((map_iter->second).getPhiDistribution());
                                                curPhi = move_loop->getResidue(j/2)->GetPhi();
                                                cm.degrees = curPhi - phi;
                                        }
                                        else { // Psi
                                                vector<double> psiDistri = (map_iter->second).getPsiDistribution(phi);
                                                psi = sampleAngle(psiDistri);
                                                curPsi = move_loop->getResidue(j/2)->GetPsi();
                                                cm.degrees = curPsi - psi;
                                        }


                                cms.push_back(cm);
                        }
                        move_loop->MultiRotate(cms);

                        // Close the loop
                        bool no_sol = PSampMethods::IKClose(move_loop,endPriorG,endG,endNextG);
                        // Check collision and avoid collision if less enough
                        if (no_sol) {
                                continue;
                        }
                        if (!move_loop->InAnyCollision()) {
                                result.push_back(loop->Clone());
                                break;
                        }               
                } // end for
        } // end while          
                                
        if (split)                      
                delete sr;      
        delete protein;         
                                
        return result;  
}

// For AddSidechain
    class equalString
    {
        public:
            size_t operator()(string const &s1, string const &s2) const
            {
                return (s1 == s2);
            }
    };

    class hashString
    {
        public:
            size_t operator()(string const &str) const
            {
                __gnu_cxx::hash<char const *>
                    h;

                return (h(str.c_str()));
            }
    };

void PSampMethods::AddSidechain (string protein_input, int addStart, int addEnd, string scwrl3_path, string protein_output) {
	// Create the sequence file
	__gnu_cxx::hash_map<string,char,hashString,equalString> ResNameMap;
	ResNameMap["ALA"] = 'a';
	ResNameMap["ARG"] = 'r';
	ResNameMap["ASN"] = 'n';
	ResNameMap["ASP"] = 'd';
	ResNameMap["CYS"] = 'c';
	ResNameMap["GLU"] = 'e';
	ResNameMap["GLN"] = 'q';
	ResNameMap["GLY"] = 'g';
	ResNameMap["HIS"] = 'h';
	ResNameMap["ILE"] = 'i';
	ResNameMap["LEU"] = 'l';
	ResNameMap["LYS"] = 'k';
	ResNameMap["MET"] = 'm';
	ResNameMap["PHE"] = 'f';
	ResNameMap["PRO"] = 'p';
	ResNameMap["SER"] = 's';
	ResNameMap["THR"] = 't';
	ResNameMap["TRP"] = 'w';
	ResNameMap["TYR"] = 'y';
	ResNameMap["VAL"] = 'v';
	ofstream seq_out("seq.out",ios::out);
	PProtein *p = PDBIO::readFromFile(protein_input);
	for (int i=0; i<p->size(); ++i) {
		PResidue *res = p->getResidue(i);
		string res_string_name = res->getName();
		char res_c_name = ResNameMap[res_string_name];
		int res_c_name_int = int(res_c_name); 
		if (res_c_name_int<65 || (res_c_name_int>90&&res_c_name_int<97) || res_c_name_int>122)
			res_c_name = 'o'; 
		int pdbid = res->getPdbId();
		if (pdbid>=addStart && pdbid<=addEnd) {
			res_c_name = toupper(res_c_name);
		}
		seq_out << res_c_name;
	}
	delete p;
	seq_out << endl;
	seq_out.close();

	// Run SCWRL3	
	string command = scwrl3_path+"/scwrl3 -i "+protein_input+" -o temp.out -s seq.out > scwrl3.log";
	system(command.c_str());
	
	// Re-format the output file so that it is in proper PDB format
        ifstream fin("temp.out");
        ofstream fout(protein_output.c_str());
        string s, s1, s2=" 1.00                  ";
        while (!fin.eof()) {
                getline(fin,s);
                if (s.size()<30)
                        break;
                s.erase(30,1);
                s1 = s.substr(0,55);
                fout << s1 << s2 << endl;
        }
        fin.close();
        fout.close();

	system("rm seq.out");
	system("rm temp.out");
	system("rm scwrl3.log");
}
		
		

vector<PProtein*> PSampMethods::SeedSampleBackboneWithSidechainLoopOnly (PProtein* original_protein, int loopSid, int loopEid, string scwrl3_path, int num_wanted) {

        static int MIN_MOVE_LOOP_SIZE=4;
        static int MAX_TRIAL_PER_PAIR=50;
        static int C_THRESHOLD=0, MAX_ENDS_NUM=max(num_wanted*2,20);
        static double MIDDLE_SIZE_RATIO=0.5;
        static double MAX_LENGTH_DISCOUNT_RATIO = 1;
        int cfFrontEndLength=0, cfBackEndLength=0;
        PProtein *frontLoop, *backLoop, *loop, *move_loop, *protein;
        Vector3 endG, endPriorG, endNextG, frontLastAtomPos, backFirstAtomPos;
        SpaceRelationship *sr;
        bool workFrontEnd=false, okEndPair=true;
        double endsDistanceThreshold=0;
        vector<PProtein*> result;
        vector<Real> *angles;

        protein = PTools::CreateSlimProtein(original_protein,loopSid,loopEid);

        int loopSize = loopEid-loopSid+1;
        int middleSize = (int)(floor(MIDDLE_SIZE_RATIO*loopSize));
        bool split = middleSize>=MIN_MOVE_LOOP_SIZE;

        loop = new PProtein(protein,loopSid,loopEid);
        if (split) {
                cfFrontEndLength = (loopSize-middleSize)/2;
                cfBackEndLength = loopSize-middleSize-cfFrontEndLength;
                frontLoop = new PProtein(loop,0,loop->size()-cfBackEndLength-1);
                backLoop = new PProtein(loop,loop->size()-cfBackEndLength,loop->size()-1);
                move_loop = new PProtein(frontLoop,cfFrontEndLength,frontLoop->size()-1);
                endsDistanceThreshold = computeMaxLength(middleSize) * MAX_LENGTH_DISCOUNT_RATIO;
                sr = new SpaceRelationship(
                                original_protein->getResidue(loopEid-cfBackEndLength),
                                original_protein->getResidue(loopEid-cfBackEndLength+1));
        }
        else {
                move_loop = loop;
                PResidue *res = protein->getResidue(loopEid);
                endPriorG = res->getAtomPosition("CA");
                endG = res->getAtomPosition("C");
                endNextG = res->getAtomPosition("O");
        }

        while (result.size() < num_wanted) {
                if (split) {
                                // Randomize the loop
                                for (int i=0; i<2*frontLoop->size(); ++i)
                                        frontLoop->RotateBackbone(i,forward,rand()%360);
                                for (int i=0; i<2*backLoop->size(); ++i)
                                        backLoop->RotateBackbone(i,backward,rand()%360);

                                // Generate new initial forward and backward ends
                                loop->inactivateResidue(0,loop->size()-1);
                                bool gotends=false;                                                               
                                while (!gotends) {
                                        while (true) {                                                            
                                                angles = generateForwardOpenLoop(frontLoop,cfFrontEndLength);
                                                if (!angles->empty()) {                                           
							angles->clear();
                                                        delete angles;                                            
                                                        break;                                                    
                                                }                                                                 
                                        }                                                                         
                                        for (int i=0; i<10; ++i) {                                                
                                                angles = generateBackwardOpenLoop(backLoop);
                                                if (!angles->empty()) {                                           
                                                        gotends = true;                                          
							angles->clear();
							delete angles; 
                                                        break;                                                    
                                                }                                                                 
                                        }                                                                         
                                }

                                frontLastAtomPos = frontLoop->getResidue(cfFrontEndLength-1)->getAtomPosition("C");
                                backFirstAtomPos = backLoop->getResidue(0)->getAtomPosition("N");
                        okEndPair = true;
                        if (frontLastAtomPos.distance(backFirstAtomPos) > endsDistanceThreshold) {
                                okEndPair = false;
                        }
                         // Check if these two ends collide
                        if (okEndPair) {
                                move_loop->inactivateResidue(0,move_loop->size()-1);
                                okEndPair = !(backLoop->InAnyCollision());
                                move_loop->activateResidue(0,move_loop->size()-1);
                        }

                        if (okEndPair) {
                                PResidue *backHeadResInProtein = protein->getResidue(loopEid-cfBackEndLength+1);
                                PSampMethods::computeGoal(backHeadResInProtein,sr,&endPriorG,&endG,&endNextG);
                        }
                } // end if

                for (int tryNum=0; okEndPair && tryNum<MAX_TRIAL_PER_PAIR; ++tryNum) {
                        // Randomize backbone DOFs
                        vector<ChainMove> cms;
                        ChainMove cm;
                        for (int j=0; j<2*move_loop->size(); ++j) {
                                cm.blockType = PID::BACKBONE;
                                cm.dir = forward;
                                cm.DOF_index = j;
                                cm.degrees = rand()%360;
                                cms.push_back(cm);
                        }
                        move_loop->MultiRotate(cms);

                        // Close the loop
                        bool no_sol = PSampMethods::IKClose(move_loop,endPriorG,endG,endNextG);
                        // Check collision and avoid collision if less enough
                        if (no_sol) {
                                continue;
                        }       
                        if (!move_loop->InAnyCollision()) {
                        	// Add sidechain
                        	PDBIO::writeToFile(protein,"temp.pdb");
                        	int pdbSid = protein->getResidue(loopSid)->getPdbId();
                        	int pdbEid = protein->getResidue(loopEid)->getPdbId();
                        	AddSidechain("temp.pdb",pdbSid,pdbEid,scwrl3_path,"tempSC.pdb");
                        	// Check for collision
                        	FILE *file = fopen("tempSC.pdb","r");
                        	if (file!=NULL) {
                        		fclose(file);
                        		PProtein *seed_with_sc = PDBIO::readFromFile("tempSC.pdb");
                        		PProtein *loop_with_sc = new PProtein(seed_with_sc,loopSid,loopEid);
                        		if (!loop_with_sc->InAnyCollision()) {
                        			result.push_back(loop_with_sc->Clone());
						seed_with_sc->Obliterate();         			
						system("rm tempSC.pdb");
						system("rm temp.pdb");
                        			break;
                        		}
                        		seed_with_sc->Obliterate();
                        		system("rm tempSC.pdb");
                        	}
                        	system("rm temp.pdb");
                        }
                } // end for
        } // end while  
                        
        if (split)      
                delete sr;      
        delete protein;         

        return result;
}

vector<PProtein*> PSampMethods::SeedSampleBackboneWithSidechainLoopOnly (PProtein* original_protein, int loopSid, int loopEid, map<string,PPhiPsiDistribution> &distri_map, string scwrl3_path, int num_wanted){

        if (distri_map.size()==20) {
                if (distri_map.count("ALA")==0 || distri_map.count("ARG")==0 || distri_map.count("ASN")==0
                || distri_map.count("ASP")==0 || distri_map.count("CYS")==0 || distri_map.count("GLN")==0
                || distri_map.count("GLU")==0 || distri_map.count("GLY")==0 || distri_map.count("HIS")==0
                || distri_map.count("ILE")==0 || distri_map.count("LEU")==0 || distri_map.count("LYS")==0
                || distri_map.count("MET")==0 || distri_map.count("PHE")==0 || distri_map.count("PRO")==0
                || distri_map.count("SER")==0 || distri_map.count("THR")==0 || distri_map.count("TRP")==0
                || distri_map.count("TYR")==0 || distri_map.count("VAL")==0 ) {
                        cerr << "Wrong distribution names. Please read the documentation carefull." << endl;
                        exit(1);
                }
        }

        static int MIN_MOVE_LOOP_SIZE=4;
        static int MAX_TRIAL_PER_PAIR=50;
        static int C_THRESHOLD=0, MAX_ENDS_NUM=max(num_wanted*2,20);
        static double MIDDLE_SIZE_RATIO=0.5;
        static double MAX_LENGTH_DISCOUNT_RATIO = 1;
        int cfFrontEndLength=0, cfBackEndLength=0;
        PProtein *frontLoop, *backLoop, *loop, *move_loop, *protein;
        Vector3 endG, endPriorG, endNextG, frontLastAtomPos, backFirstAtomPos;
        SpaceRelationship *sr;
        bool workFrontEnd=false, okEndPair=true;
        double endsDistanceThreshold=0;
        vector<PProtein*> result;
        vector<Real> *angles;

        protein = PTools::CreateSlimProtein(original_protein,loopSid,loopEid);
        
        int loopSize = loopEid-loopSid+1;
        int middleSize = (int)(floor(MIDDLE_SIZE_RATIO*loopSize));
//        bool split = middleSize>=(MIN_MOVE_LOOP_SIZE+1);
        bool split = middleSize>=MIN_MOVE_LOOP_SIZE;

	vector<string> front_aa_names, back_aa_names, move_aa_names;
	map<string,PPhiPsiDistribution>::iterator map_iter;
	string aa_name;
	double phi, psi, curPhi, curPsi;
        
        loop = new PProtein(protein,loopSid,loopEid);
        if (split) {
                cfFrontEndLength = (loopSize-middleSize)/2;
                cfBackEndLength = loopSize-middleSize-cfFrontEndLength;
                frontLoop = new PProtein(loop,0,loop->size()-cfBackEndLength-1);
                backLoop = new PProtein(loop,loop->size()-cfBackEndLength,loop->size()-1);
                move_loop = new PProtein(frontLoop,cfFrontEndLength,frontLoop->size()-1);
                endsDistanceThreshold = computeMaxLength(middleSize) * MAX_LENGTH_DISCOUNT_RATIO;
                sr = new SpaceRelationship(
                                original_protein->getResidue(loopEid-cfBackEndLength),
                                original_protein->getResidue(loopEid-cfBackEndLength+1));

        for (int i=0; i<cfFrontEndLength; ++i) {
                front_aa_names.push_back(original_protein->getResidue(loopSid+i)->getResourceName());
        }
        for (int i=0; i<cfBackEndLength; ++i) {
                back_aa_names.push_back(original_protein->getResidue(loopSid+frontLoop->size()+i)->getResourceName());
        }

        }
        else {
                move_loop = loop;
                PResidue *res = protein->getResidue(loopEid);
                endPriorG = res->getAtomPosition("CA");
                endG = res->getAtomPosition("C");
                endNextG = res->getAtomPosition("O");
        }

        for (int i=0; i<move_loop->size(); ++i) {
                move_aa_names.push_back(original_protein->getResidue(loopSid+cfFrontEndLength+i)->getResourceName());
        }

        while (result.size() < num_wanted) {
                if (split) {
                                // Randomize the loop
                                for (int i=0; i<2*frontLoop->size(); ++i)
                                        frontLoop->RotateBackbone(i,forward,rand()%360);
                                for (int i=0; i<2*backLoop->size(); ++i)
                                        backLoop->RotateBackbone(i,backward,rand()%360);

                                // Generate new initial forward and backward ends
				loop->inactivateResidue(0,loop->size()-1);
				bool gotends=false;
				while (!gotends) {
					while (true) {
						angles = generateForwardOpenLoop(frontLoop,cfFrontEndLength,distri_map,front_aa_names);
						if (!angles->empty()) {
							angles->clear();
							delete angles;
							break;
						}
					}
					for (int i=0; i<10; ++i) {
						angles = generateBackwardOpenLoop(backLoop,distri_map,back_aa_names);
						if (!angles->empty()) {
							gotends = true;
							angles->clear();
							delete angles;
							break;
						}
					}
				}

			
                                frontLastAtomPos = frontLoop->getResidue(cfFrontEndLength-1)->getAtomPosition("C");
                                backFirstAtomPos = backLoop->getResidue(0)->getAtomPosition("N");
                        okEndPair = true;
                        if (frontLastAtomPos.distance(backFirstAtomPos) > endsDistanceThreshold) {
                                okEndPair = false;
                        }
                         // Check if these two ends collide
                        if (okEndPair) {
                                move_loop->inactivateResidue(0,move_loop->size()-1);
                                okEndPair = !(backLoop->InAnyCollision());
                                move_loop->activateResidue(0,move_loop->size()-1);
                        }

                        if (okEndPair) {
                                PResidue *backHeadResInProtein = protein->getResidue(loopEid-cfBackEndLength+1);
                                PSampMethods::computeGoal(backHeadResInProtein,sr,&endPriorG,&endG,&endNextG);
                        }

                } // end if

                for (int tryNum=0; okEndPair && tryNum<MAX_TRIAL_PER_PAIR; ++tryNum) {
                        // Randomize backbone DOFs
                        vector<ChainMove> cms;
                        ChainMove cm;
                        for (int j=0; j<2*move_loop->size(); ++j) {
                                cm.blockType = PID::BACKBONE;
                                cm.dir = forward;
                                cm.DOF_index = j;
//                                cm.degrees = rand()%360;
                                        aa_name = move_aa_names[j/2];
                                        map_iter = distri_map.find(aa_name);
				        if (j%2==0) { // Phi 
                                                phi = sampleAngle((map_iter->second).getPhiDistribution());
                                                curPhi = move_loop->getResidue(j/2)->GetPhi();
                                                cm.degrees = curPhi - phi;
                                        }
                                        else { // Psi
                                                vector<double> psiDistri = (map_iter->second).getPsiDistribution(phi);
                                                psi = sampleAngle(psiDistri);
                                                curPsi = move_loop->getResidue(j/2)->GetPsi();
                                                cm.degrees = curPsi - psi;
                                        }


                                cms.push_back(cm);
                        }
                        move_loop->MultiRotate(cms);

                        // Close the loop
                        bool no_sol = PSampMethods::IKClose(move_loop,endPriorG,endG,endNextG);
                        // Check collision and avoid collision if less enough
                        if (no_sol) {
                                continue;
                        }
                        if (!move_loop->InAnyCollision()) {
                        	// Add sidechain
                        	PDBIO::writeToFile(protein,"temp.pdb");
                        	int pdbSid = protein->getResidue(loopSid)->getPdbId();
                        	int pdbEid = protein->getResidue(loopEid)->getPdbId();
                        	AddSidechain("temp.pdb",pdbSid,pdbEid,scwrl3_path,"tempSC.pdb");
                        	// Check for collision
                        	FILE *file = fopen("tempSC.pdb","r");
                        	if (file!=NULL) {
                        		fclose(file);
                        		PProtein *seed_with_sc = PDBIO::readFromFile("tempSC.pdb");
                        		PProtein *loop_with_sc = new PProtein(seed_with_sc,loopSid,loopEid);
                        		if (!loop_with_sc->InAnyCollision()) {
                        			result.push_back(loop_with_sc->Clone());
						seed_with_sc->Obliterate();
						system("rm tempSC.pdb");
						system("rm temp.pdb");
                        			break;
                        		} 
                        		seed_with_sc->Obliterate();
                        		system("rm tempSC.pdb");
                        	}
                        	system("rm temp.pdb");
                        }
                } // end for
        } // end while          
                                
        if (split)                      
                delete sr;      
        delete protein;         
                                
        return result;  
}

vector<PProtein*> PSampMethods::SeedSampleBackboneWithSidechain (PProtein *protein, int loopSid, int loopEid, string scwrl3_path, int num_wanted) {
        vector<PProtein*> loops = SeedSampleBackboneWithSidechainLoopOnly(protein,loopSid,loopEid,scwrl3_path,num_wanted);
        vector<PProtein*> result;
        int loopsize = loopEid-loopSid+1;
        for (int i=0; i<loops.size(); ++i) {
		PProtein *p = protein->Clone();
                for (int j=0; j<loopsize; ++j) {
                        p->getResidue(loopSid+j)->setPositions(loops[i]->getResidue(j)->getSpec().Atom_Positions);
                }
                result.push_back(p);
                loops[i]->Obliterate();
        }
        return result;
}

vector<PProtein*> PSampMethods::SeedSampleBackboneWithSidechain (PProtein *protein, int loopSid, int loopEid, map<string,PPhiPsiDistribution> &distri_map, string scwrl3_path, int num_wanted) {
        vector<PProtein*> loops = SeedSampleBackboneWithSidechainLoopOnly(protein,loopSid,loopEid,distri_map,scwrl3_path,num_wanted);      
        vector<PProtein*> result;                                                                                 
        int loopsize = loopEid-loopSid+1;                                                                         
        for (int i=0; i<loops.size(); ++i) {   
		PProtein *p = protein->Clone();                                                                   
                for (int j=0; j<loopsize; ++j) {                                                                  
                        p->getResidue(loopSid+j)->setPositions(loops[i]->getResidue(j)->getSpec().Atom_Positions);
                }                                                                                                 
                result.push_back(p);                                                                              
                loops[i]->Obliterate();                                                                           
        }                                                                                                         
        return result;                                                                                            
}

PProtein* PSampMethods::FillMissingLoop (PProtein *original_p, int start_pdb_id, int end_pdb_id, vector<string> loop_seq) {
        // Missing loop cannot be the beginning or the end of the protein
        if (original_p->getResidue(original_p->size()-1)->getPdbId()<=end_pdb_id ||
            original_p->getResidue(0)->getPdbId()>=start_pdb_id)
                PUtilities::AbortProgram("Missing loop cannot be the beginning or the end of the protein");

	// Build a temp loop with wrong O-C-N angle due to the mistake in PProtein::AddResidue 
	PProtein *temp = new PProtein();
	PResidue *res;
	PResidueSpec spec;
	// Add the residue before missing loop
	int headEndRid = original_p->pdbIndexToLocalIndex(start_pdb_id-1);
        int tailStartRid = original_p->pdbIndexToLocalIndex(end_pdb_id+1);
	res = original_p->getResidue(headEndRid);
	spec = res->getSpec();
	temp->AddResidue(res->getResourceName(),spec);
	// Model the missing loop (with wrong O-C-N angle)
	for (int i=0; i<loop_seq.size(); ++i) {
		temp->AddResidue(loop_seq.at(i));
		temp->getResidue(temp->size()-1)->setPdbId(start_pdb_id+i);
	}
        // Add one more residue
        res = original_p->getResidue(tailStartRid);
        temp->AddResidue(res->getResourceName());
        temp->getResidue(temp->size()-1)->setPdbId(end_pdb_id+1);
	temp->finalize();

	// Build the head part in output protein new_p
	PProtein *new_p = new PProtein();
        for (int i=0; i<=headEndRid; ++i) {
                res = original_p->getResidue(i);
                spec = res->getSpec();
                new_p->AddResidue(res->getResourceName(),spec);
        }	

	// Rotate the Psi angles to make angle O-C-N in protein temp correct and copy to new_p
	Vector3 c,o,n,c_n,c_o;
	Real cos_o_c_n;
	double angle_o_c_n, diff, best_diff;
	double best_angle;
	for (int i=0; i<temp->size()-1; ++i) {
		// res is the residue that has correct C and O positions
		res = temp->getResidue(i);
		c = res->getAtomPosition("C");
		o = res->getAtomPosition("O");
		// compute the current angle
		PResidue *res2 = temp->getResidue(i+1);
		n = res2->getAtomPosition("N");
		c_n = n-c;
		c_o = o-c;
		cos_o_c_n = c_o.dot(c_n) / (c_o.norm()*c_n.norm());
		angle_o_c_n = acos(cos_o_c_n)*180/PI;
		best_diff = PMath::abs(double(angle_o_c_n-ANGLE_O_C_N));
		best_angle = angle_o_c_n;
		int best_rotation_angle = 0;
		for (int j=0; j<360; ++j) {
			res->getDOF(PID::BACKBONE,PID::C_ALPHA,PID::C)->Rotate(forward,1);
			n = res2->getAtomPosition("N");
			c_n = n-c;
			c_o = o-c;
			cos_o_c_n = c_o.dot(c_n) / (c_o.norm()*c_n.norm());
			angle_o_c_n = acos(cos_o_c_n)*180/PI;
			diff = PMath::abs(angle_o_c_n-ANGLE_O_C_N);
			if (diff < best_diff) {
				spec = res2->getSpec();
				best_diff = diff;
				best_rotation_angle = j+1;
				best_angle = angle_o_c_n;
			}
		}
		new_p->AddResidue(res2->getResourceName(),spec);
		res->getDOF(PID::BACKBONE,PID::C_ALPHA,PID::C)->Rotate(forward,best_rotation_angle);
	}

        // Copy the tail part except the first residue after the missing loop
        for (int i=tailStartRid+1; i<original_p->size(); ++i) {
                res = original_p->getResidue(i);
                spec = res->getSpec();
                new_p->AddResidue(res->getResourceName(),spec);
        }
        new_p->finalize();

        PProteinResidue *anchor = original_p->getResidue(tailStartRid);
        // Get the position of goal end frame
        Vector3 endPriorG = anchor->getAtomPosition("CA");
        Vector3 endG = anchor->getAtomPosition("C");
        Vector3 endNextG = anchor->getAtomPosition("O");
        // Get the positions of N and Cb of the anchor
        Vector3 anchor_N = anchor->getAtomPosition("N");
        Vector3 anchor_Cb = anchor->getAtomPosition("CB");

        // Randomize and close the loop
        PProtein *loop = new PProtein(new_p,headEndRid+1,headEndRid+1+end_pdb_id+1-start_pdb_id);
PResidue *lastRes = loop->getResidue(loop->size()-1);
        bool not_found = true;
        while (not_found) {
                // Randomize the loop	
		for (int i=0; i<2*loop->size(); ++i) 
			loop->RotateBackbone(i,forward,rand()%360);
                // Find the way to close the loop while keep N and Cb as close to target as possible
                IKSolutions min_sol;
                int res_indices[3];
                res_indices[2] = loop->size()-1;
                for (int i=0; not_found && i<loop->size()-2; ++i) {
                        res_indices[0] = i;
                        for (int j=i+1; not_found && j<loop->size()-1; ++j) {
                                res_indices[1] = j;
                                IKSolutions sols = PExactIKSolver::FindSolutions(loop,res_indices,&endPriorG,&endG,&endNextG);
                                for (int k=0; not_found && k<sols.size(); ++k) {
                                        loop->MultiRotate(sols[k]);
                                        Vector3 N_pos = loop->getResidue(loop->size()-1)->getAtomPosition("N");
                                        Vector3 Cb_pos = loop->getResidue(loop->size()-1)->getAtomPosition("CB");
                                        double distance = pow(N_pos.distance(anchor_N),2) + pow(Cb_pos.distance(anchor_Cb),2);
                                        if (distance<=0.1) {
                                                not_found = false;
                                                break;
                                        }
                                        loop->AntiMultiRotate(sols[k]);
                                }
                        }
                }
        }
	
	delete temp;
        return new_p;
}

