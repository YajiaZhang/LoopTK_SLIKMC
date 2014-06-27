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

void PProteinResidue::SaveSideChain(){
  originalAngle.clear();
  //see how many chi angles in this given residue.
  unsigned chiMax = PResources::numChiIndices(this->getName());
  //store the name of atoms defining the chi angle.
  vector<string> rotList;

  //store side chain angles
  for(unsigned i = 1; i <= chiMax; i++){
    rotList = PResources::GetChiIndex(this->getName(), i);
    originalAngle.push_back(-PMath::AngleBetweenPlanes(getAtomPosition(rotList[0]),
                                                       getAtomPosition(rotList[1]),
                                                       getAtomPosition(rotList[2]),
                                                       getAtomPosition(rotList[3])));
  }
}

void PProteinResidue::ResetSideChain(){
  if (originalAngle.size()!=0){
    TryApplyRotamer(originalAngle);
    currRot=0;
  }
}

bool PProteinResidue::ApplyRotamer(int index){
  //checks if rotamer can be applied to a given residue.
  string resName = this->getName();
  if (!PResources::ContainsRotamer(resName)) {
    cerr<<resName<<" has no rotamer."<<endl;
    return false;
  }
  
  //store all rotamers for given residue.
  vector<vector<Real> > rotamerAngles = PResources::GetRotamer(resName);
  
  //check index range.
  if (index>rotamerAngles.size()){
    cerr<<"Rotamer index out of range."<<endl;
    return false;
  }
  
  TryApplyRotamer(rotamerAngles[index]);
  return true;
}

//void PProteinResidue::applyRotamer(const int& index, const vector<double>& rotamer) {
//
//	cout << "calling PProteinResidue::applyRotamer" << endl;
//	this->angles_sidechain.resize(rotamer.size());
//	for (int i = 0; i < rotamer.size(); i++) {
//		this->angles_sidechain[i] = rotamer[i];
//	}
//	this->type_sidechain = index;
//	//see how many chi angles in this give residue.
//	unsigned chiMax = PResources::numChiIndices(this->getName());
//
////	  assert( chiMax == rotamer.size());
//	if (!(chiMax == rotamer.size())) {
//		cout << "chiMax:" << chiMax << endl;
//		cout << "rotamer.size:" << rotamer.size() << endl;
//		abort();
//	}
//
//	//store the name of atoms that defines the chi angle.
//	vector<string> rotList;
//	//current dihedral angl.
//	Real curDihedral;
//	//amount to rotate.
//	Real toRotate;
//
//	//apply rotamer angles
//	for (unsigned i = 1; i <= chiMax; i++) {
//		//get the atoms used to define dihedral angles
//		rotList = PResources::GetChiIndex(this->getName(), i);
//
////		for( int j = 0; j < rotList.size(); j++) {
////			cout << rotList[j] << "\t";
////		}
////		cout << endl;
//
//		//calculate current dihedral angle (chi angle) and calculate the delta needed to achive the goal.
//		curDihedral = PMath::AngleBetweenPlanes(
//				this->getAtomPosition(rotList[0]),
//				this->getAtomPosition(rotList[1]),
//				this->getAtomPosition(rotList[2]),
//				this->getAtomPosition(rotList[3]));
//		/************************************/
//		//convert to positive
////		if (curDihedral < 0) {
////			curDihedral = 360 + curDihedral;
////			//apply rotation to chi angle..
////			this->getDOF("sidechain", rotList[1], rotList[2])->Rotate_side(forward, (rotamer[i - 1] - curDihedral), NULL);
////		}
////
////		//convert to positive
////		if (-rotamer[i - 1] < 0) {
////			toRotate = ((360 - rotamer[i - 1]) - curDihedral);
////		}
////		else {
////			toRotate = ((-rotamer[i - 1]) - curDihedral);
////		}
////		cout << "  chi " << i << " toRotate:" << toRotate << endl;
////		this->getDOF("sidechain", rotList[1], rotList[2])->Rotate_side(forward, (-toRotate), NULL);
////
////		if (toRotate > 0) {
////			//apply rotation to chi angle..
////			this->getDOF("sidechain", rotList[1], rotList[2])->Rotate_side(backward, (-toRotate), NULL);
////		}
////		else {
////			//apply rotation to chi angle..
////			this->getDOF("sidechain", rotList[1], rotList[2])->Rotate_side(backward, (-toRotate), NULL);
////		}
//		/************************************/
//
////		//convert to positive
////		if (curDihedral < 0) {
////			curDihedral = 360 + curDihedral;
////			//apply rotation to chi angle..
////			this->getDOF("sidechain", rotList[1], rotList[2])->Rotate(forward, (rotamer[i - 1] - curDihedral), this->getChain());
////		}
////
////		//convert to positive
////		if (-rotamer[i - 1] < 0) {
////			toRotate = ((360 - rotamer[i - 1]) - curDihedral);
////		}
////		else {
////			toRotate = ((-rotamer[i - 1]) - curDihedral);
////		}
////
////		cout << "  chi " << i << " toRotate:" << toRotate << endl;
////		this->getDOF("sidechain", rotList[1], rotList[2])->Rotate(forward, (-toRotate), this->getChain());
////
////		if (toRotate > 0) {
////			//apply rotation to chi angle..
////			this->getDOF("sidechain", rotList[1], rotList[2])->Rotate(backward, (-toRotate), this->getChain());
////		}
////		else {
////			//apply rotation to chi angle..
////			this->getDOF("sidechain", rotList[1], rotList[2])->Rotate(backward, (-toRotate), this->getChain());
////		}
//		/************************************/
//
//	    this->getDOF("sidechain",rotList[1],rotList[2])->Rotate_noGridUpdate(forward,(curDihedral - rotamer[i-1]), this->getChain());
//	}
//	return;
//
//}

bool PProteinResidue::TryApplyRotamer(const vector<Real> &rotamer)
{
  //see how many chi angles in this give residue.
  unsigned chiMax = PResources::numChiIndices(this->getName());
  //store the name of atoms that defines the chi angle.
  vector<string> rotList;
  //current dihedral angl.
  Real curDihedral;
  //amount to rotate.
  Real toRotate;

  //apply rotamer angles
  for(unsigned i = 1; i<=chiMax; i++){
    //get the atoms used to define dihedral angles
    rotList = PResources::GetChiIndex(this->getName(), i);
 
    //calculate current dihedral angle (chi angle) and calculate the delta needed to achive the goal.
    curDihedral = PMath::AngleBetweenPlanes(this->getAtomPosition(rotList[0]), this->getAtomPosition(rotList[1]), this->getAtomPosition(rotList[2]), this->getAtomPosition(rotList[3]));

  //convert to positive
  if (curDihedral<0){
    curDihedral=360+curDihedral;
    //apply rotation to chi angle..
    this->getDOF("sidechain",rotList[1],rotList[2])->Rotate(forward,(rotamer[i-1]-curDihedral),NULL);
  }

  //convert to positive
  if (-rotamer[i-1]<0){
    toRotate = ((360-rotamer[i-1])-curDihedral);
  }
  else{
    toRotate = ((-rotamer[i-1])-curDihedral);
  }
  
  this->getDOF("sidechain",rotList[1],rotList[2])->Rotate(forward,(-toRotate),NULL);
    
  if (toRotate > 0){
      //apply rotation to chi angle..
      this->getDOF("sidechain",rotList[1],rotList[2])->Rotate(backward,(-toRotate),NULL);
    }
    else{
      //apply rotation to chi angle..
      this->getDOF("sidechain",rotList[1],rotList[2])->Rotate(backward,(-toRotate),NULL);
    }
  }
}

bool PProteinResidue::ApplyRotamer(){
  //checks if rotamer can be applied to a given residue.
  string resName = this->getName();
  if (!PResources::ContainsRotamer(resName)) {
    cerr<<resName<<" has no rotamer"<<endl;
    return true;
  }

  //PRO will be implemented.
  if (this->getName() == PID::PRO){
    return true;
  }

  //cache original position.
  SaveSideChain();
  //store all rotamers for given residue.
  vector<vector<Real> > rotamerAngles = PResources::GetRotamer(resName);
  
  //cycle through romaters and see if collision free; keeptrack of current position.
  for(unsigned i = 0; i < rotamerAngles.size(); i++) {
    TryApplyRotamer(rotamerAngles[(i+currRot)%rotamerAngles.size()]);
    if (!this->InAnyCollision()) {
      currRot=((i+currRot)+1)%rotamerAngles.size();
      return true;
    }
  }
  
  ResetSideChain();
  cerr<<"No rotamers are collision-free. Reset to original position."<<endl;
  return false;

}

// C(t-1) -> N(t) -> Ca(t) -> C(t)
//void PProteinResidue::SetPhi(Real phi) {
 // this->getDOF(PID::BACKBONE,PID::N,PID::C_ALPHA)->Rotate(forward,(phi-GetPhi()));
//}

Real PProteinResidue::GetPhi(){
  return PMath::AngleBetweenPlanes(this->PreviousResidue()->getAtomPosition(PID::C),
                                   this->getAtomPosition(PID::N),
                                   this->getAtomPosition(PID::C_ALPHA),
                                   this->getAtomPosition(PID::C));
}

// N(t) -> Ca(t) -> C(t) -> N(t+1)
//void PProteinResidue::SetPsi(Real psi){
  //this->getDOF(PID::BACKBONE,PID::C_ALPHA,PID::C)->Rotate(forward,(psi-GetPsi()));
//}

Real PProteinResidue::GetPsi(){
  return PMath::AngleBetweenPlanes(this->getAtomPosition(PID::N),
                                   this->getAtomPosition(PID::C_ALPHA),
                                   this->getAtomPosition(PID::C),
                                   this->NextResidue()->getAtomPosition(PID::N));
}

void PProteinResidue::SetChi(Real chi, int index){
  vector<string> rotList = PResources::GetChiIndex(this->getName(), index);
  //this->getDOF("sidechain",rotList[1],rotList[2])->Rotate(forward,(chi-GetChi(index)));
}

Real PProteinResidue::GetChi(int index){
  vector<string> rotList = PResources::GetChiIndex(this->getName(), index);
  return PMath::AngleBetweenPlanes(this->getAtomPosition(rotList[0]), this->getAtomPosition(rotList[1]), this->getAtomPosition(rotList[2]), this->getAtomPosition(rotList[3]));
}


