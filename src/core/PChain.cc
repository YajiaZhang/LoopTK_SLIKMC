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

#include <queue>
#include "PBasic.h"
#include "PResources.h"
#include "PConstants.h"
/*I changed here!*/
#include <memory>
using namespace std;
/*I changed here!*/

PChain::PChain() {
  InitChain();
}

void PChain::Obliterate() {
  delete getTopLevelChain();
}

void PChain::InitChain() {
  m_parentChain = NULL;
  m_grid = new PGrid();  // create a new grid for this chain because it has no parent
  m_residues = new vector<PResidue *>;
  m_rotationEvents = new list<PRotateEventHandler *>;
  m_isFinalized = false;
}

PChain::PChain(const string &firstResidueName) {
  InitChain();
  m_residues->push_back(CreateResidue(firstResidueName));
}

PChain::~PChain() {
  list<PChain *> ch_copy = m_children;

  for (list<PChain *>::iterator it = ch_copy.begin(); it!=ch_copy.end();++it) {
    delete *it;
  }
  if (m_parentChain==NULL) {
    for (int i=0;i<m_residues->size();i++) {
      (*m_residues)[i]->DestroyBonds();
    }
    for (int i=0;i<m_residues->size();i++) {
      delete (*m_residues)[i];
    }
    delete m_grid;
    delete m_residues;
    delete m_rotationEvents;
  } else {
    m_parentChain->m_children.remove(this);
    for (int i=0;i<size();i++) {
      getResidue(i)->SetChain(m_parentChain);
    }
  }
}




PChain::PChain(const string &firstResidueName, PResidueSpec &firstResidueSpec) {
  InitChain();
  m_residues->push_back(CreateResidue(firstResidueName,firstResidueSpec));
}


PChain::PChain(PChain *protein,int resStartIndex, int resEndIndex) {
  protein->CheckFinalized();
  m_parentChain = protein;
  m_grid = protein->m_grid;
  m_residues = m_parentChain->m_residues;

  if (m_parentChain->m_parentChain==NULL) {
    m_startIndex = resStartIndex;
    m_endIndex = resEndIndex;
  }
  else {
    m_startIndex = m_parentChain->m_startIndex + resStartIndex;
    m_endIndex = m_parentChain->m_startIndex + resEndIndex;
  }
  m_rotationEvents = protein->m_rotationEvents;
  m_isFinalized = true;
  protein->m_children.push_back(this);
  for (int i=0;i<size();i++) {
    if (getResidue(i)->getChain()!=protein) 
    {
    	//PUtilities::AbortProgram("Invalid subchain created - a subchain cannot be made on top of another subchain.");
    	cout << "Warning: overlapping subchain!" << endl;
    }
    getResidue(i)->SetChain(this);
  }
  getResidue(0)->CacheDOF(m_dofs);
  getResidue(0)->CacheAtoms(m_atomCache);
  return;
}

PChain *PChain::getParent() const {
  return m_parentChain;
}

PChain *PChain::getTopLevelChain() {
  PChain *currChain = this;
  while (true) {
    if (currChain->getParent()==NULL) {
      return currChain;
    }
    currChain = currChain->getParent();
  }
  return NULL;
}



int PChain::NumAtoms(const string &blockType) const {
  CheckFinalized();
  
  //this prevents memory leaks from querying for non-existent blockTypes
  Atom_Cache::const_iterator iter = m_atomCache.find(blockType);
  if (iter != m_atomCache.end()) {
    return iter->second.size();
  } else {
    return 0;
  }
}

int PChain::NumDOF(const string &blockType) const {
  CheckFinalized();
  //this prevents memory leaks from querying for non-existent blockTypes
  DOF_Cache::const_iterator iter = m_dofs.find(blockType);

  if (iter != m_dofs.end()) {
    return iter->second.size();  
  } else {
    return 0;
  }
}

PResidue *PChain::getResidue(int localIndex) {
  if (localIndex<0||localIndex>=m_residues->size())
    PUtilities::AbortProgram(PUtilities::toStr(localIndex) + ": Invalid index into protein chain.");
  if (m_parentChain==NULL) return (*m_residues)[localIndex];
  else {
    int index = m_startIndex+localIndex;
    if (index<=m_endIndex) return (*m_residues)[index];
    else PUtilities::AbortProgram(PUtilities::toStr(localIndex) + ": Invalid index into local protein chain.");
    }
  return NULL; //to quiet compiler
}

int PChain::getResidueIndex(const PAtom* atom){ 
 for (int i=0; i<size(); i++){
    vector<PAtom *> atoms = *(getResidue(i)->getAtoms());
    for (int j=0; j<atoms.size(); j++){
        PAtom* curAtom = atoms[j];
        if (atom == curAtom) return i;
    }
  }
  return -1;
}

PBond *PChain::getBond(int res_index, const string &id1, const string &id2) {
	return getBond(res_index,id1,res_index,id2);
}

PBond *PChain::getBond(int res_index1, const string &id1, int res_index2, const string &id2) {
	PResidue *r1 = getResidue(res_index1);
	PResidue *r2 = getResidue(res_index2);
	PAtom *a = r1->getAtom(id1); 
	if(a==NULL) return NULL;
	return a->getBond(r2,id2);
}

int PChain::size() const { 
  if (m_parentChain==NULL) return m_residues->size(); 
  else return m_endIndex-m_startIndex+1; 
}

PResidue *PChain::AddResidue(const string &resName) {
  PResidue *ret;
  if (m_parentChain==NULL) {
    if (m_residues->size()==0) {
      ret = CreateResidue(resName);
    } else {
      ret = CreateResidue(resName,m_residues->back());
    }
    m_residues->push_back(ret);
    UpdateIndexRangeOnAdd(1);
  } else
    PUtilities::AbortProgram("Cannot add residues to protein subset.");
  return ret;
}

PResidue *PChain::AddResidue(const string &resName, PResidueSpec &resSpec) {
  PResidue *ret;
  if (m_parentChain==NULL) {
    if (m_residues->size()==0) {
      ret = CreateResidue(resName,resSpec);
    }
    else {
      ret = CreateResidue(resName,resSpec,m_residues->back());
    }
    m_residues->push_back(ret);
    UpdateIndexRangeOnAdd(1);
  }
  else
    PUtilities::AbortProgram("Cannot add residues to protein subset.");
  return ret;
}


//this should be removed - does nothing
void PChain::UpdateIndexRangeOnAdd(int amtAdded) {
  for (PChain *curr = this;curr->m_parentChain!=NULL;curr = curr->m_parentChain) {
    curr->m_endIndex+=amtAdded;
  }
}


void PChain::finalize() {
  assert(m_parentChain==NULL);
  /*NOTE: I changed starts: */
//  this->m_startIndex = 0;
//  this->m_endIndex = this->size() - 1;
  /*NOTE: I changed ends! */
  if (m_residues->size()<1) {
    PUtilities::AbortProgram("Cannot finalize an empty chain");
  }

  (*m_residues)[0]->CacheDOF(m_dofs);
  (*m_residues)[0]->CacheAtoms(m_atomCache);

  m_isFinalized = true;
}

void PChain::CheckFinalized() const {
  if (!m_isFinalized) {
    PUtilities::AbortProgram("Invalid operation on an un-finalized chain");
  }
}


bool PChain::IsSubChainOf(const PChain *other) const {
  for (const PChain *c = this; c != NULL; c = c->m_parentChain) {
    if (c == other) return true;
  }
  return false;
}

void PChain::traverseFromStart(AtomFunctor *atomFn, BondFunctor *bondFn)
{
	CheckFinalized();
	if (size() == 0)
	{
		PUtilities::AbortProgram("Cannot traverse -- chain has no residues.");
	}
	else
	{
		PResidue *firstRes = getResidue(0);
		if (firstRes == NULL) {
			PUtilities::AbortProgram("Cannot traverse -- first residue is NULL.");
		}
		else {
			PAtom *a = firstRes->m_atomMap[firstRes->m_shell->getStartAtomID()];
			if (a == NULL) {
				PUtilities::AbortProgram("Cannot traverse -- first residue's start atom is NULL.");
			}
			else {
				a->traverseChain(atomFn, bondFn);
			}
		}
	}
}

void PChain::MultiRotate(vector<ChainMove> &moves) {
	for (int i = 0; i < moves.size(); i++) {
		ChainMove move = moves[i];
		RotateChain(move);
	}
}

void PChain::MultiRotate_noGridUpdate(vector<ChainMove> &moves)
{
	  //TODO: Add code to handle collision grid update
		CheckFinalized();
		assert( moves.size() > 0);
		string blockType = moves[0].blockType;
		BondDirection dir = moves[0].dir;
		for( int i = 0; i < moves.size(); i++)
		{
			assert( moves[i].blockType == blockType);
			assert( moves[i].dir == dir);
			assert( !isNAN( moves[i].degrees));
		}
		vector<PBond*> &dofs = GetDOFs( blockType);

		vector<AtomSet*> atom_sets;
		for( int i = 0; i < moves.size(); i++)
		{
			int DOF_index = moves[i].DOF_index;
			if (DOF_index < 0 || DOF_index >= dofs.size())
			{
				PUtilities::AbortProgram(PUtilities::toStr(DOF_index) + ": Invalid index into dofs for block type: " + blockType + ".");
			}
			PBond* bond = dofs[ DOF_index];
			pair< PAtom*, PAtom*> atom_pair = bond->getAtomPair( dir);
			atom_sets.push_back( atom_pair.first->internalTraverse( atom_pair, this));
		}

		for( int i = 0; i < atom_sets.size(); i++)
		{
			AtomSet* atomset = atom_sets[i];
			int DOF_index = moves[i].DOF_index;
			PBond* bond = dofs[ DOF_index];
			Rotater* rotater = bond->getBondRotater( dir, moves[i].degrees);
			for( AtomSet::iterator it = atomset->begin(); it != atomset->end(); it++)
			{
				rotater->rotateAtom_nonGridUpdate( *it);
			}
			delete rotater;
		}

		for( int i = 0; i < atom_sets.size(); i++)
		{
			delete atom_sets[i];
		}
		return;
}


void PChain::updateAtomsGrid()
{
	int residue_size = this->size();
	for( int i = 0; i < residue_size; i++)
	{
		PResidue* residue = this->getResidue(i);
		vector<PAtom*>* atoms = residue->getAtoms();
		int atom_size = (*atoms).size();
		for( int j = 0; j < atom_size; j++)
		{
			PAtom* atom = (*atoms)[j];
			atom->updateGrid();
		}
	}
}

void PChain::AntiMultiRotate(vector<ChainMove> &moves) {
  for (int i=moves.size()-1;i>=0;i--) {
    ChainMove move = moves[i];
    move.degrees = -move.degrees;
    RotateChain(move);
  }
}

void PChain::RotateChain(const ChainMove &move) {
	CheckFinalized();
	string blockType = move.blockType;
	int DOFindex = move.DOF_index;
	BondDirection dir = move.dir;
	float degrees = move.degrees;
	if (isNAN(degrees)) {
		PUtilities::AbortProgram("Cannot rotate chain with nan degrees.");
	}
	vector<PBond *> &dofs = GetDOFs(blockType);
	if (DOFindex < 0 || DOFindex >= dofs.size()) {
		PUtilities::AbortProgram(PUtilities::toStr(DOFindex)+ ": Invalid index into dofs for block type: "+ blockType + ".");
	}
	PBond *dof = dofs[DOFindex];
	dof->Rotate(dir, degrees, this);
	for (list<PRotateEventHandler *>::iterator it = m_rotationEvents->begin(); it != m_rotationEvents->end(); it++) {
		(*it)->HandleRotation(this, move);
	}
	return;
}

void PChain::RotateChain_noGridUpdate(const ChainMove &move)
{
	CheckFinalized();
	string blockType = move.blockType;
	int DOFindex = move.DOF_index;
	BondDirection dir = move.dir;
	float degrees = move.degrees;
	if (isNAN(degrees)) {
	    PUtilities::AbortProgram("Cannot rotate chain with nan degrees.");
	}
	vector<PBond *> &dofs = GetDOFs(blockType);
	if (DOFindex<0||DOFindex>=dofs.size()) {
	    PUtilities::AbortProgram(PUtilities::toStr(DOFindex) + ": Invalid index into dofs for block type: " + blockType + ".");
	}
	PBond *dof = dofs[DOFindex];
	dof->Rotate_noGridUpdate(dir,degrees,this);
	for (list<PRotateEventHandler *>::iterator it = m_rotationEvents->begin();it!=m_rotationEvents->end();it++) {
		(*it)->HandleRotation(this,move);
	}
	return;
}

void PChain::RotateChain(const string &blockType, int DOFindex, BondDirection dir, float degrees) {
  ChainMove move;
  move.blockType = blockType;
  move.DOF_index = DOFindex;
  move.dir = dir;
  move.degrees = degrees;
  RotateChain(move);
}

void PChain::RotateChain_noGridUpdate(const string &blockType, int DOFindex, BondDirection dir, float degrees) {
  ChainMove move;
  move.blockType = blockType;
  move.DOF_index = DOFindex;
  move.dir = dir;
  move.degrees = degrees;
  RotateChain_noGridUpdate(move);
}

void PChain::AddRotateEventHandler(PRotateEventHandler *handler) {
  CheckFinalized();
  m_rotationEvents->push_back(handler);
}

//aborts if it's not in the list, removes at most one copy
void PChain::RemoveRotateEventHandler(PRotateEventHandler *handler) {
  bool foundIt = false;
  for (list<PRotateEventHandler *>::iterator it = m_rotationEvents->begin();it!=m_rotationEvents->end();it++) {
    if ((*it)==handler) {
      m_rotationEvents->erase(it);
      foundIt = true;
      break;
    }
  }
  if (!foundIt) {
    PUtilities::AbortProgram("Couldn't find rotate event handler to remove.");
  }
}

vector<PBond *>& PChain::GetDOFs(const string &blockType) {
  CheckFinalized();
  return m_dofs[blockType];
}

vector<string> PChain::GetBlockTypes() const {
  CheckFinalized();
  vector<string> ret;
  for (DOF_Cache::const_iterator it = m_dofs.begin(); it != m_dofs.end(); ++it) {
    ret.push_back(it->first);
  }
  return ret;
}

void PChain::RandomizeDOFs(BondDirection dir)
{
  for (DOF_Cache::const_iterator it = m_dofs.begin(); it != m_dofs.end(); ++it) {
    vector<PBond *>& dofs = GetDOFs(it->first);

    for (unsigned i = 0; i < dofs.size(); i++) {
      RotateChain(it->first, i, dir, rand() % 360);
    }
  }
}


void PChain::DetachBlocks(const string &blockType, const string &blockTypeToDetachFrom) {
  DetachBlocks(blockType,blockTypeToDetachFrom,0,size()-1);
}


void PChain::DetachBlocks(const string &blockType, const string &blockTypeToDetachFrom, int resStartIndex, int resEndIndex) {
  CheckFinalized();
  for (int i=resStartIndex;i<=resEndIndex;i++) {
    PResidue *res = getResidue(i);
    //cerr<<i<<"/"<<resEndIndex<<endl;
    res->DetachBlocks(blockType,blockTypeToDetachFrom);
  }
}

void PChain::ReattachBlocks(const string &blockType) {
  ReattachBlocks(blockType,0,size()-1);
}

void PChain::ReattachBlocks(const string &blockType, int resStartIndex, int resEndIndex) {
  CheckFinalized();  
  for (int i=resStartIndex;i<=resEndIndex;i++) {
    PResidue *res = getResidue(i);
    res->ReattachBlocks(blockType);
  }
}

void PChain::ReattachAllBlocks() {
  ReattachAllBlocks(0,size()-1);
}  


void PChain::ReattachAllBlocks(int resStartIndex, int resEndIndex) {
  CheckFinalized();
  for (int i=resStartIndex;i<=resEndIndex;i++) {
    PResidue *res = getResidue(i);
    res->ReattachAllBlocks();
  }
}

bool PChain::InStaticCollision()
{
  if (m_parentChain==NULL)
	  return InStaticCollision(0, m_residues->size()-1);
  else
	  return InStaticCollision(0, m_endIndex - m_startIndex);
}

bool PChain::InStaticCollision(int resIndex1, int resIndex2)
{
  return InCollision(resIndex1, resIndex2, STATIC);
}

bool PChain::InSelfCollision()
{
  if (m_parentChain==NULL)
	  return InSelfCollision(0, m_residues->size()-1);
  else
	  return InSelfCollision(0, m_endIndex - m_startIndex);
}

bool PChain::InSelfCollision(int resIndex1, int resIndex2)
{
  return InCollision(resIndex1, resIndex2, SELF);
}

bool PChain::InAnyCollision()
{
  if (m_parentChain==NULL) return InAnyCollision(0, m_residues->size()-1);
  else return InAnyCollision(0, m_endIndex - m_startIndex);
}

bool PChain::InAnyCollision(int resIndex1, int resIndex2)
{
  return InCollision(resIndex1, resIndex2, EITHER);
}

pair<PAtom *, PAtom *> PChain::FindStaticCollision()
{
  if (m_parentChain==NULL) return FindStaticCollision(0, m_residues->size()-1);
  else return FindStaticCollision(0, m_endIndex - m_startIndex);
}

pair<PAtom *, PAtom *> PChain::FindStaticCollision(int resIndex1, int resIndex2)
{
  return FindCollision(resIndex1, resIndex2, STATIC);
}

pair<PAtom *, PAtom *> PChain::FindSelfCollision()
{
  if (m_parentChain==NULL) return FindSelfCollision(0, m_residues->size()-1);
  else return FindSelfCollision(0, m_endIndex - m_startIndex);
}

pair<PAtom *, PAtom *> PChain::FindSelfCollision(int resIndex1, int resIndex2)
{
  return FindCollision(resIndex1, resIndex2, SELF);
}

pair<PAtom *, PAtom *> PChain::FindAnyCollision()
{
  if (m_parentChain==NULL) return FindAnyCollision(0, m_residues->size()-1);
  else return FindAnyCollision(0, m_endIndex - m_startIndex);
}

pair<PAtom *, PAtom *> PChain::FindAnyCollision(int resIndex1, int resIndex2)
{
  return FindCollision(resIndex1, resIndex2, EITHER);
}

AtomCollisions* PChain::getAllCollidingStatic()
{
  if (m_parentChain==NULL) return getAllCollidingStatic(0, m_residues->size()-1);
  else return getAllCollidingStatic(0, m_endIndex - m_startIndex);
}

AtomCollisions* PChain::getAllCollidingStatic(int resIndex1, int resIndex2)
{
  return getAllColliding(resIndex1, resIndex2, STATIC);
}

AtomCollisions* PChain::getAllCollidingSelf()
{
  if (m_parentChain==NULL) return getAllCollidingSelf(0, m_residues->size()-1);
  else return getAllCollidingSelf(0, m_endIndex - m_startIndex);
}

AtomCollisions* PChain::getAllCollidingSelf(int resIndex1, int resIndex2)
{
  return getAllColliding(resIndex1, resIndex2, SELF);
}

AtomCollisions* PChain::getAllCollidingEither()
{
  if (m_parentChain==NULL) return getAllCollidingEither(0, m_residues->size()-1);
  else return getAllCollidingEither(0, m_endIndex - m_startIndex);
}

AtomCollisions* PChain::getAllCollidingEither(int resIndex1, int resIndex2)
{
  return getAllColliding(resIndex1, resIndex2, EITHER);
}

bool PChain::InCollision(int resIndex1, int resIndex2, CollisionType type)
{
  int resMin = min(resIndex1,resIndex2), resMax = max(resIndex1, resIndex2);

  for (int i = resMin; i <= resMax; i++) {
    if (getResidue(i)->InCollision(type))
    	return true;
  }
  return false;
}

pair<PAtom *, PAtom *> PChain::FindCollision(int resIndex1, int resIndex2, CollisionType type)
{
  int resMin = min(resIndex1, resIndex2), resMax = max(resIndex1, resIndex2);

  for (int i = resMin; i <= resMax; i++) {
    PResidue *curResidue = getResidue(i);
    pair<PAtom *, PAtom *> toTry = curResidue->FindCollision(type);

    if (toTry.first != NULL) {
      assert(toTry.second != NULL);
      return toTry;
    }
  }

  return make_pair<PAtom *, PAtom *>(NULL, NULL);
}

PChainState* PChain::saveChainState()
{
	return this->savePartialChainState( 0, this->size() - 1);
}

PChainState* PChain::savePartialChainState(int index_start, int index_end) {
	assert( index_start >= 0 && index_end < this->size());
	PChainState* state = new PChainState();
	state->index_start = index_start;
	state->index_end = index_end;
	for( int i = index_start; i <= index_end; i++)
	{
		PResidue* residue = this->getResidue(i);
		vector<PAtom*>* atoms = residue->getAtoms();
		vector<Vector3> atom_pos;
		for( int j = 0; j < atoms->size(); j++)
		{
			PAtom* atom = (*atoms)[j];
			atom_pos.push_back( atom->getPos());
		}
		state->atoms_pos.push_back( atom_pos);
	}
	return state;
}

void PChain::restoreChainState( PChainState* state)
{
	assert( state != NULL);
	for( int i = state->index_start; i <= state->index_end; i++)
	{
		PResidue* residue = this->getResidue(i);
		vector<PAtom*>* atoms = residue->getAtoms();
		for( int j = 0; j < atoms->size(); j++)
		{
			PAtom* atom = (*atoms)[j];
			atom->changePosition( state->atoms_pos[ i - state->index_start][j]);
		}
	}
	return;
}

void PChain::restoreChainState_noGridUpdate( PChainState* state)
{
	assert( state != NULL);
	for( int i = state->index_start; i <= state->index_end; i++)
	{
		PResidue* residue = this->getResidue(i);
		vector<PAtom*>* atoms = residue->getAtoms();
		for( int j = 0; j < atoms->size(); j++)
		{
			PAtom* atom = (*atoms)[j];
			atom->changePosition_nonGridUpdate( state->atoms_pos[ i - state->index_start][j]);
		}
	}
	return;
}

void PChain::attachResidues() {
	for( int i = 0; i < this->size(); i++) {
		this->getResidue(i)->SetChain(this);
	}
}

void PChain::attachResidues(int start, int end) {
	assert( start >= 0 && end < this->size() && start <= end);
	for( int i = start; i <= end; i++)
	{
		this->getResidue(i)->SetChain(this);
	}
}

void PChain::detachResidues() {
	for( int i = 0; i < this->size(); i++)
	{
		this->getResidue(i)->SetChain( this->m_parentChain);
	}
}


void PChain::detachResidues(int start, int end) {
	assert( start >= 0 && end < this->size() && start <= end);
	for( int i = start; i <= end; i++)
	{
		this->getResidue(i)->SetChain( this->m_parentChain);
	}
}


bool PChain::areResiduesFullyControl() {
	for( int i = 0; i < this->size(); i++)
	{
		if( this->getResidue(i)->getChain() != this)
			return false;
	}
	return true;
}

void PChain::attachResidue(int index) {
	assert( index >= 0 && index < this->size());
	this->getResidue(index)->SetChain(this);
}

void PChain::detachResidue(int index) {
	assert( index >= 0 && index < this->size());
	this->getResidue(index)->SetChain( this->m_parentChain);
}

//void PChain::rotateBackboneTo(int index_DOF, BondDirection dir, double degree) {
//	int index_residue = index_DOF / 2;
//	double degree_current = this->getDihedralAngle( index_DOF % 2, index_residue);
//	//NOTE:calcualte how many degree to go!
//	double d = degree_current - degree;
//	this->RotateChain("backbone", index_DOF, dir, d);
//	return;
//}
//
//void PChain::rotateBackboneTo_noGridUpdate(int index_DOF, BondDirection dir, double degree)
//{
//	int index_residue = index_DOF / 2;
//	double degree_current = this->getDihedralAngle( index_DOF % 2, index_residue);
//	//NOTE:calcualte how many degree to go!
//	double d = degree_current - degree;
//	this->RotateChain_noGridUpdate("backbone", index_DOF, dir, d);
//	return;
//}

//double PChain::getDihedralAngle(int type, int index_res) {
//	Vector3 b1;
//	Vector3 b2;
//	Vector3 b3;
//	if( PHI == type)
//	{
//		//Rotate Phi, cannot be the Phi in first residue since it just simply rotate the whole chain
//		//residue_prev provides C; residue provides N, Ca, C
//		PResidue* res_prev = NULL;
//		if( index_res == 0)
//		{
//			assert( this->m_parentChain != NULL && this->m_startIndex != 0);
//			res_prev = this->m_parentChain->getResidue( this->m_startIndex - 1);
//		}
//		else
//		{
//			res_prev = this->getResidue(index_res - 1);
//		}
//
//		PResidue* res = this->getResidue(index_res);
//		Vector3 C_prev = res_prev->getAtomPosition( PID::C);
//		Vector3 N = res->getAtomPosition( PID::N);
//		Vector3 Ca = res->getAtomPosition( PID::C_ALPHA);
//		Vector3 C = res->getAtomPosition( PID::C);
//		b1.set( N - C_prev);
//		b2.set( Ca - N);
//		b3.set( C - Ca);
//	}
//	else if( PSI == type)
//	{
//		//Rotate Psi, cannot be the Psi in the last residue
//		//residue provides N, Ca, C, residue_next provides N
////		assert( index_res != this->size() - 1);
//		PResidue* res = this->getResidue(index_res);
//
//		PResidue* res_next = NULL;
//		if( index_res == this->size() - 1)
//		{
//			assert( this->m_parentChain != NULL && this->m_endIndex != (this->m_parentChain->size() - 1));
//			res_next = this->m_parentChain->getResidue( this->m_endIndex + 1);
//		}
//		else
//		{
//			res_next = this->getResidue(index_res + 1);
//		}
//
//		Vector3 N = res->getAtomPosition( PID::N);
//		Vector3 Ca = res->getAtomPosition( PID::C_ALPHA);
//		Vector3 C = res->getAtomPosition( PID::C);
//		Vector3 N_next = res_next->getAtomPosition( PID::N);
//		b1.set( Ca - N);
//		b2.set( C - Ca);
//		b3.set( N_next - C);
//	}
//	else
//	{
//		PUtilities::AbortProgram("Error in selecting type of dihedral angles");
//	}
//	double radian = atan2( b2.norm() * b1.dot(cross( b2, b3)), cross( b1, b2).dot(cross( b2, b3)));
//	double degree = radian / PI * 180.0;
//	return degree;
//}

DihedralAngle* PChain::getDihedralAngleAtResidue(int index_local) {
	//TODO: Should add some error checking mechanism.
	assert( this == this->getTopLevelChain() || this->m_parentChain == this->getTopLevelChain());
	DihedralAngle* pair = new DihedralAngle();

	PChain* chain_toplevel = this;
	int index = index_local;

	if( chain_toplevel->m_parentChain != NULL)
	{
		chain_toplevel = m_parentChain;
		index = index + this->m_startIndex;
	}

    int residue_size = chain_toplevel->size();
	Vector3 b1;
	Vector3 b2;
	Vector3 b3;
	Vector3 b4;
	if (index > 0 && index < residue_size - 1) {
		PResidue* residue_prev = chain_toplevel->getResidue(index - 1);
		PResidue* residue = chain_toplevel->getResidue(index);
		PResidue* residue_next = chain_toplevel->getResidue(index + 1);
		Vector3 C_prev = residue_prev->getAtomPosition(PID::C);
		Vector3 N = residue->getAtomPosition(PID::N);
		Vector3 Ca = residue->getAtomPosition(PID::C_ALPHA);
		Vector3 C = residue->getAtomPosition(PID::C);
		Vector3 N_next = residue_next->getAtomPosition(PID::N);

		//NOTE: pay attention to the direction.
		b1.set(N - C_prev);
		b2.set(Ca - N);
		b3.set(C - Ca);
		b4.set(N_next - C);

		double phi_radian = atan2(b2.norm() * b1.dot(cross(b2, b3)), cross(b1, b2).dot(cross(b2, b3)));
		double psi_radian = atan2(b3.norm() * b2.dot(cross(b3, b4)), cross(b2, b3).dot(cross(b3, b4)));

		pair->phi = phi_radian / PI * 180.0;
		pair->psi = psi_radian / PI * 180.0;
//		pair->residue_name = residue->getName();
	}
	else if (index == 0) {
		//No Phi angle, set to be 360
		PResidue* residue = chain_toplevel->getResidue(index);
		PResidue* residue_next = chain_toplevel->getResidue(index + 1);
		Vector3 N = residue->getAtomPosition(PID::N);
		Vector3 Ca = residue->getAtomPosition(PID::C_ALPHA);
		Vector3 C = residue->getAtomPosition(PID::C);
		Vector3 N_next = residue_next->getAtomPosition(PID::N);
		b2.set(Ca - N);
		b3.set(C - Ca);
		b4.set(N_next - C);
		double psi_radian = atan2(b3.norm() * b2.dot(cross(b3, b4)), cross(b2, b3).dot(cross(b3, b4)));
		//NOTE: change for sidechains
		pair->phi = 360.0;
//		pair->phi = 180.0;
		pair->psi = psi_radian / PI * 180.0;
//		pair->residue_name = residue->getName();
	}
	else if (index == (residue_size - 1)) {
		//No Psi angle, set to be 360
		PResidue* residue_prev = chain_toplevel->getResidue(index - 1);
		PResidue* residue = chain_toplevel->getResidue(index);
		Vector3 C_prev = residue_prev->getAtomPosition(PID::C);
		Vector3 N = residue->getAtomPosition(PID::N);
		Vector3 Ca = residue->getAtomPosition(PID::C_ALPHA);
		Vector3 C = residue->getAtomPosition(PID::C);
		b1.set(N - C_prev);
		b2.set(Ca - N);
		b3.set(C - Ca);
		double phi_radian = atan2(b2.norm() * b1.dot(cross(b2, b3)), cross(b1, b2).dot(cross(b2, b3)));
		pair->phi = phi_radian / PI * 180.0;
		pair->psi = 360;
//		pair->psi = 180.0;
//		pair->residue_name = residue->getName();
	}
	else {
		PUtilities::AbortProgram("Something Wrong Here!");
	}
	return pair;
}

void PChain::getDihedralAngles(vector<DihedralAngle>& angles)
{
	for( int i = 0; i < this->size(); i++)
	{
		angles.push_back( *this->getDihedralAngleAtResidue(i));
	}
	return;
}



void PChain::printDihedralAngles() {
	for( int i = 0; i < this->size(); i++)
	{
		cout << this->getDihedralAngleAtResidue(i)->phi << "\t" << this->getDihedralAngleAtResidue(i)->psi << endl;
	}
	return;
}

AtomCollisions* PChain::getAllColliding(int resIndex1, int resIndex2, CollisionType type)
{
  int resMin = min(resIndex1, resIndex2), resMax = max(resIndex1, resIndex2);
  AtomCollisions *ret = new AtomCollisions();
  
  for (int i = resMin; i <= resMax; i++) {
    PResidue *curResidue = getResidue(i);
    auto_ptr<AtomCollisions> curColl(curResidue->getAllColliding(type));

    ret->insert(curColl->begin(), curColl->end());
  }

  return ret;
}

PAtom* PChain::getAtomAtRes(const string &atomID, int resNum)
{
  PResidue *res = getResidue(resNum);
  if (res == NULL) {
    return NULL;
  }
 
  return res->getAtom(atomID);
}

PAtom *PChain::getAtom(const string &blockType,int index) {
  CheckFinalized();
  vector<PAtom *> &atomVec  = m_atomCache[blockType];
  if (atomVec.size()==0) {
    PUtilities::AbortProgram("No atoms defined for given block type within this chain.");
  }
  return atomVec[index];
}


Vector3 PChain::getAtomPos(const string &blockType,int index) {
  return getAtom(blockType,index)->getPos();
}


pair<int, int> PChain::getTopLevelIndices() const {
  if (m_parentChain==NULL) {
    PUtilities::AbortProgram("No indices into top level chain for master protein.");

  }
  return make_pair(m_startIndex, m_endIndex);
}


PChain *PChain::Clone() {
  PChain *ret = new PChain();
  CloneResiduesIntoChain(ret);
  if (m_parentChain!=NULL) {
    pair<int,int> indices = getTopLevelIndices();
    ret = new PChain(ret,indices.first,indices.second);
  }
  return ret;
}

void PChain::CloneResiduesIntoChain(PChain *other) {
  CheckFinalized();
  for (int i=0;i<m_residues->size();i++) {
    PResidueSpec spec = (*m_residues)[i]->getSpec();
    PResidue *newRes = other->AddResidue((*m_residues)[i]->getResourceName(),spec);
    newRes->setName((*m_residues)[i]->getName());
  }
  other->finalize();
}

int PChain::getResidueLocalIndex(PResidue *res) {
  for (int i=m_startIndex; i<=m_endIndex; i++) {
    if (m_residues->at(i) == res) {
      return i-m_startIndex;
    }
  }
  return -1;
}

int PChain::getResidueGlobalIndex(PResidue *res) {
  for (int i=0; i<m_residues->size(); i++) {
    if (m_residues->at(i)==res) {
      return i;
    }
  }
  return -1;
}

vector<const PAtom*> PChain::extractPath(AtomNode* leaveNode){
    vector<const PAtom*> vec;
    if (leaveNode == NULL) return vec;
    else{
        AtomNode *curNode = leaveNode;
        while(curNode!=NULL && curNode->atom!=NULL){
            vec.push_back(curNode->atom);
            curNode = curNode->parent;
        }
        return vec;
    }
}

vector<const PAtom*> PChain::getShortestPath(const PAtom *a1, ConstAtomSet a2Set, int maxLength, int extendNum){
    AtomNode *rootNode = new AtomNode();
    rootNode->atom = a1;
    rootNode->parent = NULL;
    rootNode->level = 1;

    queue<AtomNode* > bfsQueue;
    queue<AtomNode* > invalidQueue;

    ConstAtomSet atomsSearched;

    AtomNode* curNode;
    const vector<PBond *> *curBonds;
    int curLength = 0;

    bfsQueue.push(rootNode);
    while(!bfsQueue.empty()) {
        if (curLength >= maxLength) return extractPath(NULL);

        curNode = bfsQueue.front();
        bfsQueue.pop();

        /* Did we find the target atom? */
        if (a2Set.find(curNode->atom) != a2Set.end()){
            break;
        }else{/* Enqueue all atoms bonded to the current atom. */
            curBonds = curNode->atom->getBonds();
            for (int i=0; i<curBonds->size(); i++){
                const PAtom *bondedAtom =
                    (PAtom*)PUtilities::PointerThatIsNot((*curBonds)[i]->getAtom1(),
                    (*curBonds)[i]->getAtom2(), curNode->atom);
                if (atomsSearched.find(bondedAtom) == atomsSearched.end()){
                    atomsSearched.insert(bondedAtom);
                    AtomNode* newNode = new AtomNode();
                    newNode->atom = bondedAtom;
                    newNode->parent = curNode;
                    newNode->level = curNode->level+1;
                    bfsQueue.push(newNode);
                  }
            }
        }
        curLength ++;
    }
    vector<const PAtom*> path = extractPath(curNode);
    ConstAtomSet atomsInPath;
    for (int i=0; i<path.size(); i++)
        atomsInPath.insert(path[i]);
    
    //extend front
    if (path.size()==0) return path;
    const PAtom *refAtom = path[0];
    vector<const PAtom*> frontPath;

    bool added=false;
    for (int i=0; i<extendNum; i++){
        curBonds = refAtom->getBonds();
        for (int j=0; j<curBonds->size(); j++){
            const PAtom *bondedAtom = 
                (PAtom*)PUtilities::PointerThatIsNot((*curBonds)[j]->getAtom1(),
                (*curBonds)[j]->getAtom2(), refAtom);
            if (atomsInPath.find(bondedAtom) == atomsInPath.end()){
              atomsInPath.insert(bondedAtom); 
              frontPath.push_back(bondedAtom);
              refAtom = bondedAtom;
              added = true;
              break;
            }
        }
        if (added == false)
            break;
        added = false;
    }

    //extend end
    refAtom = path[path.size()-1];
    vector<const PAtom*> endPath;

    added=false;
    for (int i=0; i<extendNum; i++){
        curBonds = refAtom->getBonds();
        for (int j=0; j<curBonds->size(); j++){

            const PAtom *bondedAtom = 
                (PAtom*)PUtilities::PointerThatIsNot((*curBonds)[j]->getAtom1(),
                (*curBonds)[j]->getAtom2(), refAtom);
            if (atomsInPath.find(bondedAtom) == atomsInPath.end()){
              atomsInPath.insert(bondedAtom); 
              endPath.push_back(bondedAtom);
              refAtom = bondedAtom;
              added = true;
              break;
            }
        }

        if (added == false)
            break;
        added = false;
    }
   
    //add to the path
    path.insert(path.begin(), frontPath.rbegin(), frontPath.rend());
    path.insert(path.end(), endPath.begin(), endPath.end());

 return path;

}

vector<const PAtom*> PChain::getShortestPath(const PAtom *a1, ConstAtomSet a2Set, int maxLength){
    AtomNode *rootNode = new AtomNode();
    rootNode->atom = a1;
    rootNode->parent = NULL;
    rootNode->level = 1;

    queue<AtomNode* > bfsQueue;
    queue<AtomNode* > invalidQueue;

    ConstAtomSet atomsSearched;

    AtomNode* curNode;
    const vector<PBond *> *curBonds;
    int curLength = 0;

    bfsQueue.push(rootNode);
    while(!bfsQueue.empty()) {
        if (curLength >= maxLength) return extractPath(NULL);

        curNode = bfsQueue.front();
        bfsQueue.pop();

        /* Did we find the target atom? */
        if (a2Set.find(curNode->atom) != a2Set.end()){
            return extractPath(curNode);
        }else{/* Enqueue all atoms bonded to the current atom. */
            curBonds = curNode->atom->getBonds();
            for (int i=0; i<curBonds->size(); i++){
                const PAtom *bondedAtom =
                    (PAtom*)PUtilities::PointerThatIsNot((*curBonds)[i]->getAtom1(),
                    (*curBonds)[i]->getAtom2(), curNode->atom);
                if (atomsSearched.find(bondedAtom) == atomsSearched.end()){
                    atomsSearched.insert(bondedAtom);
                   // AtomNode* newNode = new AtomNode();
                   // newNode->atom = bondedAtom;
                   // newNode->parent = curNode;
                   // newNode->level = curNode->level+1;
                        AtomNode* newNode = new AtomNode();
                        newNode->atom = bondedAtom;
                        newNode->parent = curNode;
                        newNode->level = curNode->level+1;

                bfsQueue.push(newNode);

                    /*
                     * @todo: modified the distane code
                     */    
                    double dx = curNode->atom->getPos().x-bondedAtom->getPos().x;
                    double dy = curNode->atom->getPos().y-bondedAtom->getPos().y;
                    double dz = curNode->atom->getPos().z-bondedAtom->getPos().z;
                    double dis = sqrt(dx*dx+dy*dy+dz*dz );
                   if (dis>2){
                       // bfsQueue.push(newNode);
                        cout << "getID:"<< (bondedAtom->getID())<<endl;
                        cout << "getName:"<<(bondedAtom->getName())<<endl;

                    }
                    

                  }
            }
        }
        curLength ++;
    }
}


vector<const PAtom*> PChain::getShortestPath(const PAtom *a1, const PAtom *a2, int maxLength){
    AtomNode *rootNode = new AtomNode();
    rootNode->atom = a1;
    rootNode->parent = NULL;
    rootNode->level = 1;

    queue<AtomNode* > bfsQueue;
    ConstAtomSet atomsSearched;

    AtomNode* curNode;
    const vector<PBond *> *curBonds;
    int curLength = 0;

    bfsQueue.push(rootNode);
    while(!bfsQueue.empty()) {
        if (curLength >= maxLength) return extractPath(NULL);

        curNode = bfsQueue.front();
        bfsQueue.pop();

        /* Did we find the target atom? */
        if (curNode->atom == a2){
            return extractPath(curNode);
        }else{/* Enqueue all atoms bonded to the current atom. */
            curBonds = curNode->atom->getBonds();
            for (int i=0; i<curBonds->size(); i++){
                const PAtom *bondedAtom =
                    (PAtom*)PUtilities::PointerThatIsNot((*curBonds)[i]->getAtom1(),
                    (*curBonds)[i]->getAtom2(), curNode->atom);
                if (atomsSearched.find(bondedAtom) == atomsSearched.end()){
                    atomsSearched.insert(bondedAtom);
                    AtomNode* newNode = new AtomNode();
                    newNode->atom = bondedAtom;
                    newNode->parent = curNode;
                    newNode->level = curNode->level+1;
                    bfsQueue.push(newNode);
                }
            }
        }
        curLength ++;
    }
}

void PChain::inactivateResidue(int Rid1, int Rid2) {
  PResidue *res;
  for (int rid=Rid1; rid<=Rid2; ++rid) {
    res = getResidue(rid);
    res->inactivate();
  }
}

void PChain::activateResidue(int Rid1, int Rid2) {
  PResidue *res;
  for (int rid=Rid1; rid<=Rid2; ++rid) {
    res = getResidue(rid);
    res->activate();
  }
}

void PChain::getBackbonePositions( vector<Vector3>& positions)
{
	for( int i = 0; i < this->size(); i++)
	{
		Vector3 N = this->getAtomAtRes( PID::N, i)->getPos();
		Vector3 Ca = this->getAtomAtRes( PID::C_ALPHA, i)->getPos();
		Vector3 C = this->getAtomAtRes( PID::C, i)->getPos();
		positions.push_back( N);
		positions.push_back( Ca);
		positions.push_back( C);
	}
	return;
}

void PChain::printBackbonePosition()
{
	for( int i = 0; i < this->size(); i++)
	{
		Vector3 N = this->getAtomAtRes( PID::N, i)->getPos();
		Vector3 Ca = this->getAtomAtRes( PID::C_ALPHA, i)->getPos();
		Vector3 C = this->getAtomAtRes( PID::C, i)->getPos();
		cout << "N:\t" << N << endl;
		cout << "Ca:\t" << Ca << endl;
		cout << "C:\t" << C << endl;
	}
	return;
}

void PChain::printBackboneCollisionGrid()
{
	for( int i = 0; i < this->size(); i++)
	{
		Vector3 N = this->getAtomAtRes( PID::N, i)->getGridPos();
		Vector3 Ca = this->getAtomAtRes( PID::C_ALPHA, i)->getGridPos();
		Vector3 C = this->getAtomAtRes( PID::C, i)->getGridPos();
		cout << "\t" << N << endl;
		cout << "\t" << Ca << endl;
		cout << "\t" << C << endl;
	}
}

vector<Vector3>* PChain::getBackboneCollisionGrid()
{
	vector<Vector3>* grids = new vector<Vector3>();
	for( int i = 0; i < this->size(); i++)
	{
		Vector3 N = this->getAtomAtRes( PID::N, i)->getGridPos();
		Vector3 Ca = this->getAtomAtRes( PID::C_ALPHA, i)->getGridPos();
		Vector3 C = this->getAtomAtRes( PID::C, i)->getGridPos();
//		cout << "\t" << N << endl;
//		cout << "\t" << Ca << endl;
//		cout << "\t" << C << endl;

		grids->push_back( N);
		grids->push_back( Ca);
		grids->push_back( C);
	}
	return grids;
}
