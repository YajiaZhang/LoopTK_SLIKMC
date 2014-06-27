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
/*I changed here!*/
#include <memory>
using namespace std;
/*I changed here!*/

PBlockShell::PBlockShell(const string &name, const string &type) {
  m_name = name;
  m_type = type;
}

void PBlockShell::addAtom(const string &id, const string &name, Vector3 relativePosition) { 
  StringPair toInsert = make_pair(id, name);

  pair<HASH_MAP_STR(string)::iterator,bool> ret = m_atoms.insert(toInsert);
  if (!ret.second)
    PUtilities::AbortProgram("Atom ID clash in block \"" + m_name + "\" of atom ID \"" + id + "\"");
  else
    m_defaultRelPositions.addAtom(id,relativePosition);
}


void PBlockShell::addBond(const string &id1, const string &id2, bool isDOF) {
  addBond(PBondShell(id1,id2,isDOF));
}

void PBlockShell::addBond(const PBondShell &pbs) {
  m_bondSpecs.push_back(pbs);
}


void PAtomPositionsSpec::addAtom(const string &id, const Vector3 &position) {
  pair<string, Vector3> toInsert = make_pair(id, position);
  pair<AtomPositions::iterator,bool> ret = m_atomPositions.insert(toInsert);

  if (!ret.second)
    PUtilities::AbortProgram("Atom ID clash in block spec for atomID \"" + id + "\"");
}

Vector3 PAtomPositionsSpec::getAtomPosition(const string &id) const {
  AtomPositions::const_iterator found = m_atomPositions.find(id);

  if (found == m_atomPositions.end()) {
    PUtilities::AbortProgram(id + " does not exist in block spec!");
  } else {
    return found->second;
  }
}

//use relative coordinates as global coordinates for atom positions
PBlock::PBlock(PResidue *myRes, const string &shellName) {
  m_shell = PResources::GetBlockShell(shellName);
  CreateBlock(myRes,*(m_shell->getDefaultPositions()));
}

PBlock::PBlock(PResidue *myRes, const string &shellName, const PAtomPositionsSpec &spec) {
  m_shell = PResources::GetBlockShell(shellName);
  CreateBlock(myRes,spec);
}


PBlock::PBlock(PResidue *myRes, const string &shellName, const PAtomPositionsSpec &spec, PBlock *toConnect) {
  m_shell = PResources::GetBlockShell(shellName);
  m_residue = myRes;
  CreateBlock(myRes,spec);
  PBlockConnection *connector = PResources::GetBlockConnection(toConnect->getName(),m_shell->getName());
  connector->BondBlocksTogether(toConnect,this);
}


//toConnect's positions are already defined
//this function should look up the necessary block connection and connect itself
PBlock::PBlock(PResidue *myRes, const string &shellName, PBlock *toConnect) {
  m_shell = PResources::GetBlockShell(shellName);
  m_residue = myRes;
  m_isOn = true;
  PBlockConnection *connector = PResources::GetBlockConnection(toConnect->getName(),m_shell->getName());
  connector->GenerateNewBlock(toConnect,this);
  
}


void PBlock::ApplyTransform(Matrix4 transform) {
  for(HASH_MAP_STR(PAtom *)::iterator it = m_atoms.begin();it!=m_atoms.end();++it) {
    it->second->ApplyTransform(transform);
  }
}

void PBlock::CreateBlock(PResidue *myRes, const PAtomPositionsSpec &spec) {
  m_residue = myRes;
  m_isOn = true;
  HASH_MAP_STR(string) *atoms = m_shell->getAtoms();
  for(HASH_MAP_STR(string)::iterator it = atoms->begin();it!=atoms->end();++it) {
    string newAtomShell = it->second;
    Vector3 position = spec.getAtomPosition(it->first);
    PAtom *newAtom = new PAtom(this,newAtomShell,position,it->first);
    m_atoms[it->first] = newAtom;
  }
  vector<PBondShell> *bonds = m_shell->getBonds();
  for(int i=0;i<bonds->size();i++) {
    string atom1 = (*bonds)[i].getAtomId1();
    string atom2 = (*bonds)[i].getAtomId2();
    new PBond(m_atoms[atom1],m_atoms[atom2],(*bonds)[i].isDOF());
  }
}

void PBlock::AddAtomsToMap(HASH_MAP_STR(PAtom *) &atomCache) const {
  for(HASH_MAP_STR(PAtom *)::const_iterator it = m_atoms.begin(); it != m_atoms.end(); ++it) {
    atomCache[it->first] = it->second;
  }
}

PBlock::~PBlock() {
  for(HASH_MAP_STR(PAtom *)::iterator it = m_atoms.begin();it!=m_atoms.end();++it) {
    delete it->second;
  }
}

bool PBlock::InStaticCollision() const {
  return InCollision(STATIC);
}

bool PBlock::InSelfCollision() const {
  return InCollision(SELF);
}

bool PBlock::InAnyCollision() const {
  return InCollision(EITHER);
}


pair<PAtom *, PAtom *> PBlock::FindStaticCollision() const {
  return FindCollision(STATIC);
}

pair<PAtom *, PAtom *> PBlock::FindSelfCollision() const {
  return FindCollision(SELF);
}

pair<PAtom *, PAtom *> PBlock::FindAnyCollision() const {
  return FindCollision(EITHER);
}

AtomCollisions* PBlock::getAllCollidingStatic() const {
  return getAllColliding(STATIC);
}

AtomCollisions* PBlock::getAllCollidingSelf() const {
  return getAllColliding(SELF);
}

AtomCollisions* PBlock::getAllCollidingEither() const {
  return getAllColliding(EITHER);
}
  
PChain* PBlock::getChain() const {
  assert(m_residue!=NULL);
  return m_residue->getChain();
}




void PBlock::DeactivateBlock(PBlock *block, PBlock *from) {
  block->m_isOn = false;
  block->m_blockReconnector = new PBlockReconnector(from,block);
  for(HASH_MAP_STR(PAtom *)::iterator it=block->m_atoms.begin();it!=block->m_atoms.end();++it) {
    it->second->removeMeFromGrid();
  }

}

void PBlock::ActivateBlock(PBlock *block, PBlock *from) {
  block->m_isOn = true;
  for(HASH_MAP_STR(PAtom *)::iterator it=block->m_atoms.begin();it!=block->m_atoms.end();++it) {
    it->second->insertMeIntoGrid();
  }
}

void PBlock::ReattachBlock(PBlock *block, PBlock *from ) {
  block->ActivateAndDestroyReattachment();
}


void PBlock::ActivateAndDestroyReattachment() {
  //assert(!isOn());
  m_blockReconnector->ReconnectBlocks();
  delete m_blockReconnector;
}

bool PBlock::InCollision(CollisionType type) const {
  pair<PAtom *,PAtom *> coll = FindCollision(type);

  if (coll.first==NULL) {
    assert(coll.second==NULL);
    return false;
  } else {
    return true;
  }
}

pair<PAtom *, PAtom *> PBlock::FindCollision(CollisionType type) const {
  if (isOn()) {
    for(HASH_MAP_STR(PAtom *)::const_iterator it = m_atoms.begin(); it != m_atoms.end(); ++it) {
      PAtom *a = it->second, *someAtom;

      if (type == STATIC)
    	  someAtom = a->FindStaticCollision();
      else
    	  if (type == SELF)
    		  someAtom = a->FindSelfCollision();
      else
    	  someAtom = a->FindAnyCollision();

      if (someAtom != NULL && someAtom->getParentBlock()->isOn())
    	  return make_pair(a,someAtom);
    }
  }
  return make_pair<PAtom *, PAtom *>(NULL,NULL);
}

AtomCollisions* PBlock::getAllColliding(CollisionType type) const {
  AtomCollisions *ret = new AtomCollisions();
  AtomSet *cur_collisions;

  if (isOn()) {
    for(HASH_MAP_STR(PAtom *)::const_iterator it=m_atoms.begin();it!=m_atoms.end();it++) {
      PAtom *a = it->second;

      if (type == STATIC) cur_collisions = a->getAllCollidingStatic();
      else if (type == SELF) cur_collisions = a->getAllCollidingSelf();
      else cur_collisions = a->getAllCollidingEither();

      auto_ptr<AtomSet> cc(cur_collisions);

      for(AtomSet::const_iterator it = cc->begin(); it != cc->end(); it++) {
        ret->insert(make_pair(a,*it));
      }
    }
  }

  return ret;
}


PBlock *PBlock::FindFirstBlockOfType(const string &blockType) {
  for(int i=0;i<m_connectedBlocks.size();i++) {
    PBlock *tester = m_connectedBlocks[i];
    if (tester->getType()==blockType&&getParentResidue()==tester->getParentResidue()) return tester;
  }
  return NULL;
}

void PBlock::Detach(const string &blockType) {
  if (!isOn()) return;   // Can only detach once.
  PBlock *toDetachFrom = FindFirstBlockOfType(blockType);
  if (toDetachFrom==NULL) {
    return;
  }
  TraverseBlocks(toDetachFrom,&DeactivateBlock,true);
}


void PBlock::Reattach() {
  if (isOn()) return;
  PBlock *toAttachTo = m_blockReconnector->getAnchor();
  if (!toAttachTo->isOn()) {
    toAttachTo->Reattach();
    return;
  }
  TraverseBlocks(toAttachTo,&ActivateBlock,false);
  TraverseBlocks(toAttachTo,&ReattachBlock,false);
}


void PBlock::TraverseBlocks(PBlock *from, void (*ManipulateBlockFunc)(PBlock *block, PBlock *from), bool stopAtOffBlock) {
  ManipulateBlockFunc(this,from);
  for(int i=0;i<m_connectedBlocks.size();i++) {
    PBlock *toTraverse = m_connectedBlocks[i];
    if (getParentResidue()==toTraverse->getParentResidue()&&toTraverse!=from&&(toTraverse->isOn()||!stopAtOffBlock)) {
      toTraverse->TraverseBlocks(this,ManipulateBlockFunc,stopAtOffBlock);
    }
  }
}

const PSpaceManager* PBlock::getSpaceManager() const {
  return getChain()->getSpaceManager();
}
