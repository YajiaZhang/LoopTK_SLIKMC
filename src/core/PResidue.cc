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

PResidueShell::PResidueShell(const string &name, const string &headBlockId, const string &headBlockName, 
  const string &tailBlockId, const string &tailBlockName, const string &startAtomID) {
  ConstructorHelper(name,headBlockId,headBlockName,tailBlockId,tailBlockName,startAtomID);
} 

PResidueShell::PResidueShell(const string &name,
                             const string &coreBlockId,
                             const string &coreBlockName,
                             const string &startAtomID) {
  ConstructorHelper(name,coreBlockId,coreBlockName,coreBlockId,coreBlockName,startAtomID);
} 


void PResidueShell::ConstructorHelper(const string &name,
                                      const string &headBlockId,
                                      const string &headBlockName, 
                                      const string &tailBlockId,
                                      const string &tailBlockName,
                                      const string &startAtomID) {
  m_name = name;
  m_startAtomID = startAtomID;
  m_headBlock = headBlockId;
  m_tailBlock = tailBlockId;
  m_blocks[headBlockId] = headBlockName;
  if (headBlockId == tailBlockId) { 
    if (headBlockName != tailBlockName) {
      PUtilities::AbortProgram("Head and tail block Id's are same in PResidueShell \"" +
                               name + "\", but their atom types are different");
    } else {
      m_blocks[tailBlockId] = tailBlockName;
    }
  }

}

void PResidueShell::addBlock(const string &id, const string &name) {
  StringPair toInsert = make_pair(id, name);
  pair<HASH_MAP_STR(string)::iterator, bool> ret = m_blocks.insert(toInsert);

  if (!ret.second) {
    PUtilities::AbortProgram("Block ID clash in block \""+m_name+"\" of block ID \"" + id + "\"");
  }
}


//connection idToDefine TO idDefined
void PResidueShell::addBlockConnection(const string &idDefined, const string &idToDefine) {
  m_blockBonds.push_back(StringPair(idDefined,idToDefine));
}

//uses default positions (starts with default for the head block, then connects everything to it)
PResidue::PResidue(PChain *loop, const string &shellName) {

	this->type_sidechain = -1;

  SetChain(loop);  
  m_shell = PResources::GetResidueShell(shellName);
  m_blocks[m_shell->getHeadId()] = new PBlock(this,m_shell->getHeadBlockShell());
  ConstructAroundHead(loop);
  EstablishLink(NULL);
  setName(getResourceName());
  finalize();
}



PResidue::PResidue(PChain *loop, const string &shellName, PResidueSpec &spec) {
  m_shell = PResources::GetResidueShell(shellName);
  BuildWithPositions(loop,spec);
  EstablishLink(NULL);
  setName(getResourceName());
  finalize();
  AddAuxInfo(spec);
}


PResidue::PResidue(PChain *loop, const string &shellName, PResidueSpec &spec, PResidue *toConnect) {
  m_shell = PResources::GetResidueShell(shellName);
  BuildWithPositions(loop,spec);
  EstablishLink(toConnect);
  PBlock *other = toConnect->getTailBlock();
  PBlock *myHead = m_blocks[m_shell->getHeadId()];
  PBlockConnection *connector = PResources::GetBlockConnection(other->getName(),myHead->getName());
  connector->BondBlocksTogether(other,myHead);
  setName(getResourceName());
  finalize();
  AddAuxInfo(spec);
}

void PResidue::DestroyBonds() {
  for(HASH_MAP_STR(PAtom *)::iterator it = m_atomMap.begin();it!=m_atomMap.end();it++) {
    it->second->DestroyBonds();
  }
}

void PResidue::BuildWithPositions(PChain *loop, const PResidueSpec &spec)
{
  /* Set this residue's parent chain. */
  SetChain(loop);
  m_pdb_id = spec.pdb_id;

  if (!IsFullyDefined(m_shell->getHeadBlockShell(),spec)) {
    PUtilities::AbortProgram(m_shell->getHeadBlockShell() + " is not fully defined in a residue, and the chain cannot be constructed");
  } else {
    m_blocks[m_shell->getHeadId()] = new PBlock(this,m_shell->getHeadBlockShell(),spec.Atom_Positions);
  }
  
  HASH_MAP_STR(string) *blockShells = m_shell->getBlocks();
  vector<StringPair> *blockBonds = m_shell->getBlockBonds();
  while(true) {
    bool somethingHappenned = false; //in case not all block connections were properly defined, it won't loop forever (just connects all possible blocks)
    for(int i=0;i<blockBonds->size();i++) {
      StringPair bond = (*blockBonds)[i];
      PBlock *defined = m_blocks[bond.first];
      if (defined!=NULL&&m_blocks[bond.second]==NULL) {
        string shellName = (*blockShells)[bond.second];
        if (!IsFullyDefined(shellName,spec)) {
          LoopTK::DisplayWarning(shellName + " block of residue " + m_shell->getName() +
              " is not fully defined; using canonical representation.");
          m_blocks[bond.second] = new PBlock(this,shellName,defined);
        } else {
          m_blocks[bond.second] = new PBlock(this,shellName,spec.Atom_Positions,defined);
        }
        somethingHappenned = true;
      }
    }
    if (!somethingHappenned) break;
  }

  
  /*
    StringMap *blockShells = m_shell->getBlocks();
    for(StringMap::const_iterator it = blockShells->begin(); it != blockShells->end(); ++it) {
    string newBlockShell = it->second;
    PBlockShell *curShell = PResources::GetBlockShell(it->second);
  
  if (IsFullyDefined(curShell, spec)) {
  m_blocks[it->first] = new PBlock(this, newBlockShell, spec.Atom_Positions);
  } else {
  LoopTK::DisplayWarning(it->second + " block of residue is not fully defined; using default positions.");
  m_blocks[it->first] = new PBlock(this, newBlockShell);
  }
    }
  */
  
  
}

PResidue::PResidue(PChain *loop, const string &shellName, PResidue *toConnect) {
  m_shell = PResources::GetResidueShell(shellName);
  PBlock *tail = toConnect->getTailBlock();
  SetChain(loop);
  EstablishLink(toConnect);
  m_blocks[m_shell->getHeadId()] = new PBlock(this,m_shell->getHeadBlockShell(),tail);
  ConstructAroundHead(loop);
  setName(getResourceName());
  finalize();
}

void PResidue::AddAuxInfo(PResidueSpec &spec) {
	vector<PAtom *> *atoms = getAtoms();
	for(unsigned i=0;i<atoms->size();i++) {
		PAtom *a = (*atoms)[i];
		if(spec.Atom_Occ.find(a->getID())!=spec.Atom_Occ.end()) {
			a->setOccupancy(spec.Atom_Occ[a->getID()]);
		}
		if(spec.Atom_Temp.find(a->getID())!=spec.Atom_Temp.end()) {
			a->setTempFactor(spec.Atom_Temp[a->getID()]);
		}
	}
}

void PResidue::EstablishLink(PResidue *prior) {
  if (prior!=NULL) {
    prior->m_nextRes = this;
    m_prevRes = prior;
  } else {
    m_prevRes = NULL;
  }
}

void PResidue::ConstructAroundHead(PChain *loop) {
  SetChain(loop);
  HASH_MAP_STR(string) *blockShells = m_shell->getBlocks();
  vector<StringPair> *blockBonds = m_shell->getBlockBonds();
  while(true) {
    bool somethingHappenned = false; //in case not all block connections were properly defined, it won't loop forever (just connects all possible blocks)
    for(int i=0;i<blockBonds->size();i++) {
      StringPair bond = (*blockBonds)[i];
      PBlock *defined = m_blocks[bond.first];
      if (defined!=NULL&&m_blocks[bond.second]==NULL) {
  //  cerr<<bond.first<<endl<<bond.second<<endl;
  m_blocks[bond.second] = new PBlock(this,(*blockShells)[bond.second],defined);
  somethingHappenned = true;
      }
    }
    if (!somethingHappenned) break;
  }
}


void PResidue::finalize() {
  for(HASH_MAP_STR(PBlock *)::iterator it = m_blocks.begin();it!=m_blocks.end();++it) {
    it->second->AddAtomsToMap(m_atomMap);
    m_blockTypeCache[it->second->getType()].push_back(it->second);
  }
  for(HASH_MAP_STR(PAtom *)::iterator it = m_atomMap.begin();it!=m_atomMap.end();++it) {
    m_atomCache.push_back(it->second);
  }
}

PResidue::~PResidue() {
  for(HASH_MAP_STR(PBlock *)::iterator it = m_blocks.begin();it!=m_blocks.end();++it) {
    delete it->second;
  }
}

bool PResidue::InStaticCollision() const
{
  return InCollision(STATIC);
}

bool PResidue::InSelfCollision() const
{
  return InCollision(SELF);
}

bool PResidue::InAnyCollision() const
{
  return InCollision(EITHER);
}

pair<PAtom *, PAtom *> PResidue::FindStaticCollision() const
{
  return FindCollision(STATIC);
}

pair<PAtom *, PAtom *> PResidue::FindSelfCollision() const
{
  return FindCollision(SELF);
}

pair<PAtom *, PAtom *> PResidue::FindAnyCollision() const
{
  return FindCollision(EITHER);
}

AtomCollisions* PResidue::getAllCollidingStatic() const
{
  return getAllColliding(STATIC);
}

AtomCollisions* PResidue::getAllCollidingSelf() const
{
  return getAllColliding(SELF);
}

AtomCollisions* PResidue::getAllCollidingEither() const
{
  return getAllColliding(EITHER);
}

PBond *PResidue::getBond(const string &id1, const string &id2) {
	PAtom *a1 = getAtom(id1);
	if(a1==NULL) return NULL;
	return a1->getBond(this,id2);
}

struct DOFCacher:BondFunctor {
public:
  DOFCacher(DOF_Cache *d) { DOFcache = d; }
  void operator()(PBond *bond, PAtom *atomFrom, PAtom *atomTo) {
    bond->setDirection(atomTo);
    if (bond->isDOF()) {
      (*DOFcache)[atomTo->getParentBlock()->getType()].push_back(bond);
      if (atomFrom->getParentResidue()==atomTo->getParentResidue()) {
  atomFrom->getParentResidue()->CacheSingleDOF(atomTo->getParentBlock()->getType(),bond);
      }
    }
  }

  DOF_Cache *DOFcache;
};

bool PResidue::applyRotamer(const int& index, const vector<double>& rotamer) {

//	cout << "calling PProteinResidueinResidue::applyRotamer" << endl;
	this->angles_sidechain.resize(rotamer.size());
	for (int i = 0; i < rotamer.size(); i++) {
		this->angles_sidechain[i] = rotamer[i];
	}
	this->type_sidechain = index;
	//see how many chi angles in this give residue.
	unsigned chiMax = PResources::numChiIndices(this->getName());

//	  assert( chiMax == rotamer.size());
	if (!(chiMax == rotamer.size())) {
		cout << "chiMax:" << chiMax << endl;
		cout << "rotamer.size:" << rotamer.size() << endl;
		abort();
	}

	//store the name of atoms that defines the chi angle.
	vector<string> rotList;

	//apply rotamer angles
	for (unsigned i = 1; i <= chiMax; i++) {
		//get the atoms used to define dihedral angles
		rotList = PResources::GetChiIndex(this->getName(), i);
		//calculate current dihedral angle (chi angle) and calculate the delta needed to achive the goal.
		Real curDihedral = PMath::AngleBetweenPlanes(
				this->getAtomPosition(rotList[0]),
				this->getAtomPosition(rotList[1]),
				this->getAtomPosition(rotList[2]),
				this->getAtomPosition(rotList[3]));

	    this->getDOF("sidechain",rotList[1],rotList[2])->Rotate_noGridUpdate(forward,(curDihedral - rotamer[i-1]), this->getChain());
	}
	return true;
}

void PResidue::getChiAngles( vector<double>& chis) {
	assert( chis.size() == 0);
	unsigned chiMax = PResources::numChiIndices(this->getName());
	vector<string> rotList;
	for (unsigned i = 1; i <= chiMax; i++) {
		//get the atoms used to define dihedral angles
		rotList = PResources::GetChiIndex(this->getName(), i);
		//calculate current dihedral angle (chi angle) and calculate the delta needed to achive the goal.
		Real curDihedral = PMath::AngleBetweenPlanes(
				this->getAtomPosition(rotList[0]),
				this->getAtomPosition(rotList[1]),
				this->getAtomPosition(rotList[2]),
				this->getAtomPosition(rotList[3]));

		chis.push_back( curDihedral);
	}
}


void PResidue::calChi(vector<double>& chis) {
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
	    chis.push_back( curDihedral);
	  }
}


void PResidue::getSideChainAngles( int& index, vector<double>& angles){
	index = this->type_sidechain;
	for( int i = 0; i < this->angles_sidechain.size(); i++) {
		angles.push_back( this->angles_sidechain[i]);
	}
	return;
}


void PResidue::CacheSingleDOF(const string &blockType,PBond *bond) {
  pair<string,string> query;
  query.first = bond->getAtom1()->getID();
  query.second = bond->getAtom2()->getID();
  m_dofs[blockType][query] = bond;
}

PBond *PResidue::getDOF(const string &blockType, const string &atomId1, const string &atomId2) {
  if (m_dofs.find(blockType)==m_dofs.end()) return NULL;
  else
    {
      pair<string,string> query;
      query.first = atomId1;
      query.second = atomId2;
      if (m_dofs[blockType].find(query)==m_dofs[blockType].end()) return NULL;
      else
  return m_dofs[blockType][query];
    }
}

void PResidue::RandomizeDOFs(const string &blockType) {
  RandomizeDOFs(blockType,getChain());
}

void PResidue::RandomizeDOFs(const string &blockType, PChain *rootChain) {
  if (m_dofs.find(blockType)==m_dofs.end()) return;
  for(ID_To_DOF_Map::iterator it = m_dofs[blockType].begin();it!=m_dofs[blockType].end();++it) {
    it->second->Rotate(forward,rand()%360,rootChain);
  }
}

const ID_To_DOF_Map *PResidue::getDOFsForBlock(const string &blockType) {
  if (m_dofs.find(blockType)==m_dofs.end()) {
    return NULL;
  } else return &(m_dofs[blockType]);
}

void PResidue::CacheDOF(DOF_Cache &DOFcache) {
  PAtom *a = m_atomMap[m_shell->getStartAtomID()];
  DOFCacher dcf(&DOFcache);
  a->traverseBonds(&dcf);
}

struct AtomCacher:AtomFunctor {
public:
  AtomCacher(Atom_Cache *a) { AtomCache = a; }

  void operator()(PAtom *atom, PBond *bondFrom) {
    (*AtomCache)[atom->getParentBlock()->getType()].push_back(atom);
  }

  Atom_Cache *AtomCache;

};


void PResidue::CacheAtoms(Atom_Cache &AtomCache) {
  PAtom *a = m_atomMap[m_shell->getStartAtomID()];
  AtomCacher acf(&AtomCache);
  a->traverseAtoms(&acf);
}

PChain* PResidue::getChain() const {
  return m_loop;
}


void PResidue::DetachBlocks(const string &blockType, const string &blockToDetachFrom) {
  if (m_blockTypeCache.find(blockType)==m_blockTypeCache.end()) return;
  vector<PBlock *> &blocks = m_blockTypeCache[blockType];
  for(int i=0;i<blocks.size();i++) {
    //cerr<<"atoms:"<<blocks[i]->size()<<","<<blocks[i]->getName()<<endl;
    blocks[i]->Detach(blockToDetachFrom);
  }
}

void PResidue::ReattachBlocks(const string &blockType) {
  if (m_blockTypeCache.find(blockType)==m_blockTypeCache.end()) return;
  vector<PBlock *> &blocks = m_blockTypeCache[blockType];
  for(int i=0;i<blocks.size();i++) {
    blocks[i]->Reattach();
  }
}

void PResidue::ReattachAllBlocks() {
  for(BlockTypeCache::iterator it = m_blockTypeCache.begin();it!=m_blockTypeCache.end();++it) {
    ReattachBlocks(it->first);
  }
}




/**
 * Collision detection helper
 * methods
 */

bool PResidue::InCollision(CollisionType type) const {
  pair<PAtom *,PAtom *> coll = FindCollision(type);
  if (coll.first==NULL) {
    assert(coll.second==NULL);
    return false;
  }
  else
    return true;
}

pair<PAtom *,PAtom *> PResidue::FindCollision(CollisionType type) const {
  for(HASH_MAP_STR(PBlock *)::const_iterator it = m_blocks.begin();it!=m_blocks.end();++it) {
    PBlock *b = it->second;
    pair<PAtom *, PAtom *> toTry;

    if (type == STATIC)
    	toTry = b->FindStaticCollision();
    else if (type == SELF)
    	toTry = b->FindSelfCollision();
    else
    	toTry = b->FindAnyCollision();

    if (toTry.first != NULL)
    	return toTry;
  }
  return make_pair<PAtom *,PAtom *>(NULL,NULL);
}


AtomCollisions* PResidue::getAllColliding(CollisionType type) const {  
  AtomCollisions *ret = new AtomCollisions(), *cur_collisions;

  for(HASH_MAP_STR(PBlock *)::const_iterator it = m_blocks.begin();it!=m_blocks.end();++it) {
    PBlock *cur_block = it->second;

    if (type == STATIC) {
      cur_collisions = cur_block->getAllCollidingStatic();
    } else if (type == SELF) {
      cur_collisions = cur_block->getAllCollidingSelf();
    } else {
      cur_collisions = cur_block->getAllCollidingEither();
    }

    auto_ptr<AtomCollisions> cc(cur_collisions);
    ret->insert(cc->begin(), cc->end());
  }
  return ret;
}


const PSpaceManager* PResidue::getSpaceManager() const {
  return getChain()->getSpaceManager();
}


PResidueSpec PResidue::getSpec() {
  PResidueSpec ret;
  ret.pdb_id = m_pdb_id;

  for(HASH_MAP_STR(PAtom *)::iterator it = m_atomMap.begin();it!=m_atomMap.end();++it) {
    ret.Atom_Positions.addAtom(it->first,it->second->getPos());
    ret.Atom_Occ[it->first] = it->second->getOccupancy();
    ret.Atom_Temp[it->first] = it->second->getTempFactor();
  }
  return ret;
}

void PResidue::setPositions(const PAtomPositionsSpec &spec) {
  for(HASH_MAP_STR(PAtom *)::iterator it = m_atomMap.begin(); it!=m_atomMap.end(); ++it) {
    PAtom *atom = it->second;
    if (spec.contains(it->first)) {
      atom->changePosition(spec.getAtomPosition(it->first));
    }
  }
}

bool PResidue::IsFullyDefined(const string &shellName, const PResidueSpec &spec)
{
  PBlockShell *shell = PResources::GetBlockShell(shellName);
  StringMap *blockAtoms = shell->getAtoms();

  for(StringMap::const_iterator it = blockAtoms->begin(); it != blockAtoms->end(); ++it) {
    if (!spec.Atom_Positions.contains(it->first)) {
      return false;
    }
  }
  return true;
}

void PResidue::inactivate() {
  vector<PAtom*> *atoms = getAtoms();
  vector<PAtom*>::iterator it;
  for (it=atoms->begin(); it!=atoms->end(); ++it) {
    (*it)->inactivate();
  }
}

void PResidue::activate() {
        vector<PAtom*> *atoms = getAtoms();
        vector<PAtom*>::iterator it;
        for (it=atoms->begin(); it!=atoms->end(); ++it) {
                (*it)->activate();
        }
}

// returns true if there is at least one atom in the residue is active
bool PResidue::isActive () {
	vector<PAtom*> *atoms = getAtoms();
	vector<PAtom*>::iterator it;
	for (it=atoms->begin(); it!=atoms->end(); ++it) {
		if ( (*it)->isActive() ) 
			return true;
	}
	return false;
}

void PResidue::printCollision() {
	vector<PAtom*> *atoms = getAtoms();
	vector<PAtom*>::iterator it;
	for (it=atoms->begin(); it!=atoms->end(); ++it) {
		(*it)->printCollision();
	}
}
