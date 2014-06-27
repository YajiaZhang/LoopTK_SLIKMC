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
#include "PUtilities.h"
#include "PConstants.h"
#include "PBond.h"

#include <queue>
#include <math/math.h>
using namespace std;

#define    MIN_ALLOC    10000

/*
 * Constructors and destructors
 */

PAtom::PAtom(PBlock *block, const string &name, Vector3 position, const string &id) {
  m_atomBlock = block;
  m_atomPos = position;

  m_atomShell = PResources::GetAtomShell(name);
  m_colorSet = false;
  m_id = id;

  m_occupancy = 1.00;
  m_tempFactor = 0.00;

  PResidue *res = block->getParentResidue();
  PChain *chain = res->getChain();
  m_grid = chain->m_grid;

//  insertMeIntoGrid();
  this->m_grid->addAtom( this);

  m_gridPos = m_grid->scaleToGrid(m_atomPos);
  m_active = true;
}

PAtom::~PAtom() {
}

/*
 * Accessor methods
 */

PBond *PAtom::getBond(PResidue *res, const string &id) {
	for(unsigned i=0;i<m_bonds.size();i++) {
		PBond *b = m_bonds[i];
		PAtom *bondedAtom = (PAtom *) PUtilities::PointerThatIsNot(b->getAtom1(),b->getAtom2(),this);
		if(bondedAtom->getParentResidue()==res&&bondedAtom->getID()==id) return b;
	}
	return NULL;
}

//Vector3 PAtom::getGridPos() const {
//  return m_grid->scaleToGrid(m_atomPos);
//}

void PAtom::DestroyBonds() {
  for(unsigned i=0;i<m_bonds.size();i++) {
    PBond *b = m_bonds[i];
    PAtom *bondedAtom = (PAtom *) PUtilities::PointerThatIsNot(b->getAtom1(),b->getAtom2(),this);

    for(vector<PBond *>::iterator it=bondedAtom->m_bonds.begin();it!=bondedAtom->m_bonds.end();it++) {
      if (*it == b) {
        bondedAtom->m_bonds.erase(it);
        break;
      }
    }
    delete b;
  }
}


PChain* PAtom::getChain() const {
  if (m_atomBlock == NULL) {
    PUtilities::AbortProgram("Error: atom " + getName() + " is not encapsulated in a block!");
  }

  return m_atomBlock->getChain();
}

bool PAtom::WithinActiveBlock() const {
  if (m_atomBlock == NULL) {
    PUtilities::AbortProgram("Error: atom " + getName() + " is not encapsulated in a block!");
  }

  return m_atomBlock->isOn();
}

bool PAtom::isBonded(const PAtom *other) const {
  const vector<PBond *> *bonds = getBonds();
  PAtom *bondFirst, *bondSecond;

  for (vector<PBond *>::const_iterator it = bonds->begin(); it != bonds->end(); ++it) {
    bondFirst = (*it)->getAtom1();
    bondSecond = (*it)->getAtom2();

    if ((bondFirst == this && bondSecond == other) ||
        (bondFirst == other && bondSecond == this)) {
      return true;
    }
  }

  return false;
}

bool PAtom::isOnBackbone() const {
  return (m_atomBlock->getType() == "backbone");
}

/*
 * Collision detection methods
 */

bool PAtom::InStaticCollision() const {
  if (!m_active) return false;

  return m_grid->InStaticCollision(this);
}

bool PAtom::InSelfCollision() const {
  if (!m_active) return false;

  return m_grid->InSelfCollision(this);
}

bool PAtom::InAnyCollision() const {
  if (!m_active) return false;

  return m_grid->InAnyCollision(this);
}

AtomSet* PAtom::getAllCollidingStatic() const {
  if (!m_active) return new AtomSet();

  return m_grid->getAllCollidingStatic(this);
}

AtomSet* PAtom::getAllCollidingSelf() const {
  if (!m_active) return new AtomSet();

  return m_grid->getAllCollidingSelf(this);
}

AtomSet* PAtom::getAllCollidingEither() const {
  if (!m_active) return new AtomSet();

  return m_grid->getAllCollidingEither(this);
}

PAtom* PAtom::FindStaticCollision() const {
  if (!m_active) return NULL;

  return m_grid->getStaticCollidingAtom(this);
}

PAtom* PAtom::FindSelfCollision() const {
  if (!m_active) return NULL;

  return m_grid->getSelfCollidingAtom(this);
}

PAtom* PAtom::FindAnyCollision() const {
  if (!m_active) return NULL;

  return m_grid->getAnyCollidingAtom(this);
}

/*
 * Position and grid mutating methods
 */

void PAtom::changePosition(const Vector3 &newPosition) {
	if (WithinActiveBlock())
	{
		Vector3 curGridPos = getGridPos();
		Vector3 newGridPos = m_grid->scaleToGrid(newPosition);
		if (curGridPos != newGridPos)
		{
			removeMeFromGrid();
			m_atomPos = newPosition;
			insertMeIntoGrid();
		}
		else
		{
			m_atomPos = newPosition;
		}

//		Vector3 newGridPos = m_grid->scaleToGrid(newPosition);
//		if ( m_gridPos != newGridPos)
//		{
//			this->m_grid->removeAtom( this);
//			m_atomPos = newPosition;
//			this->m_grid->addAtom( this);
//			m_gridPos = newGridPos;
//		}
//		else
//		{
//			m_atomPos = newPosition;
//		}
	}
	else
	{
		m_atomPos = newPosition;
	}
}

void PAtom::changePosition_nonGridUpdate(const Vector3& newPosition) {
	m_atomPos = newPosition;
}

void PAtom::updateGrid()
{
	if( WithinActiveBlock())
	{
		Vector3 newGridPos = m_grid->scaleToGrid( m_atomPos);
		if( newGridPos != m_gridPos)
		{
			//Inconsistency happens, then change to new grid.
			m_grid->changeAtomPos( this, m_gridPos, newGridPos);
			m_gridPos = newGridPos;
		}
	}
}

Vector3 PAtom::getGridPos()
{
	Vector3 gridPos = m_grid->scaleToGrid( m_atomPos);
	return gridPos;
}

void PAtom::ApplyTransform(const Matrix4 &transform) {
  Vector4 myPos(getPos());
  myPos[3] = 1;
  Vector4 newPos4 = transform*myPos;
  Vector3 newPos3;
  newPos3[0] = newPos4[0]/newPos4[3];
  newPos3[1] = newPos4[1]/newPos4[3];
  newPos3[2] = newPos4[2]/newPos4[3];

  changePosition(newPos3);
}

void PAtom::insertMeIntoGrid() {
  m_grid->addAtom(this);
}

void PAtom::removeMeFromGrid() {
  m_grid->removeAtom(this);
}

/*
 * Atom-bond graph traversal methods
 */

int PAtom::shortestBondPath(const PAtom *a1, const PAtom *a2, int threshold) {
  queue<pair<const PAtom*, int> > bfsQueue;
  ConstAtomSet atomsSearched;

  pair<const PAtom*, int> curNode;
  const vector<PBond *> *curBonds;

  bfsQueue.push(make_pair(a1, 0));
  while(!bfsQueue.empty()) {
    curNode = bfsQueue.front();
    bfsQueue.pop();

    /* Did we find the target atom? */
    if (curNode.first == a2 || (threshold >= 0 && curNode.second >= threshold)) {
      return curNode.second;
    }

    /* Enqueue all atoms bonded to the current atom. */
    curBonds = curNode.first->getBonds();
    for(vector<PBond *>::const_iterator it = curBonds->begin(); it != curBonds->end(); ++it) {
      const PAtom *bondedAtom = (PAtom *)PUtilities::PointerThatIsNot((*it)->getAtom1(), (*it)->getAtom2(), curNode.first);
      if (atomsSearched.find(bondedAtom) == atomsSearched.end()) {
        atomsSearched.insert(bondedAtom);
        bfsQueue.push(make_pair(bondedAtom, curNode.second + 1));
      }
    }
  }

  return -1;
}

void PAtom::traverseChain(AtomFunctor *atomFn, BondFunctor *bondFn) {
  AtomSet traversed = AtomSet(MIN_ALLOC);
  internalTraverse(atomFn, bondFn, traversed, getChain(), NULL);
}

PResidue *PAtom::getParentResidue() { 
  return getParentBlock()->getParentResidue(); 
}

////NOTE: Yajia added
AtomSet* PAtom::internalTraverse(pair<PAtom*, PAtom*>& atom_pair,PChain* rootChain) {
	AtomSet* atomsTraversed = new AtomSet();
	atomsTraversed->insert( atom_pair.first);
	atom_pair.second->internalTraverse( atomsTraversed, rootChain);
	return atomsTraversed;
}


//NOTE: Yajia added
void PAtom::internalTraverse(AtomSet* atomsTraversed, PChain* rootChain) {
	if (WithinActiveBlock() && getChain()->IsSubChainOf(rootChain) && atomsTraversed->find(this) == atomsTraversed->end())
	{
		atomsTraversed->insert( this);
		for( int i = 0; i < m_bonds.size(); i++)
		{
			PAtom* atom1 = m_bonds[i]->getAtom1();
			PAtom* atom2 = m_bonds[i]->getAtom2();
			PAtom* opposingAtom = (atom1 != this ? atom1: atom2);
			opposingAtom->internalTraverse( atomsTraversed, rootChain);
		}
	}
}

void PAtom::internalTraverse_noGridUpdate(AtomFunctor *atomFn, AtomSet &atomsTraversed, PChain *rootChain)
{
	if (WithinActiveBlock() && getChain()->IsSubChainOf(rootChain) && atomsTraversed.find(this) == atomsTraversed.end())
	{
		((Rotater*)atomFn)->rotateAtom_nonGridUpdate( this);
		atomsTraversed.insert(this);

		for (unsigned i = 0; i < m_bonds.size(); i++) {
			PAtom* atom1 = m_bonds[i]->getAtom1();
			PAtom* atom2 = m_bonds[i]->getAtom2();
			PAtom* opposingAtom = (atom1 != this ? atom1: atom2);
			opposingAtom->internalTraverse_noGridUpdate(atomFn, atomsTraversed, rootChain);
		}
	}
}
//NOTE: This is how this method called
//nextAtom->internalTraverse(atomFn, bondFn, atomsTraversed, rootChain, this);
//bondFn: NULL;
//atomFn: Rotator
//
void PAtom::internalTraverse(AtomFunctor *atomFn, BondFunctor *bondFn, AtomSet &atomsTraversed, PChain *rootChain, PBond *bondFrom) {
	if (WithinActiveBlock() && getChain()->IsSubChainOf(rootChain) && atomsTraversed.find(this) == atomsTraversed.end())
	{
		if (atomFn != NULL)
			(*atomFn)(this, bondFrom);
		atomsTraversed.insert(this);

		//NOTE: Yajia changed here!
		if (bondFn != NULL)
		{
			for (unsigned i = 0; i < m_bonds.size(); i++)
			{
				PAtom* atom1 = m_bonds[i]->getAtom1();
				PAtom* atom2 = m_bonds[i]->getAtom2();
//				PAtom *opposingAtom = (PAtom *) PUtilities::PointerThatIsNot(m_bonds[i]->getAtom1(), m_bonds[i]->getAtom2(), this);
				PAtom* opposingAtom = (atom1 == this ? atom2: atom1);
				if (atomsTraversed.find(opposingAtom) == atomsTraversed.end() && opposingAtom->WithinActiveBlock())
				{
					(*bondFn)(m_bonds[i], this, opposingAtom);
				}
			}
		}

		for (unsigned i = 0; i < m_bonds.size(); i++) {
			PAtom* atom1 = m_bonds[i]->getAtom1();
			PAtom* atom2 = m_bonds[i]->getAtom2();
			PAtom* opposingAtom = (atom1 != this ? atom1: atom2);
			opposingAtom->internalTraverse(atomFn, bondFn, atomsTraversed, rootChain, m_bonds[i]);
		}
	}
}

const PSpaceManager* PAtom::getSpaceManager() const { 
  return getChain()->getSpaceManager(); 
}

void PAtom::printCollision() {
	PAtom* a2 = m_grid->getAnyCollidingAtom(this);
	if (a2!=NULL) {
		PResidue *r1 = this->getParentResidue();
		PResidue *r2 = a2->getParentResidue();

		Real radiusSum = this->getVanDerWaalsRadius() + a2->getVanDerWaalsRadius();
		Real centerDist = this->getPos().distance(a2->getPos());
		Real threshold = COLLISION_THRESHOLD * radiusSum;

		cout << "Colliding: " << r1->getPdbId() << " " << r1->getName() << " " << this->getID() << " with ";
		cout << r2->getPdbId() << " " << r2->getName() << " " << a2->getID();
		cout << " distance=" << centerDist << ", threshold=" << threshold << endl;
	}
}

