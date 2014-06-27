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
#include "PExtension.h"
#include "PResources.h"
#include "PTools.h"

using namespace std;

Rotater::Rotater(Vector3 orig, Vector3 axis, float degrees) {
    rotMat = PMath::FindRotationMatrix(axis, DtoR(degrees));
    origin = orig;
  }

void Rotater::operator()(PAtom *atom, PBond *bondFrom)
{
    Vector3 currPos = atom->getPos()-origin;
    atom->changePosition(rotMat*currPos+origin);
  }

void Rotater::rotateAtom_nonGridUpdate(PAtom* atom) {
	Vector3 currPos = atom->getPos()-origin;
	atom->changePosition_nonGridUpdate(rotMat*currPos+origin);
}


//struct Rotater: AtomFunctor {
//public:
//  Rotater(Vector3 orig, Vector3 axis, float degrees) {
//    rotMat = PMath::FindRotationMatrix(axis, DtoR(degrees));
//    origin = orig;
//  }
//
//  void operator()(PAtom *atom, PBond *bondFrom) {
//    Vector3 currPos = atom->getPos()-origin;
//    atom->changePosition(rotMat*currPos+origin);
//  }
//
//private:
//  Matrix3 rotMat;
//  Vector3 origin;
//};

PAtom *PBond::getAtom(BondDirection dir) const {
  if (dir == forward) {
    return m_forwardDirection;
  } else {
    return (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
  }
}
//NOTE: How this called!
//traverseChain(dir,&r,NULL,chain);
//r: Rotator: AtomFunctor, BondFunctor* bondFn remains NULL.
void PBond::traverseChain(BondDirection dir, AtomFunctor *atomFn, BondFunctor *bondFn, PChain *rootChain)
{
  PAtom *prevAtom, *nextAtom;
  AtomSet atomsTraversed;

  //NOTE: Yajia changed here! Rewrite the code without calling the function PUtilities::PointerThatIsNot()
  if (dir == forward) {
//    prevAtom = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
	prevAtom = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
    nextAtom = m_forwardDirection;
  }
  else {
    prevAtom = m_forwardDirection;
//    nextAtom = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
    nextAtom = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
  }

  if (bondFn != NULL) {
    (*bondFn)(this, prevAtom, nextAtom);
  }

  atomsTraversed.insert(prevAtom);
  nextAtom->internalTraverse(atomFn, bondFn, atomsTraversed, rootChain, this);
}


void PBond::traverseChain_noGridUpdate(BondDirection dir, Rotater* atomFn, BondFunctor* bondFn, PChain* rootChain) {
	  PAtom *prevAtom, *nextAtom;
	  AtomSet atomsTraversed;

	  //NOTE: Yajia changed here! Rewrite the code without calling the function PUtilities::PointerThatIsNot()
	  if (dir == forward) {
	//    prevAtom = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
		prevAtom = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	    nextAtom = m_forwardDirection;
	  }
	  else {
	    prevAtom = m_forwardDirection;
	//    nextAtom = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
	    nextAtom = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	  }

	  if (bondFn != NULL) {
	    (*bondFn)(this, prevAtom, nextAtom);
	  }

	  atomsTraversed.insert(prevAtom);
	  nextAtom->internalTraverse_noGridUpdate(atomFn, atomsTraversed, rootChain);
}

pair<PAtom*, PAtom*> PBond::getAtomPair(BondDirection dir) {
	  assert(isDOF());
	  PAtom *prevAtom, *nextAtom;

	  //NOTE: Yajia changed here! Rewrite the code without calling the function PUtilities::PointerThatIsNot()
	  if (dir == forward) {
	//    prevAtom = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
		prevAtom = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	    nextAtom = m_forwardDirection;
	  }
	  else {
	    prevAtom = m_forwardDirection;
	//    nextAtom = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
	    nextAtom = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	  }
	  return make_pair<PAtom*, PAtom*>(prevAtom, nextAtom);
}

void PBond::Rotate(BondDirection dir, float degrees) {
	PChain *c1 = getAtom1()->getChain();
	PChain *c2 = getAtom2()->getChain();
	PChain *common = PTools::LowestCommonChain(c1,c2);
	Rotate(dir,degrees,common);
}

void PBond::Rotate(BondDirection dir, float degrees, PChain *chain) {
	  assert(isDOF());
	  PAtom *startV;
	  PAtom *endV;

	  if (dir==forward) {
		endV = m_forwardDirection;
		//NOTE: Yajia changed here! Rewrite the code without calling the function PUtilities::PointerThatIsNot()
	//    startV = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
		startV = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	  }
	  else {
		startV = m_forwardDirection;
	//    endV = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
		endV = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	  }

	  Vector3 axis = endV->getPos()-startV->getPos();
	  Rotater r(startV->getPos(),axis,degrees);
	  traverseChain(dir,&r,NULL,chain);
}

void PBond::Rotate_noGridUpdate(BondDirection dir, float degrees, PChain *chain)
{
	  assert(isDOF());
	  PAtom *startV;
	  PAtom *endV;

	  if (dir==forward) {
		endV = m_forwardDirection;
		//NOTE: Yajia changed here! Rewrite the code without calling the function PUtilities::PointerThatIsNot()
	//    startV = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
		startV = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	  }
	  else {
		startV = m_forwardDirection;
	//    endV = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
		endV = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	  }

	  Vector3 axis = endV->getPos()-startV->getPos();
	  Rotater r(startV->getPos(),axis,degrees);
//	  traverseChain(dir,&r,NULL,chain);
	  this->traverseChain_noGridUpdate( dir, &r, NULL, chain);
	  return;
}

Rotater* PBond::getBondRotater(BondDirection dir, float degrees) {
	  assert(isDOF());
	  PAtom *startV;
	  PAtom *endV;

	  if (dir==forward) {
	    endV = m_forwardDirection;
	    //NOTE: Yajia changed here! Rewrite the code without calling the function PUtilities::PointerThatIsNot()
	//    startV = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
	    startV = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	  }
	  else {
	    startV = m_forwardDirection;
	//    endV = (PAtom *) PUtilities::PointerThatIsNot(m_atom1, m_atom2, m_forwardDirection);
	    endV = (m_atom1 == m_forwardDirection ? m_atom2: m_atom1);
	  }

	  Vector3 axis = endV->getPos()-startV->getPos();
	  Rotater* r = new Rotater(startV->getPos(),axis,degrees);
	  return r;
}

Real PBond::getLength() {
  return m_atom1->getPos().distance(m_atom2->getPos());
}





