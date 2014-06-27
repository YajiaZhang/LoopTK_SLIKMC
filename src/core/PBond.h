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

#ifndef __P_BOND_H
#define __P_BOND_H

#include "PEnums.h"

struct Rotater: AtomFunctor {
public:
  Rotater(Vector3 orig, Vector3 axis, float degrees);

  void operator()(PAtom *atom, PBond *bondFrom);

  void rotateAtom_nonGridUpdate( PAtom* atom);

private:
  Matrix3 rotMat;
  Vector3 origin;
};

/**
 * Encapsulates a bond between two atoms. Contains
 * pointers to the two bonded atoms and a binary
 * bit indicating whether the bond is a degree of
 * freedom.
 */
class PBond {
 public:
  friend struct DOFCacher;

  /**
   * Constructs a new PBond.
   */
  PBond(PAtom *a1, PAtom *a2, bool isDOF) { 
    m_atom1 = a1; 
    m_atom2 = a2;
    m_isDOF = isDOF;
    m_forwardDirection = NULL;
    a1->m_bonds.push_back(this);
    a2->m_bonds.push_back(this);
 }


  /**
   * Returns the length of this bond.
  */
  Real getLength();

  /**
   * Initiates a traversal of the graph of atoms and bonds
   * in the specified direction, calling the given atom functor
   * at each atom and the given bond functor at each bond.
   */
  void traverseChain(BondDirection dir, AtomFunctor *atomFn, BondFunctor *bondFn, PChain *rootChain);

  //Yajia Zhang added!
  void traverseChain_noGridUpdate(BondDirection dir, Rotater *atomFn, BondFunctor *bondFn, PChain *rootChain);
  void Rotate_noGridUpdate(BondDirection dir, float degrees, PChain *chain);
//  void traverseChain_side(BondDirection dir, AtomFunctor *atomFn, BondFunctor *bondFn, PChain *rootChain);

  /**
   * Initiates a traversal of the graph of atoms and bonds
   * in the specified direction, calling the given atom functor
   * at each atom.
   */
  void traverseAtoms(BondDirection dir, AtomFunctor *atomFn, PChain *rootChain) {
    traverseChain(dir, atomFn, NULL, rootChain);
  }

  /**
   * Initiates a traversal of the graph of atoms and bonds
   * in the specified direction, calling the given bond functor
   * at each bond.
   */
  void traverseBonds(BondDirection dir, BondFunctor *bondFn, PChain *rootChain) {
    traverseChain(dir, NULL, bondFn, rootChain);
  }

  /**
   * Rotates the chain starting from this bond by
   * the specified number of degrees, moving in
   * the direction indicated by dir, limiting the rotation to "chain".
   */
  void Rotate(BondDirection dir, float degrees, PChain *chain);

  /**
   * Does a rotation on the smallest common subchain between the two atoms.
   */ 
  void Rotate(BondDirection dir, float degrees);

  /**
   * Returns the first atom in the bond.
   */
  PAtom *getAtom1() const { return m_atom1; }

  /**
   * Returns the second atom in the bond.
   */
  PAtom *getAtom2() const { return m_atom2; }

  /**
   * Returns the bonded atom in the specified
   * direction.
   */
  PAtom *getAtom(BondDirection dir) const;

  /**
   * Returns true if this bond is a degree
   * of freedom, false otherwise.
   */
  bool isDOF() const { return m_isDOF; }

  //Yajia Zhang added this:
  Rotater* getBondRotater( BondDirection dir, float degrees);
  pair< PAtom*, PAtom*> getAtomPair( BondDirection dir);
 private:

  /*
   * Sets the specified atom to be the
   * "forward" atom of this bond. (for DOFCacher)
   */
  void setDirection(PAtom *front) { m_forwardDirection = front; }

  PAtom *m_atom1;
  PAtom *m_atom2;
  PAtom *m_forwardDirection;
  bool m_isDOF;
};

#endif  // __P_BOND_H
