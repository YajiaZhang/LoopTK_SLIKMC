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

#ifndef __P_ATOM_H
#define __P_ATOM_H

#include <string>
using std::string;

#include <math3d/primitives.h>
#include <GLdraw/GLMaterial.h>
using namespace Math3D;

#include "PAtomShell.h"
#include "PFunctors.h"
#include "PHashing.h"
//#include "PBond.h"

class PBlock;
class PChain;
class PGrid;
class PSpaceManager;
class PResidue;

/**
 * Encapsulates a single atom in the protein loop. PAtom
 * contains methods to translate and apply matrix
 * transformations to the atom's position; find other
 * atom(s) in collision with it; and perform traversals
 * of the graph represented by atoms and their bonds. Makes
 * an assumption that the atoms of the same type are the same
 * size.
 */

class PAtom {
  public:
    friend class PBond;		// so it can add itself to the m_bonds array
    friend class PResidue;	// so that it can call CacheDOF
    friend class PBlock;		// so it can remove atoms on deactivation

    /**
     * Constructs a new PAtom in the specified block, with
     * the given name and global coordinates.
     */
    PAtom(PBlock * block, const string &name, Vector3 position, const string &id);

    /**
     * Destroys a PAtom, frees any memory associated
     * with it, and removes it from the collision grid.
     */
    ~PAtom();

    /**
     * Returns the ID this atom has within its residue.
     */
    string getID()const { return m_id; }

    /**
     * Moves this atom to the specified global coordinates
     * and updates its position in the collision grid.
     */
    void changePosition(const Vector3 &newPosition);

    //NOTE: Yajia added!
    void changePosition_nonGridUpdate( const Vector3 &newPosition);
    void updateGrid();

    /**
     * Applies the specified matrix transformation to
     * this atom's position.
     */
    void ApplyTransform(const Matrix4 &transform);

    /**
     * Returns the set of all atoms in static collision
     * with this atom.
     */
    AtomSet* getAllCollidingStatic() const;

    /**
     * Returns the set of all atoms in self collision
     * with this atom.
     */
    AtomSet* getAllCollidingSelf() const;

    /**
     * Returns the set of atoms in any collision (self
     * or static) with this atom.
     */
    AtomSet* getAllCollidingEither() const;

    /**
     * Returns bond to atom in residue res with id id.
    */
    PBond *getBond(PResidue *res, const string &id);
    
    /**
     * Returns any atom in static collision with this
     * atom, or NULL if no such atom exists.
     */
    PAtom *FindStaticCollision() const;

    /**
     * Returns any atom in self collision with this
     * atom, or NULL if no such atom exists.
     */
    PAtom *FindSelfCollision() const;

    /**
     * Returns any atom in either static or self collision
     * with this atom, or NULL if no such atom exists.
     */
    PAtom *FindAnyCollision() const;

    /**
     * Returns true if any atom is in static collision with
     * this atom, false otherwise.
     */
    bool InStaticCollision() const;

    /**
     * Returns true if any atom is in self collision with
     * this atom, false otherwise.
     */
    bool InSelfCollision() const;

    /**
     * Returns true if any atom is in either static or self
     * collision with this atom, false otherwise.
     */
    bool InAnyCollision() const;

    /**
     * Returns true if this atom is on the backbone
     */
    bool isOnBackbone() const;
      
    /**
     * Returns the atom's parent block.
     */
    PBlock *getParentBlock() { return m_atomBlock; }

    /**
     * Returns the atom's parent residue.
    */
    PResidue *getParentResidue();

    /**
     * Returns the atom's parent chain.
     */
    PChain *getChain() const;

    /**
     * Returns the atom's global coordinates.
     */
    Vector3 getPos() const { return m_atomPos; }

    /**
     * Returns all the bonds that include
     * this atom.
     */
    const vector<PBond *> *getBonds() const { return &m_bonds; }

    /**
     * Returns the atom's covalent radius.
     */
    Real getCovalentRadius() const { return m_atomShell->getCovalentRadius(); }

    /**
     * Returns the atom's van der Waals radius.
     */
    Real getVanDerWaalsRadius() const { return m_atomShell->getVanDerWaalsRadius(); }

    /**
     * Returns the atom's occupancy.
     */
    Real getOccupancy() const { return m_occupancy; }

    /**
     * Sets the atom's occupancy to <code>newVal</code>.
     */
    void setOccupancy(Real newVal) { m_occupancy = newVal; }

    /**
     * Returns the atom's temperature factor.
     */
    Real getTempFactor() const { return m_tempFactor; }

    /**
     * Sets the atom's temperature factor to <code>newVal</code>.
     */
    void setTempFactor(Real newVal) { m_tempFactor = newVal; }

    /**
     * Sets the color to be drawn when this atom
     * is displayed graphically. 
     * If not used, the default color is from PAtomShell.
     */
    void setColor(const GLColor setColor) {m_colorSet=true; m_atomColor = setColor;}

    /**
     * Reverts to the default color.
     */
    void revertColor(){ m_colorSet=false; }

    /**
     * Returns the color to be drawn when this atom
     * is displayed graphically.
     */
    GLColor getColor() const {
      if (m_colorSet) { 
        return m_atomColor;
      } else { 
        return m_atomShell->getColor();
      }
    }

    /**
     * Returns this atom's name.
     */
    string getName() const { return m_atomShell->getName(); }

    /**
     * Returns this atom's coordinates in the collision
     * grid.
     */
//    Vector3 getGridPos() const;

    /**
     * Returns true if this atom's parent block is
     * currently active, false otherwise.
     */
    bool WithinActiveBlock() const;

    /**
     * Returns true if this atom is bonded to <code>other</code>,
     * false otherwise.
     */
    bool isBonded(const PAtom *other) const;

    /**
     * Returns the length of the shortest bond path between <code>a1</code>
     * and <code>a2</code>.  If <code>a1 == a2</code>, returns 0; if
     * <code>a1</code> is bonded to <code>a2</code>, returns 1.
     *
     * The optional parameter <code>threshold</code> allows the caller to
     * specify a "cutoff" for the maximum length of the path.  If
     * <code>threshold >= 0</code> and the shortest path is found to be at
     * least <code>threshold</code> bonds long, this value will be returned
     * immediately.  The default value <code>threshold = -1</code> specifies
     * no cutoff.
     */
    static int shortestBondPath(const PAtom *a1, const PAtom *a2, int threshold = -1);

    /**
     * Initiates a traversal of the graph of atoms (nodes)
     * and their bonds (vertices) at this atom, using the
     * specified functors to be executed on each atom and bond.
     */
    void traverseChain(AtomFunctor *atomFn, BondFunctor *bondFn);

    /**
     * Initiates a traversal of the graph of atoms (nodes)
     * and their bonds (vertices) at this atom, using the
     * specified functor to be executed on each atom.
     */
    void traverseAtoms(AtomFunctor *atomFn) { traverseChain(atomFn, NULL); }

    /**
     * Initiates a traversal of the graph of atoms (nodes)
     * and their bonds (vertices) at this atom, using the
     * specified functor to be executed on each bond.
     */
    void traverseBonds(BondFunctor *bondFn) { traverseChain(NULL, bondFn); }

    /**
     * Returns the space manager for the chain this atom is in.
    */
    const PSpaceManager* getSpaceManager() const;

    /**
     * Returns true if this atom is active for collision checking,
     * false otherwise.
     */
    bool isActive() const { return m_active; }

    /**
     * Makes this atom inactive for collision checking.
     */
    void inactivate() { m_active = false; }

    /**
     * Make the atom active for collision checking.
     */
    void activate() { m_active = true; }

    /**
     * Print one atom which is in collision with this atom.
     */
    void printCollision();

    //NOTE: Yajia
    static AtomSet* internalTraverse( pair<PAtom*, PAtom*>& atom_pair, PChain* rootChain);

    Vector3 getGridPos();

  private:

    /* Helper functions. */
    void insertMeIntoGrid();
    void removeMeFromGrid();
    
    /* Called by PResidue. */
    void DestroyBonds();

    Vector3 m_atomPos;
    Vector3 m_gridPos;

    string m_id;
    PAtomShell *m_atomShell;
    vector<PBond *> m_bonds;
    PBlock *m_atomBlock;
    PGrid *m_grid;
    bool m_colorSet;
    GLColor m_atomColor;
    Real m_tempFactor, m_occupancy;
   
    void internalTraverse(AtomFunctor *atomFn,
                          BondFunctor *bondFn,
                          AtomSet &atomsTraversed,
                          PChain *rootChain,
                          PBond *bondFrom);
    //NOTE: Yajia
    void internalTraverse( AtomSet* atomsTraversed, PChain* rootChain);
    void internalTraverse_noGridUpdate( AtomFunctor *atomFn, AtomSet &atomsTraversed, PChain *rootChain);

    bool m_active;
};

#endif  // __P_ATOM_H
