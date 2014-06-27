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

#ifndef __P_CHAIN_H
#define __P_CHAIN_H

#include "PAtom.h"
#include "PAtomNode.h"
#include "PGrid.h"
#include "PLightChain.h"
#include "PResidue.h"
#include "PStructs.h"
#include "PChainState.h"
#include "DihedralAngle.h"

#include <list>
using std::list;

//@package Main Infrastructure
/**
 *
 * Encapsulates a linearly connected chain of residues.
 * Methods are provided to add and access residues and check for collisions
 * within the chain. This class makes no assumptions 
 * on the structure of the residues and atoms within.
 * Use PProtein for extended functionality if building chains 
 * using the configuration data provided with %LoopTK.
 * 
 *  Usage of a PChain must follow the following sequence of steps:
 *
 * -# \b Construction - all residues are added to chain. Use the
 *       PChain::AddResidue() methods to do this.
 * -# \b Finalization - finishes construction of the chain and prepares it
 *       for manipulation. Use PChain::finalize() to do this.
 * -# \b Analysis - the full range of functionality of PChain can be used
 *       at this point. However, no more residues can be added.
 *
 */
class PChain: public PLightChain {
 public:
  friend class PAtom; // so that it can access the grid

  /**
   * Constructs a new \c PChain with no residues.
   */
  PChain();

  /**
   * Constructs a new \c PChain with one residue, using
   * the default positions for its atoms.
   */
  PChain(const string &firstResidueName);

  /**
   * Constructs a new \c PChain with one residue (\c firstResidueName), using
   * the positions specified by \c firstResidueSpec for its atoms.
   */
  PChain(const string &firstResidueName, PResidueSpec &firstResidueSpec);

  /**
   * Constructs a new \c PChain as a subchain of \c protein,
   * from residue \c resStartIndex to residue \c resEndIndex.
   */
  PChain(PChain *protein, int resStartIndex, int resEndIndex);

  /**
   * Destroys this \c PChain and its subchains, freeing
   * any memory associated with them.
   */
  ~PChain();

  /**
   * Obliterates this chain and all parents, brethren, and children, 
   * freeing all memory associated with them.
   */
  void Obliterate();


  /**
   * Adds the residue \c resName to the end
   * of the chain, using its default atom
   * positions. Returns the created residue.
   */
  PResidue *AddResidue(const string &resName);

  /**
   * Adds the specified residue to the end
   * of the chain, using the specified atom
   * positions. Returns the created residue.
   */
  PResidue *AddResidue(const string &resName, PResidueSpec &resSpec);

  /**
   * Returns true if this \c PChain is a subchain
   * of \c other, false otherwise.
   */
  bool IsSubChainOf(const PChain *other) const;

  /**
   * Gets parent chain. Returns \c NULL if this is the
   * top level chain.
   */
  PChain *getParent() const;

  /**
   * Gets top level chain (the one for which
   * <tt>parent == NULL</tt>).
   */
  PChain *getTopLevelChain();

  /**
   * Returns the bond between atoms \c id1 and \c id2, assuming both
   * are in residue \c res_index.
   *
   * Returns \c NULL if the bond cannot be found.
   */
  PBond *getBond(int res_index, const string &id1, const string &id2);

  /**
   * Returns the bond between two atoms in different residues: the first has ID
   * \c id1 in residue \c res_index1, the second has ID \c id2 in residue
   * \c res_index2.
   *
   * Returns \c NULL if the bond cannot be found.
   */
  PBond *getBond(int res_index1, const string &id1, int res_index2, const string &id2);

  /**
   * Returns the number of residues in this \c PChain.
   */
  int size() const;

  /**
   * Returns the residue at the specified index. The
   * index is "local" to this PChain; for example, if
   * this chain is a subchain of some other chain,
   * accessing residue 0 in this chain might return
   * the 20th residue in the parent chain.
   */
  PResidue *getResidue(int localIndex);

  /**
   * Returns the residue index of the specified atom.
   */
  int getResidueIndex(const PAtom* atom);


  /**
   * Convenience method: returns the atom with the
   * specified ID in the specified residue, or NULL
   * if such an atom does not exist.
   */
  PAtom* getAtomAtRes(const string &atomID, int resNum);

  /**
   * Sets color of the atom (displayed by PChainNavigator) whose ID
   * is \c atomID in residue \c resNum.
   */
  void setAtomColorAtRes(const string &atomID, int resNum, GLColor setColor) {
    getAtomAtRes(atomID,resNum)->setColor(setColor);
  }

  /**
   * Reverts the color of the atom whose ID is \c atomID,
   * in residue \c resNum, to its default value.
   */
  void revertAtomColorAtRes(const string &atomID, int resNum) {
    getAtomAtRes(atomID,resNum)->revertColor();
  }

  /**
   * Retrieves the atom at the given index for the block type
   * (indices are in order of how the chain would be traversed
   * from beginning to end).
   */
  PAtom *getAtom(const string &blockType, int index);
  
  /**
   * Gets the atom position at the given index (see PChain::getAtom).
   */
  Vector3 getAtomPos(const string &blockType,int index);

  /**
   * Returns the number of atoms in this chain of the given \c blockType.
   */
  int NumAtoms(const string &blockType) const;

  /**
   * Returns \c true if any atom in this \c PChain is in
   * static collision, \c false otherwise.
   */
  bool InStaticCollision();

  /**
   * Returns \c true if any atom in this \c PChain between
   * the specified residue indices is in static
   * collision, \c false otherwise.
   */
  bool InStaticCollision(int resIndex1, int resIndex2);

  /**
   * Returns \c true if any atom in this \c PChain is in
   * self collision, \c false otherwise.
   */
  bool InSelfCollision();

  /**
   * Returns \c true if any atom in this \c PChain between
   * the specified residue indices is in self
   * collision, \c false otherwise.
   */
  bool InSelfCollision(int resIndex1, int resIndex2);

  /**
   * Returns \c true if any atom in this PChain is in
   * static or self collision, \c false otherwise.
   */
  bool InAnyCollision();

  /**
   * Returns \c true if any atom in this \c PChain between
   * the specified residue indices is in static or
   * self collision, \c false otherwise.
   */
  bool InAnyCollision(int resIndex1, int resIndex2);

  /**
   * Returns a pair of atoms in static collision within
   * the chain, or <tt>(NULL, NULL)</tt> if no
   * such atoms exist.
   */
  pair<PAtom *, PAtom *> FindStaticCollision();

  /**
   * Returns a pair of atoms in static collision between
   * the specified residue indices, or <tt>(NULL, NULL)</tt>
   * if no such atoms exist.
   */
  pair<PAtom *, PAtom *> FindStaticCollision(int resIndex1, int resIndex2);

  /**
   * Returns a pair of atoms in self collision within
   * the chain, or <tt>(NULL, NULL)</tt> if no such
   * atoms exist.
   */
  pair<PAtom *, PAtom *> FindSelfCollision();

  /**
   * Returns a pair of atoms in self collision between
   * the specified residue indices, or <tt>(NULL, NULL)</tt>
   * if no such atoms exist.
   */
  pair<PAtom *, PAtom *> FindSelfCollision(int resIndex1, int resIndex2);

  /**
   * Returns a pair of atoms in either static or self collision
   * in the chain, or <code>(NULL, NULL)</code> if no such
   * atoms exist.
   */
  pair<PAtom *, PAtom *> FindAnyCollision();

  /**
   * Returns a pair of atoms in either static or self collision
   * between the specified residue indices, or <code>(NULL, NULL)</code>
   * if no such atoms exist.
   */
  pair<PAtom *, PAtom *> FindAnyCollision(int resIndex1, int resIndex2);

  /**
   * Returns a set of all atoms in the chain in
   * static collision.
   */
  AtomCollisions* getAllCollidingStatic();

  /**
   * Returns a set of all atoms between the specified
   * residue indices in static collision.
   */
  AtomCollisions* getAllCollidingStatic(int resIndex1, int resIndex2);

  /**
   * Returns a set of all atoms in the chain in
   * self collision.
   */
  AtomCollisions* getAllCollidingSelf();

  /**
   * Returns a set of all atoms between the specified
   * residue indices in self collision.
   */
  AtomCollisions* getAllCollidingSelf(int resIndex1, int resIndex2);

  /**
   * Returns a set of all atoms in the chain in
   * either static or self collision.
   */
  AtomCollisions* getAllCollidingEither();

  /**
   * Returns a set of all atoms between the specified
   * residue indices in either static or self collision.
   */
  AtomCollisions* getAllCollidingEither(int resIndex1, int resIndex2);

  /**
   * Returns all the degrees of freedom for the specified
   * type of block (such as backbone or sidechain) in the
   * PChain.
   */
  vector<PBond *>& GetDOFs(const string &blockType);
  
  /**
   * Returns the number of degrees of freedom (DOF) of the specified block type 
   * within this chain.
   */
  int NumDOF(const string &blockType) const;
  
  /**
   * Rotates the chain for a given blocktype, degrees of freedom (DOF), direction,
   * and angle in degrees.
   */
  void RotateChain(const string &blockType,int DOFindex, BondDirection dir, float degrees);

  /**
   * @brief (Yajia Zhang added) Do similar work as RotateChain but only virtually change the atom position without changing the atom grid map.
   */
  void RotateChain_noGridUpdate(const string &blockType, int DOFindex, BondDirection dir, float degrees);

  /**
   * Rotates the chain given a <code>ChainMove</code>.
   */
  void RotateChain(const ChainMove &move);

  /**
   * @brief (Yajia Zhang added) Do similar work as RotateChain but only virtually change the atom position without changing the atom grid map.
   */
  void RotateChain_noGridUpdate(const ChainMove &move);

  /**
   * MultiRotates the chain given a vector of <code>ChainMove</code>.
   */
  void MultiRotate(vector<ChainMove> &moves);

  /**
   * @brief (Yajia Zhang added) Do similar work as MultiRotate but only virtually change the atom position without changing the atom grid map.
   * @param moves
   */
  void MultiRotate_noGridUpdate(vector<ChainMove> &moves);
  /**
   * Does the opposite rotations in the opposite order as MultiRotate() would do.
   */  
  void AntiMultiRotate(vector<ChainMove> &moves);
  
  /**
   * Returns all the block types in this <code>PChain</code>.
   */
  vector<string> GetBlockTypes() const;

  /**
   * Detaches all \c blockType blocks attached to \c blockTypeToDetachFrom.
   */
  void DetachBlocks(const string &blockType, const string &blockTypeToDetachFrom);

  /**
   * Detaches \c blockType blocks attached to \c blockTypeToDetachFrom 
   * from \c resStartIndex to \c resEndIndex.
   */
  void DetachBlocks(const string &blockType, const string &blockTypeToDetachFrom, int resStartIndex, int resEndIndex);

  /**
   * Reattaches blocks of type \c blockType to the original position.
   */
  void ReattachBlocks(const string &blockType);

  /**
   * Returns blocks of type \c blockType and within the range of
   * \c resStartIndex to \c resEndIndex to their original positions.
   */
  void ReattachBlocks(const string &blockType, int resStartIndex, int resEndIndex);

  /**
   * Reattaches all detached blocks regardless of block type.
   */
  void ReattachAllBlocks();  

  /**
   * Reattaches all detached blocks regardless of blocktype from resStartIndex 
   * to resEndIndex.
   */
  void ReattachAllBlocks(int resStartIndex, int resEndIndex);

  /**
   * Applies a random rotation in the given direction 
   * to all degrees of freedom in the chain (or subchain).
   */
  void RandomizeDOFs(BondDirection dir);

  /**
   * Finalizes the chain and makes it ready for manipulation.
   * Must be called as soon as the chain is finished being 
   * constructed or it will be unstable.
   */
  void finalize();

  /**
   * Adds the object to the event list for rotations
   * (HandleRotation function invoked whenever there is a rotation)
   */
  void AddRotateEventHandler(PRotateEventHandler *handler);

  /**
   * Removes the given event handler from the rotation list.
   * Aborts program if it's not in the list
   */
  void RemoveRotateEventHandler(PRotateEventHandler *handler);


  /**
   * Traverses the chain from the start atom of residue 0, calling the
   * specified <code>AtomFunctor</code> at each atom.
   */
  void traverseFromStart(AtomFunctor *atomFn) { traverseFromStart(atomFn, NULL); }

  /**
   * Traverses the chain from the start atom of residue 0, calling the
   * specified <code>BondFunctor</code> at each bond.
   */
  void traverseFromStart(BondFunctor *bondFn) { traverseFromStart(NULL, bondFn); }

  /**
   * Traverses the chain from the start atom of residue 0, calling the
   * specified <code>AtomFunctor</code> at each atom <i>and</i> the
   * specified <code>BondFunctor</code> at each bond.
   */
  void traverseFromStart(AtomFunctor *atomFn, BondFunctor *bondFn);

  /**
   * Returns a <code>pair</code> of integers indicating
   * the residues in the top-level chain represented by
   * this <code>PChain</code>.  The first element of the
   * <code>pair</code> is the start residue index; the
   * second element of the <code>pair</code> is the end index.
   */
  pair<int, int> getTopLevelIndices() const;
  
  /**
   * Returns the space manager of this <code>PChain</code>.
   */
  const PSpaceManager* getSpaceManager() { return m_grid; }

  /**
   * Creates an exact copy of this <code>PChain</code>. If
   * a subchain is cloned, all its parents will be replicated
   * as well; thus, to fully deallocate the newly-cloned
   * chain's memory, invoke <code>Obliterate()</code> when
   * it is no longer needed.
   */
  PChain *Clone();

  /**
   * Return the local index of a residue in the chain. 
   * Return -1 if there is no such residue in the chain.
   */
  int getResidueLocalIndex(PResidue *res);

  /**
   * Return the global index of a residue in the top level chain.
   * Return -1 if there is no such residue.
  */
  int getResidueGlobalIndex(PResidue *res);

 /**
  * Returns a <code>AtomPath</code> which is a vector going through
  * a set of atoms. This path is generated from an atom to a set of 
  * atoms and the returning path will be the shortest one.
  */
  vector<const PAtom*> getShortestPath(const PAtom *a1, ConstAtomSet a2Set, int maxLength);
 
  /**
   * Returns a <code>AtomPath</code> which is a vector going through
   * a set of atoms. This path is generated from an atom to another
   * atom.
   */
  vector<const PAtom*> getShortestPath(const PAtom *a1, const PAtom *a2, int maxLength);    
   
  /**
   * Returns a <code>AtomPath</code> which is a vector going through
   * a set of atoms. This path is generated from an atom to another
   * atom with some additional atoms in two end. 
   */
  vector<const PAtom*> getShortestPath(const PAtom *a1, ConstAtomSet a2Set, int maxLength, int extendNum);

  /**
   * Inactivates the residues from index \c Rid1 to \c Rid2 so that these
   * residues are not considered in collision detection.
   */
  void inactivateResidue (int Rid1, int Rid2);

  /**
   * Activates the residues from index \c Rid1 to \c Rid2 so that these residues
   * are considered in collision
   */
  void activateResidue (int Rid1, int Rid2);

  /**
   * Returns number of children to this PChain.
   */
  unsigned int NumChildren() { return m_children.size(); }

  void updateAtomsGrid();

  /**
   * Returns child PChain to this PChain at the given index.
   */
  PChain *getChild(unsigned int index) { 
   for(list<PChain *>::iterator it = m_children.begin();it!=m_children.end();it++) {
    if(index==0) return *it;
    index--;
   }
   return NULL;
  }

  /*NOTE: Start my code here!*/
  PChainState* savePartialChainState( int index_start, int index_end);
  PChainState* saveChainState();
  void restoreChainState( PChainState* state);
  void restoreChainState_noGridUpdate( PChainState* state);

  //NOTE: Gain control to all the residues in the chain
  void attachResidues();
  //NOTE: Gain control to the residues with local index start to end in the chain
  void attachResidues( int start, int end);
  void attachResidue( int index);
  //NOTE: Work with attachResidues, lose control to its residues.
  //Parent chain regains the control.
  void detachResidues();
  void detachResidues( int start, int end);
  void detachResidue( int index);
  bool areResiduesFullyControl();

  void getBackbonePositions( vector<Vector3>& positions);
  DihedralAngle* getDihedralAngleAtResidue( int index);
  void getDihedralAngles( vector<DihedralAngle>& angles);
  /*NOTE: End my code here!*/

  void printDihedralAngles();
  void printBackbonePosition();
  void printBackboneCollisionGrid();
  vector<Vector3>* getBackboneCollisionGrid();
 protected:
  void CloneResiduesIntoChain(PChain *other);
  virtual PResidue *CreateResidue(const string &name) { return new PResidue(this,name); }
  virtual PResidue *CreateResidue(const string &name, PResidueSpec &spec) { return new PResidue(this,name,spec); }
  virtual PResidue *CreateResidue(const string &name, PResidue *res) { return new PResidue(this,name,res); }
  virtual PResidue *CreateResidue(const string &name, PResidueSpec &spec, PResidue *res)
  { return new PResidue(this,name,spec,res); }

 private:
  void updateDihedralAngles_subchain( PChain* child);
  void UpdateIndexRangeOnAdd(int amtAdded);
  void InitChain();
  void CheckFinalized() const;
  vector<const PAtom*> extractPath(AtomNode* leaveNode); //post-process of getShortestPath


  bool m_isFinalized;
  vector<PResidue *> *m_residues;
  PGrid *m_grid;

  int m_startIndex;		/* Start index into the parent chain, if applicable. */
  int m_endIndex;		/* End index into the parent chain, if applicable. */

  PChain *m_parentChain;	/* Parent chain, or NULL if this is the whole chain. */
  list<PChain *> m_children;	/* List of children of this PChain. */

  //QUESTION: What's this for?
  list<PRotateEventHandler *> *m_rotationEvents;

  /* A bond's block type is the block type of the atom in the forward direction. */

  DOF_Cache m_dofs;		/* Map of block types to DOF's of that type. */
  Atom_Cache m_atomCache;

  /* Collision-detection helper methods. */

  bool InCollision(int resIndex1, int resIndex2, CollisionType type);
  pair<PAtom *, PAtom *> FindCollision(int resIndex1, int resIndex2, CollisionType type);
  AtomCollisions* getAllColliding(int resIndex1, int resIndex2, CollisionType type);
};

#endif  // __P_CHAIN_H
