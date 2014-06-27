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

#ifndef __P_BLOCK_H
#define __P_BLOCK_H

#include "PBlockShell.h"
#include "PEnums.h"

class PBlockReconnector;

/**
 * Encapsulates a set of atoms and
 * bonds between them. Blocks may be
 * connected together to form residues.
 */

class PBlock {
 public:
  friend class PBlockConnection;  /* Needed so blocks can be connected together. */
  friend class PBlockReconnector; /* Needed for detached blocks to be reconnected. */
  
  /**
   * Constructs a new <code>PBlock</code> as a member of the residue
   * <code>myRes</code> and with a block shell specified by
   * <code>shellName</code>.  This constructor uses relative coordinates
   * as global positions.
   */
   
  PBlock(PResidue *myRes, const string &shellName);

  /**
   * Constructs a new <code>PBlock</code> as a member of the residue
   * <code>myRes</code> and with a block shell specified by
   * <code>shellName</code>.  This constructor uses the atom positions
   * specified by <code>spec</code> and should be called when the
   * <code>PBlock</code> is created "in isolation" (not connected to other blocks).
   */

  PBlock(PResidue *myRes, const string &shellName, const PAtomPositionsSpec &spec);

  /**
   * Constructs a new <code>PBlock</code> as a member of the residue
   * <code>myRes</code> and with a block shell specified by
   * <code>shellName</code>.  This constructor uses the atom positions
   * specified by <code>spec</code> and connects the new <code>PBlock</code>
   * to the block <code>toConnect</code>; it should be called when the
   * block's positions are known, e.g. from a PDB file, but the block
   * must also be connected to an existing block.
   */

  PBlock(PResidue *myRes, const string &shellName, const PAtomPositionsSpec &spec, PBlock *toConnect);

  /**
   * Constructs a new <code>PBlock</code> as a member of the residue
   * <code>myRes</code> and with a block shell specified by
   * <code>shellName</code>.  This constructor uses the block's 
   * default atom positions and connects the new <code>PBlock</code>
   * to the block <code>toConnect</code>; it should be called when the
   * block's positions are not known (and thus should be determined
   * randomly) but the block must still be connected to an existing block.
   */

  PBlock(PResidue *myRes, const string &shellName, PBlock *toConnect);

  ~PBlock();

  /**
   * Returns the block's name, e.g. "TRP".
   */

  string getName() const {return m_shell->getName(); }

  /**
   * Returns the block's type, e.g. "backbone" or "sidechain".
   */

  string getType() const {return m_shell->getType(); }

  /**
   * Returns <code>true</code> if this block is active,
   * <code>false</code> otherwise.
   */

  bool isOn() const { return m_isOn; }

  /**
   * Returns this block's parent residue.
   */

  PResidue *getParentResidue() { return m_residue; }

  /**
   * Returns the number of atoms within this block.
  */

  int size() { return m_atoms.size(); }

  /**
   * Returns this block's parent chain.
   */

  PChain* getChain() const;
  
  /**
   * Detaches this block from the parent <code>PChain</code>. The
   * <code>PBlock</code> "remembers" the relative coordinates to
   * the first block it can find of type <code>blockType</code>,
   * so the <code>Detach</code> method can only be called once.
   */

  void Detach(const string &blockType);

  /**
   * Reattaches this block to its parent chain.
   */

  void Reattach();

  /**
   * Returns true if any atom in this <code>PBlock</code>
   * is in static collision, false otherwise.
   */

  bool InStaticCollision() const;

  /**
   * Returns true if any atom in this <code>PBlock</code>
   * is in self collision, false otherwise.
   */

  bool InSelfCollision() const;

  /**
   * Returns true if any atom in this <code>PBlock</code>
   * is in static or self collision, false otherwise.
   */

  bool InAnyCollision() const;
  
  /**
   * Returns a pair of atoms in static collision within
   * this <code>PBlock</code>, or <code>(NULL, NULL)</code>
   * if no such atoms exist.
   */

  pair<PAtom *, PAtom *> FindStaticCollision() const;

  /**
   * Returns a pair of atoms in self collision within
   * this <code>PBlock</code>, or <code>(NULL, NULL)</code>
   * if no such atoms exist.
   */

  pair<PAtom *, PAtom *> FindSelfCollision() const;

  /**
   * Returns a pair of atoms in either static or self collision
   * within this <code>PBlock</code>, or <code>(NULL, NULL)</code>
   * if no such atoms exist.
   */

  pair<PAtom *, PAtom *> FindAnyCollision() const;

  /**
   * Returns a set of all atoms in this <code>PBlock</code>
   * in static collision.
   */

  AtomCollisions* getAllCollidingStatic() const;

  /**
   * Returns a set of all atoms in this <code>PBlock</code>
   * in self collision.
   */

  AtomCollisions* getAllCollidingSelf() const;

  /**
   * Returns a set of all atoms in this <code>PBlock</code>
   * in either static or self collision.
   */

  AtomCollisions* getAllCollidingEither() const;

  /**
   * Adds all atoms in this <code>PBlock</code> to the map
   * <code>atomCache</code>.
   */

  void AddAtomsToMap(HASH_MAP_STR(PAtom *) &atomCache) const;

 /**
   * Returns the space manager for the chain this block is in.
  */

  const PSpaceManager* getSpaceManager() const;
 private:

  void CreateBlock(PResidue *myRes, const PAtomPositionsSpec &spec);
  void ActivateAndDestroyReattachment();
  void ApplyTransform(Matrix4 transform);


  PBlockShell *m_shell;
  PResidue *m_residue; 
  HASH_MAP_STR(PAtom *) m_atoms; //atom Id's to the respective atom

  /* Collision detection helper methods. */

  bool InCollision(CollisionType type) const;
  pair<PAtom *, PAtom *> FindCollision(CollisionType type) const;
  AtomCollisions* getAllColliding(CollisionType type) const;

  //prepared by finalize()
  vector<PBond *> m_dofs;


  //selective manipulation data members:
  bool m_isOn;
  vector<PBlock *> m_connectedBlocks;
  PBlockReconnector *m_blockReconnector;

  
  //selective manipulation helper functions
  PBlock *FindFirstBlockOfType(const string &blockType);
  void TraverseBlocks(PBlock *from, void (*ManipulateBlockFunc)(PBlock *block, PBlock *from), bool skipOffBlocks);
  static void DeactivateBlock(PBlock *block, PBlock *from);
  static void ActivateBlock(PBlock *block, PBlock *from);
  static void ReattachBlock(PBlock *block, PBlock *from);
};

#endif  // __P_BLOCK_H
