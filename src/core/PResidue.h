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

#ifndef __P_RESIDUE_H
#define __P_RESIDUE_H

#include "PEnums.h"
#include "PResidueShell.h"
#include "PResidueSpec.h"

/**
 * Encapsulates a single residue in a chain as
 * a set of connected blocks.
 * It is assumed that <code>PBlock</code> are not connected cyclically
 * residues have a head and a tail and they  connect together head to tail.
*/

class PResidue {
 public:
  friend class PChain;
  friend struct DOFCacher;

  /**
   * Constructs a new <code>PResidue</code> as a member of the
   * specified PChain <code>loop</code>.  This constructor uses
   * default positions for all blocks in the residue shell specified
   * by <code>shellName</code>, starting with default positions
   * for the head block.  After construction, the new <code>PResidue</code>
   * is "isolated", i.e. it is not connected to any other residues
   */

  PResidue(PChain *loop, const string &shellName);

  /**
   * Constructs a new <code>PResidue</code> as a member of the
   * specified PChain <code>loop</code>.  This constructor specifies
   * the residue shell (<code>shellName</code>) for the new residue
   * but also uses the atom positions defined in <code>spec</code>.
   * After construction, the new <code>PResidue</code> is "isolated", i.e.
   * it is not connected to any other residues.
   */

  PResidue(PChain *loop, const string &shellName, PResidueSpec &spec); 

  /**
   * Constructs a new <code>PResidue</code> as a member of the
   * specified PChain <code>loop</code>.  This constructor uses
   * default positions for all blocks in the residue shell specified
   * by <code>shellName</code>, starting with default positions
   * for the head block.  After construction, the new <code>PResidue</code>'s
   * "head" is connected to the tail of <code>toConnect</code>.
   */

  PResidue(PChain *loop, const string &shellName, PResidue *toConnect);

  /**
   * Constructs a new <code>PResidue</code> as a member of the
   * specified PChain <code>loop</code>.  This constructor specifies
   * the residue shell (<code>shellName</code>) for the new residue
   * but also uses the atom positions defined in <code>spec</code>.
   * After construction, the new <code>PResidue</code>'s "head"
   * is connected to the tail of <code>toConnect</code>.
   */
  
  PResidue(PChain *loop, const string &shellName, PResidueSpec &spec, PResidue *toConnect);

  ~PResidue();


  /**
   * Yajia Zhang added
   * */
  bool applyRotamer( const int&, const vector<double>&);

  void getChiAngles( vector<double>& chis);

  /**
   * Yajia Zhang added
   * */
  void calChi( vector<double>& chis);

  /**
   * Returns the next residue in the chain following this
   * <code>PResidue</code>, or <code>NULL</code> if this
   * residue is not part of a chain.
   */

  PResidue *NextResidue() const { return m_nextRes; }

  /**
   * Returns the previous residue in the chain before this
   * <code>PResidue</code>, or <code>NULL</code> if this
   * residue is not part of a chain.
   */

  PResidue *PreviousResidue() const { return m_prevRes; }

  /**
   * Returns the DOF between the two given atom ID's in the specified block type.
   */

  PBond *getDOF(const string &blockType, const string &atomId1, const string &atomId2);


  /**
   * Returns a pointer to a mapping of atom Id's to bonds for all DOF's associated with the given block. 
   * Use an iterator to iterate through all the bonds.
   *
   * Returns NULL if it can't find the blockType.
   * Return should NOT be deleted.
   */

  const ID_To_DOF_Map *getDOFsForBlock(const string &blockType);

  /**
   * Randomizes the DOFs associated with the given blockType.
   */

  void RandomizeDOFs(const string &blockType, PChain *rootChain);
  void RandomizeDOFs(const string &blockType);

  void DetachBlocks(const string &blockType, const string &blockToDetachFrom);
  void ReattachBlocks(const string &blockType);
  void ReattachAllBlocks();

  /**
   * Returns bond between given atoms, or NULL if bond or either atom does not exist.
   */
  PBond *getBond(const string &id1, const string &id2);

  /**
   * Returns a pointer to this residue's parent
   * chain.
   */

  PChain* getChain() const;

  /**
   * Returns this residue's name, currently just
   * the name of its underlying <code>PResidueShell</code>.
   */

  string getResourceName() const { return m_shell->getName(); }

  /** 
   * Returns this residue's name which can be custom set
   */
   string getName() const { return m_customName; }

  /** 
   * Sets the name of the residue
   */
   void setName(string name) { m_customName = name; }

  /**
   * Returns a pointer to this residue's head block.
   */

  PBlock *getHeadBlock() { return m_blocks[m_shell->getHeadId()]; }

  /**
   * Returns a pointer to this residue's tail block.
   */

  PBlock *getTailBlock() { return m_blocks[m_shell->getTailId()]; }
  
  /**
   * Returns a pointer to a <code>vector</code> of all
   * atoms in this residue.
   */

  vector<PAtom *> *getAtoms() { return &m_atomCache; }

  /**
   * Returns a pointer to a <code>hash_map</code> mapping
   * all atom IDs to atom pointers in this residue.
   */

  HASH_MAP_STR(PAtom *) *getAtomMap() { return &m_atomMap; }

  /**
   * Returns a pointer to the <code>PAtom</code> in this residue
   * with the specified ID, or <code>NULL</code> if no atom with
   * this ID exists in the residue.
   */

  PAtom *getAtom(const string &id) { 
	if(m_atomMap.find(id)==m_atomMap.end()) return NULL;
	else
		return m_atomMap[id]; 
  }

  /**
   * Returns the <code>Vector3</code> position of the atom
   * in this residue with the specified ID.
   */

  Vector3 getAtomPosition(const string &id) { return getAtom(id)->getPos(); }

  /**
   * Returns true if any atom in this <code>PResidue</code>
   * is in static collision, false otherwise.
   */

  bool InStaticCollision() const;

  /**
   * Returns true if any atom in this <code>PResidue</code>
   * is in self collision, false otherwise.
   */

  bool InSelfCollision() const;

  /**
   * Returns true if any atom in this <code>PResidue</code>
   * is in either static or self collision, false otherwise.
   */

  bool InAnyCollision() const;
  
  /**
   * Returns a pair of atoms in static collision within
   * this <code>PResidue</code>, or <code>(NULL, NULL)</code>
   * if no such atoms exist.
   */

  pair<PAtom *, PAtom *> FindStaticCollision() const;

  /**
   * Returns a pair of atoms in self collision within
   * this <code>PResidue</code>, or <code>(NULL, NULL)</code>
   * if no such atoms exist.
   */

  pair<PAtom *, PAtom *> FindSelfCollision() const;

  /**
   * Returns a pair of atoms in either static or self collision within
   * this <code>PResidue</code>, or <code>(NULL, NULL)</code>
   * if no such atoms exist.
   */

  pair<PAtom *, PAtom *> FindAnyCollision() const;

  /**
   * Returns a set of all atoms in static collision
   * in this <code>PResidue</code>.
   */

  AtomCollisions* getAllCollidingStatic() const;

  /**
   * Returns a set of all atoms in self collision
   * in this <code>PResidue</code>.
   */

  AtomCollisions* getAllCollidingSelf() const;

  /**
   * Returns a set of all atoms in either static or
   * self collision in this <code>PResidue</code>.
   */

  AtomCollisions* getAllCollidingEither() const;

  /**
   * Returns the space manager of the chain this residue is in.
  */

  const PSpaceManager* getSpaceManager() const;

  /**
   * Gets a PResidueSpec with the information about this residue's 
   * atom positions.
   */
  
  PResidueSpec getSpec();

  /**
   * Sets the residue's atom positions from the PAtomPositionsSpec \c spec.
   * Atoms in the residue that are not found in \c spec are ignored.
   */
  void setPositions(const PAtomPositionsSpec &spec);

  /**
   * Returns the index of this residue in a PDB file.
   */
  int getPdbId() const { return m_pdb_id; }

  /**
   * Set the PDB index of this residue.
   */
  void setPdbId(int id) {m_pdb_id=id;}


  /**
   * Inactivate the residue so that it won't be considered in collision 
   */
  void inactivate();

  /**
   * Activate the residue so that it will be considered in collision
   */
  void activate();

  /**
   * Returns true if there is at least one atom in the residue is active for collision checking.
   */
  bool isActive();

  /**
   * Print out all the collision atom pairs in this residue.
  **/
  void printCollision();

  /*Yajia added: return the type index for the sidechain rotamer for probability evaluation.*/
  void getSideChainAngles( int& index, vector<double>& angles);

 private:

  void EstablishLink(PResidue *prior);
  //these two functions are called by PChain:
  void SetChain(PChain *loop) { assert(loop != NULL); m_loop = loop; }
  //caches atoms, block types
  void finalize();
  void ChangeBlockState(const string &blockType, bool status);

  //for use by PChain: //should be called after finalize is called
  void CacheDOF(DOF_Cache &DOFcache);
  void CacheAtoms(Atom_Cache &AtomCache);
  void DestroyBonds();

  void ConstructAroundHead(PChain *loop);
  void BuildWithPositions(PChain *loop, const PResidueSpec &spec);
  void AddAuxInfo(PResidueSpec &spec);

  // Returns true iff all of shell's atoms are defined in spec.
  bool IsFullyDefined(const string &shellName, const PResidueSpec &spec);


  // collision detection helper methods
  bool InCollision(CollisionType type) const;
  pair<PAtom *, PAtom *> FindCollision(CollisionType type) const;
  AtomCollisions* getAllColliding(CollisionType type) const;

  //for use by DOF_Cacher
  void CacheSingleDOF(const string &blockType,PBond *bond);

  /** Member data. **/
  PResidueShell *m_shell;
  HASH_MAP_STR(PBlock *) m_blocks;
  PChain *m_loop;
  string m_customName;

  // Index of this residue in a PDB file
  int m_pdb_id;
  
  // Atoms in this residue
  HASH_MAP_STR(PAtom *) m_atomMap;

  // Next and previous residues in the chain
  PResidue *m_nextRes, *m_prevRes;

  //caches:
  //turning blocks on and off should update this cache (or should it?)
  vector<PAtom *> m_atomCache;
  typedef HASH_MAP_STR(vector<PBlock *>) BlockTypeCache;
  BlockTypeCache m_blockTypeCache; 
  HASH_MAP_STR(HASH_MAP_STRPAIR_OR(PBond *)) m_dofs;
 public:
  //Yajia added, for indicating side chain type
  int type_sidechain;
  vector<double> angles_sidechain;
};

#endif  // __P_RESIDUE_H
