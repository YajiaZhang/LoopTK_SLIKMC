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

#ifndef PRESOURCES_H
#define PRESOURCES_H

#include "PBasic.h"
#include "PHashing.h"

//first block is block that's already defined, second is block to define
typedef HASH_MAP_STRPAIR_OR(PBlockConnection *) BlockConnectionData;

typedef HASH_MAP_STR(PAtomShell *) AtomShellData;
typedef HASH_MAP_STR(PBlockShell *) BlockShellData;
typedef HASH_MAP_STR(PResidueShell *) ResidueShellData;
//@package Resource Management
/**
 *
 * A resource manager allowing lookups and
 * insertions of atom shells, block shells,
 * block connections and residue shells by
 * name.
 */

class PResources {

 friend class PInit;

 public:

  /**
   * Empties the resource manager and frees
   * all memory associated with existing
   * resources.
   */

  static void FreeResources();

  /**
   * Returns the atom shell mapped to the specified name. Throws a
   * resource error if the atom shell has not been inserted.
   */

  static PAtomShell *GetAtomShell(const string &name);

  /**
   * Returns the block shell mapped to the specified name. Throws a
   * resource error if the block shell has not been inserted.
   */

  static PBlockShell *GetBlockShell(const string &name);

  /**
   * Returns the block connection mapped to the specified pair of
   * defined-toDefine blocks. Throws a resource error if this block
   * connection has not been inserted.
   */

  static PBlockConnection *GetBlockConnection(const string &definedName, const string &toDefineName); 

  /**
   * Returns the residue shell mapped to the specified name. Throws a
   * resource error if the block shell has not been inserted.
   */

  static PResidueShell *GetResidueShell(const string &name);

  /**
   * Returns the epsilon (van der Waals energy coefficient) value for
   * the specified element pair, e.g. <code>("C", "S")</code>.
   */

  static Real GetEpsilonValue(const StringPair &elemPair);

  /**
   * Returns the new atom ID to use, if defined, when the specified
   * <code>id</code> is read from PDB in a residue named <code>resName</code>.
   */

  static string GetAtomIDMapping(const string &id, const string &resName);

  /**
   * Returns the pair of atom IDs that form the specified chi angle
   * in the specified residue.
   */

  static vector<string> GetChiIndex(const string &resName, int chiIndex);

  /**
   * Returns all rotamer data for the specified residue.
   */

  static vector<vector<Real> > GetRotamer(const string &resName);

  /**
   * Returns the number of rotamer configurations available for the specified residue.
   */

  static int GetRotamerSize(const string &resName);

  /**
   * Returns <code>true</code> if the resource manager contains
   * the atom shell <code>name</code>, <code>false</code> otherwise.
   */

  static bool ContainsAtomShell(const string &name) {
    return(m_atoms.find(name) != m_atoms.end());
  }

  /**
   * Returns <code>true</code> if the resource manager contains
   * the block shell <code>name</code>, <code>false</code>otherwise.
   */

  static bool ContainsBlockShell(const string &name) {
    return(m_blocks.find(name) != m_blocks.end());
  }

  /**
   * Returns <code>true</code> if the resource manager contains
   * a block connection between <code>definedName</code> and
   * <code>toDefineName</code>, <code>false</code> otherwise.
   */

  static bool ContainsBlockConnection(const string &definedName, const string &toDefineName) {
    return(m_blockConnections.find(make_pair(definedName, toDefineName)) != m_blockConnections.end());
  }

  /**
   * Returns <code>true</code> if the resource manager contains
   * the residue shell <code>name</code>, <code>false</code> otherwise.
   */

  static bool ContainsResidueShell(const string &name) {
    return(m_residues.find(name) != m_residues.end());
  }

  /**
   * Returns <code>true</code> if the resource manager contains
   * an atom ID mapping to use when atom <code>id</code> is read from PDB
   * in a residue named <code>resName</code>, <code>false</code> otherwise.
   */

  static bool ContainsAtomIDMapping(const string &id, const string &resName) {
    return(m_atomIDs.find(make_pair(resName, id)) != m_atomIDs.end());
  }

  /**
   * Returns <code>true</code> if the resource manager contains
   * a mapping for the specified <code>chiIndex</code> in the
   * specified residue.
   */

  static bool ContainsChiIndex(const string &resName, int chiIndex) {
    return(m_chiIndices.find(make_pair(resName, PUtilities::toStr(chiIndex))) != m_chiIndices.end());
  }

  /**
   * Returns <code>true</code> if the resource manager contains
   * rotamer data for residue <code>resName</code>, <code>false</code>
   * otherwise.
   */

  static bool ContainsRotamer(const string &resName) {
    return(m_Rotamers.find(resName) != m_Rotamers.end());
  }

  /**
   * Returns the number of atom shells
   * currently in the resource manager.
   */

  static int numAtoms() { return m_atoms.size(); }

  /**
   * Returns the number of block shells
   * currently in the resource manager.
   */

  static int numBlocks() { return m_blocks.size(); }

  /**
   * Returns the number of block connections
   * currently in the resource manager.
   */

  static int numConnections() { return m_blockConnections.size(); }

  /**
   * Returns the number of residue shells
   * currently in the resource manager.
   */

  static int numResidues() { return m_residues.size(); }

  /**
   * Returns the number of chi indices
   * in the resource manager for residue
   * <code>resName</code>.
   */

  static int numChiIndices(const string &resName);

  /**
   * Returns the names of all atom shells
   * currently in the resource manager.
   */

  static vector<string> getAtomNames() { return m_atomNames; }

  /**
   * Returns the names of all block shells
   * currently in the resource manager.
   */

  static vector<string> getBlockNames() { return m_blockNames; }

  /**
   * Returns the names of all block connections
   * currently in the resource manager.
   */

  static vector<StringPair> getConnectionNames() { return m_connectionNames; }

  /**
   * Returns the names of all residue shells
   * currently in the resource manager.
   */

  static vector<string> getResidueNames() { return m_residueNames; }
 
 private:
  static void ResourceError(const string &id);

  static AtomShellData m_atoms;
  static BlockConnectionData m_blockConnections;
  static BlockShellData m_blocks;
  static ResidueShellData m_residues;

  static HASH_MAP_STRPAIR_EX(string) m_atomIDs;
  static HASH_MAP_STRPAIR_EX(vector<string>) m_chiIndices;
  static HASH_MAP_STRPAIR_OR(Real) m_epsilonVals;
  static HASH_MAP_STR(vector<vector<Real> >) m_Rotamers;

  static vector<string> m_atomNames, m_blockNames, m_residueNames;
  static vector<StringPair> m_connectionNames;
  
  static void AddAtomShell(PAtomShell *atomShell);
  static void AddBlockShell(PBlockShell *blockShell);
  static void AddBlockConnection(PBlockConnection *blockConnection);
  static void AddResidueShell(PResidueShell *resShell);

  static void AddAtomIDMapping(const string &resName, const string &fromID, const string &toID);
  static void AddChiIndex(const string &resName, int chiIndex, const vector<string> &bondedAtoms);
  static void AddEpsilonValue(const StringPair &atomTypes, Real val);
  static void AddRotamer(const string &resName, const vector<string> &chiDegrees);

  template <typename T>
  static void FreeResMap(T& resMap) {
    for(typename T::iterator it = resMap.begin(); it != resMap.end(); ++it) {
      delete it->second;
    }

    resMap.clear();
  }

};

#endif
