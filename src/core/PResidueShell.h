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

#ifndef __P_RESIDUE_SHELL_H
#define __P_RESIDUE_SHELL_H

//TODO: What is the purpose of having this class? Used in graph traversing?
/**
 * Encapsulates information common to all residues.
 * PResidueShell contains methods and data to store
 * and access the blocks and block connections within
 * a particular residue.
 */

class PResidueShell {
 public:

  /**
   * Constructs a new PResidueShell with the specified name,
   * head block, tail block, and ID of the atom where graph
   * traversals are to begin.  Use this constructor when the
   * head and tail blocks are not the same.
   */
  PResidueShell(const string &name,
                const string &headBlockId,
                const string &headBlockName, 
                const string &tailBlockId,
                const string &tailBlockName,
                const string &startAtomID);

  /**
   * Constructs a new PResidueShell with the specified name,
   * core block and ID of the atom where graph traversals
   * are to begin.  Use this constructor when the head and
   * tail blocks are the same.
   */
  PResidueShell(const string &name, const string &coreBlockId, const string &coreBlockName, const string &startAtomID);

  /**
   * Adds the specified block to the residue.
   */
  void addBlock(const string &id, const string &name);

  /**
   * Adds the specified block connection to
   * the residue, indicating the block already
   * defined to the block being defined.
   */
  void addBlockConnection(const string &idDefined, const string &idToDefine);

  /**
   * Returns the ID of the atom in this residue where
   * graph traversals (of atoms and their bonds) are
   * to begin.
   */
  string getStartAtomID() const { return m_startAtomID; }

  /**
   * Returns the name of this residue.
   */
  string getName() const { return m_name; }

  /**
   * Returns the ID of the head block in this residue.
   */
  string getHeadId() const { return m_headBlock; }

  /**
   * Returns the ID of the tail block in this residue.
   */
  string getTailId() const { return m_tailBlock; }

  /**
   * Returns the name of the head block in this residue.
   */
  string getHeadBlockShell() const { return m_blocks.find(m_headBlock)->second; }

  /**
   * Returns the name of the tail block in this residue.
   */
  string getTailBlockShell() { return m_blocks[m_tailBlock]; }

  /**
   * Returns a pointer to a vector of all the
   * inter-block bonds in this residue.
   */
  vector<StringPair>* getBlockBonds() { return &m_blockBonds; }

  /**
   * Returns a pointer to the map of block ID ->
   * block names for this residue.
   */
  StringMap* getBlocks() { return &m_blocks; }
  
 private:
  void ConstructorHelper(const string &name, const string &headBlockId, const string &headBlockName,
	const string &tailBlockId, const string &tailBlockName, const string &startAtomID);

  string m_name;
  string m_startAtomID;
  StringMap m_blocks;
  string m_headBlock;
  string m_tailBlock;
  vector<StringPair> m_blockBonds;
};

#endif  // __P_RESIDUE_SHELL_H
