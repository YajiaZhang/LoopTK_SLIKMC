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

#ifndef __P_BLOCK_SHELL_H
#define __P_BLOCK_SHELL_H

class PBondShell;

/**
 * Encapsulates data common to all blocks.
 */

class PBlockShell {
 public:

  /**
   * Constructs a new PBlockShell with the given
   * name and type (such as backbone or sidechain).
   */
  PBlockShell(const string &name, const string &type);

  /**
   * Adds an atom to this block: id is the atom's
   * unique identifier within the block (CA, CB,
   * CG, etc.) whereas name is the type of atom
   * (C, N, O, etc.) to add.
   */
  void addAtom(const string &id, const string &name, Vector3 relativePosition);

  /**
   * Adds a bond to this block. id1 and id2 specify
   * the identifiers of the bonded atoms. isDOF is true
   * if this bond is a degree of freedom, false otherwise.
   */
  void addBond(const string &id1, const string &id2, bool isDOF);

  /**
   * Adds the specified PBondShell to this block.
   */
  void addBond(const PBondShell &pbs);
  
  /**
   * Returns the block's name.
   */
  string getName() const { return m_name; }

  /**
   * Returns the block's type.
   */
  string getType() const { return m_type; }

  /**
   * Returns the list of bonds in this block.
   */
  vector<PBondShell>* getBonds() { return &m_bondSpecs; }

  /**
   * Returns the map of atom names to atom types
   * for this <code>PBlockShell</code>.
   */
  StringMap* getAtoms() {return &m_atoms; }

  /**
   * Returns the map of default positions for atoms
   * in this <code>PBlockShell</code>.
   */
  PAtomPositionsSpec *getDefaultPositions() {return &m_defaultRelPositions; }

 private:
  string m_name;        /* Block name, e.g. "TRP". */
  string m_type;        /* Block type, e.g. "sidechain". */
  StringMap m_atoms;
  vector<PBondShell> m_bondSpecs;
  PAtomPositionsSpec m_defaultRelPositions;
};

#endif  // __P_BLOCK_SHELL_H
