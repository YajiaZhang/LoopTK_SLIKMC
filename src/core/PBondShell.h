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

#ifndef __P_BOND_SHELL_H
#define __P_BOND_SHELL_H

/**
 * Encapsulates data common to all
 * inter-atom bonds.
 */

class PBondShell {
 public:

  /**
   * Constructs a new PBondShell in an uninitialized
   * state. Do not call this constructor; it is only
   * included to prevent compiler errors.
   */

  PBondShell() {}

  /**
   * Constructs a new PBondShell between the specified
   * atoms. The PBondShell will represent a degree of
   * freedom if and only if isDOF is true.
   */

  PBondShell(const string &atomId1, const string &atomId2, bool isDOF) {
    m_atomId1 = atomId1;
    m_atomId2 = atomId2;
    m_isDOF = isDOF;
  }

  /**
   * Returns the ID of the first atom in this bond.
   */

  string getAtomId1() const { return m_atomId1; }

  /**
   * Returns the ID of the second atom in this bond.
   */

  string getAtomId2() const { return m_atomId2; }

  /**
   * Returns true if this bond is a degree of freedom,
   * false otherwise.
   */

  bool isDOF() const { return m_isDOF; }

 private:
  string m_atomId1;
  string m_atomId2;
  bool m_isDOF;
};

#endif  // __P_BOND_SHELL_H
