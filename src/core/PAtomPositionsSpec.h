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

#ifndef __P_ATOM_POSITIONS_SPEC
#define __P_ATOM_POSITIONS_SPEC

/**
 * Encapsulates a map of atom IDs
 * to their respective positions.
 */

class PAtomPositionsSpec {
 public:

  /**
   * Adds an atom with the specified ID
   * and specified position to the map.
   */
  void addAtom(const string &id, const Vector3 &position);

  /**
   * Returns the number of atom positions
   * currently in the map.
   */
  int size() const { return m_atomPositions.size(); }

  /**
   * Clears all data in the map.
   */
  void clear() { m_atomPositions.clear(); }

  /**
   * Returns true if the map contains <code>id</code>,
   * false otherwise.
   */
  bool contains(const string &id) const {
    return m_atomPositions.find(id) != m_atomPositions.end();
  }

  /**
   * Returns the position of the atom with
   * the specified ID. Aborts the program if
   * the atom does not exist in the map.
   */
  Vector3 getAtomPosition(const string &id) const;

 private:
  AtomPositions m_atomPositions;
};

#endif  // __P_ATOM_POSITIONS_SPEC
