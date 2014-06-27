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

#ifndef __P_SPACE_MANAGER_H
#define __P_SPACE_MANAGER_H

#include <list>
using std::list;
//@package Grid
 
/**
 *
 * An abstract class allowing for queries about the
 * space in which a given chain exists.  Each chain
 * has its own <code>PSpaceManager</code>.
*/

class PSpaceManager {
  public:

    /**
     * Returns a <code>list</code> of all atoms within the
     * given <code>distance</code> from the given <code>point</code>.
     */

    virtual list<PAtom *> AtomsNearPoint(const Vector3 &point, Real distance) const = 0;

    virtual pair<Vector3, Vector3> getCellBoundingBox(const Vector3 &pos) const = 0;
    virtual pair<Vector3, Vector3> getGlobalBoundingBox() const = 0;
    virtual OccupancyMap getOccupancy() const = 0;

};

#endif  // __P_SPACE_MANASGER_H
