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

#ifndef __P_GRID_H
#define __P_GRID_H

#include "PEnums.h"
#include "PSpaceManager.h"
// @package Grid
/**
 *
 * An implementation of <code>PSpaceManager</code>, <code>PGrid</code>
 * represents a grid of atoms; it allows for efficient collision
 * detection, among other methods.
 */

typedef bool(*CollisionMetric)(const PAtom *, const PAtom *);

class PGrid: public PSpaceManager {
 public:
  friend class PAtom;
  friend class PAtomShell;
  friend class PChain;

  /**
  * Sets the atoms that are included within bond threshold for the collision test.
  */
  static const int BOND_THRESHOLD = 3;

  /**
  * Returns the length of the grid box.
  */
  Real getSideLength() const { return m_sideLength; }

  /**
  * Sets the length of the grid box.
  */
  void setSideLength(int newVal);

  /**
  * Given a point and a distance, returns list of atoms that 
  * satisfies the constraint.
  */
  list<PAtom *> AtomsNearPoint(const Vector3& point, Real distance) const;

  pair<Vector3, Vector3> getCellBoundingBox(const Vector3 &pos) const;
  pair<Vector3, Vector3> getGlobalBoundingBox() const;
  OccupancyMap getOccupancy() const;
  static int getOccupancyNeighborNumber(const pair<Vector3, Vector3> &bound, OccupancyMap &occupancy, Vector3 pt);

  //Yajia Zhang added!
  void changeAtomPos( PAtom* atom, Vector3& gridPos_prev, Vector3& gridPos_curr);
 private:
  
  PGrid();
  
  /* Private side length helper methods. */
  static void updateSide(Real atomRadius);
  void updateDelta()
  {
	  m_delta = int(ceil(double(m_defaultSideLength) / double(m_sideLength)));
  }
  
  /* Methods to add or remove atoms from the grid. */
  void addAtom(PAtom *atom);
  void removeAtom(PAtom *atom);

  /* Internal collision detection methods. */
  PAtom* getStaticCollidingAtom(const PAtom* atom) const;
  PAtom* getSelfCollidingAtom(const PAtom *atom) const;
  PAtom* getAnyCollidingAtom(const PAtom *atom) const;

  AtomSet* getAllCollidingStatic(const PAtom *atom) const;
  AtomSet* getAllCollidingSelf(const PAtom *atom) const;
  AtomSet* getAllCollidingEither(const PAtom *atom) const;

  AtomSet* findColliding(const PAtom *atom,
                         bool getOnlyOne,
                         CollisionType type,
                         CollisionMetric inCollision = vanDerWaalsCollision) const;

  bool InStaticCollision(const PAtom *atom) const;
  bool InSelfCollision(const PAtom *atom) const;
  bool InAnyCollision(const PAtom *atom) const;

  Vector3 scaleToGrid(const Vector3 &pos) const;

  /* Collision metrics. */
  static bool covalentCollision(const PAtom *a1, const PAtom *a2);
  static bool vanDerWaalsCollision(const PAtom *a1, const PAtom *a2);
  static bool atomsBonded(const PAtom *a1, const PAtom *a2);

  /* Private member variables. */
  //NOTE:The size of the grid
  static int m_defaultSideLength;   /* Default side length for new PGrids.    */
  int m_sideLength;                 /* Grid cell side length for THIS PGrid.  */

  int m_delta;                      /* How many cells in each x, y, z direction */
                                    /* to search for atom collisions; equal to  */
                                    /* ceil(m_defaultSideLength / m_sideLength. */
  //NOTE: A hash map: this finds all the atoms in the grid given the grid coordinate.
  CollisionMap m_collisionGrid;     /* The hash_map at the core of the PGrid. */

};

#endif  // __P_GRID_H
