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

#include "PBasic.h"
#include "PConstants.h"
#include "PExtension.h"
#include "PLibraries.h"
/*I changed here!*/
#include <memory>
using namespace std;
/*I changed here!*/
int PGrid::m_defaultSideLength;

/* Constructors and atom insertion/deletion functions. */

PGrid::PGrid() 
{
  m_sideLength = m_defaultSideLength;
  updateDelta();
}

void PGrid::addAtom(PAtom *atom)
{
  Vector3 gridPos = scaleToGrid(atom->getPos());
  CollisionMap::iterator found = m_collisionGrid.find(gridPos);

  /* Add to the collision grid. */
  if (found == m_collisionGrid.end()) {
    AtomSet newSet;
    newSet.insert(atom);
    m_collisionGrid[gridPos] = newSet;
  }
  else {
    found->second.insert(atom);
  }
}

void PGrid::removeAtom(PAtom *atom)
{
  if (atom == NULL) {
    PUtilities::AbortProgram("Can't remove null atom from the grid!");
  }

  Vector3 gridPos = scaleToGrid(atom->getPos());
  CollisionMap::iterator found = m_collisionGrid.find(gridPos);
  
  if (found == m_collisionGrid.end()) {
    PUtilities::AbortProgram("Didn't find requested atom (" + atom->getName() + " in the grid.");
  }

  /* Remove from the collision grid. */
  found->second.erase(atom);
  if (found->second.size() == 0) {
    m_collisionGrid.erase(found);
  }
}

void PGrid::changeAtomPos(PAtom* atom, Vector3& gridPos_prev, Vector3& gridPos_curr) {
	  if (atom == NULL) {
	    PUtilities::AbortProgram("Can't remove null atom from the grid!");
	  }
	  CollisionMap::iterator found = m_collisionGrid.find( gridPos_prev);
	  if (found == m_collisionGrid.end()) {
	    PUtilities::AbortProgram("Didn't find requested atom (" + atom->getName() + ") in the grid.");
//		PUtilities::AbortProgram("haha");

	  }
	  /* Remove from the collision grid. */
	  found->second.erase(atom);
	  if (found->second.size() == 0) {
	    m_collisionGrid.erase(found);
	  }

	  CollisionMap::iterator found2 = m_collisionGrid.find(gridPos_curr);

	  /* Add to the collision grid. */
	  if (found2 == m_collisionGrid.end()) {
	    AtomSet newSet;
	    newSet.insert(atom);
	    m_collisionGrid[gridPos_curr] = newSet;
	  }
	  else
	  {
	    found2->second.insert(atom);
	  }
	  return;
}

/* Edge length management functions. */

Vector3 PGrid::scaleToGrid(const Vector3 &pos) const
{
  int  x = int(floor(pos.x / Real(m_sideLength))),
    y = int(floor(pos.y / Real(m_sideLength))),
    z = int(floor(pos.z / Real(m_sideLength)));

  return Vector3(x, y, z);
}

void PGrid::updateSide(Real atomRadius)
{
  int atomDiameter = int(2*atomRadius) + 1;

  if (atomDiameter > m_defaultSideLength) {
    m_defaultSideLength = atomDiameter;
  }
}

void PGrid::setSideLength(int newVal)
{
  AtomSet allAtoms;

  /* Gather the set of all atoms currently in the grid. */
  for(CollisionMap::const_iterator i = m_collisionGrid.begin(); i != m_collisionGrid.end(); ++i) {
    for(AtomSet::const_iterator j = i->second.begin(); j != i->second.end(); ++j) {
      allAtoms.insert(*j);
    }
  }

  /* Rebuild the grid with the new edge length. */
  m_collisionGrid.clear();
  m_sideLength = newVal;
  updateDelta();

  for(AtomSet::const_iterator it = allAtoms.begin(); it != allAtoms.end(); ++it) {
    addAtom(*it);
  }
}

/* Collision detection functions. */

bool PGrid::atomsBonded(const PAtom *a1, const PAtom *a2)
{
  const vector<PBond *> &bonds = *(a1->getBonds());
  PAtom *bondFirst, *bondSecond;

  for(vector<PBond *>::const_iterator it = bonds.begin(); it != bonds.end(); ++it) {
    bondFirst = (*it)->getAtom1();
    bondSecond = (*it)->getAtom2();

    if ( (bondFirst == a1 && bondSecond == a2) || (bondFirst == a2 && bondSecond == a1) ) {
      return true;
    }
  }

  return false;
}

bool PGrid::InStaticCollision(const PAtom *atom) const
{
  return (getStaticCollidingAtom(atom) != NULL);
}

bool PGrid::InSelfCollision(const PAtom *atom) const
{
  return (getSelfCollidingAtom(atom) != NULL);
}

bool PGrid::InAnyCollision(const PAtom *atom) const
{
  return (InStaticCollision(atom) || InSelfCollision(atom));
}

PAtom* PGrid::getStaticCollidingAtom(const PAtom *atom) const
{
  auto_ptr<AtomSet> collSet(findColliding(atom, true, STATIC));

  if (collSet->size() == 0) {
    return NULL;
  } else {
	AtomSet::iterator it;
	for (it=collSet->begin(); it!=collSet->end(); ++it) {
		if ((*it)->isActive())
			return *it;
	}
	return NULL;

  }
}

PAtom* PGrid::getSelfCollidingAtom(const PAtom *atom) const
{
  auto_ptr<AtomSet> collSet(findColliding(atom, true, SELF));

  if (collSet->size() == 0) {
    return NULL;
  } else {
        AtomSet::iterator it;
        for (it=collSet->begin(); it!=collSet->end(); ++it) {
                if ((*it)->isActive())
                        return *it;
        }
        return NULL;

  }
}

PAtom* PGrid::getAnyCollidingAtom(const PAtom *atom) const
{
  PAtom *a = getStaticCollidingAtom(atom);
  if (a != NULL) return a;

  return getSelfCollidingAtom(atom);
}

AtomSet* PGrid::getAllCollidingStatic(const PAtom *atom) const
{
  return findColliding(atom, false, STATIC);
}

AtomSet* PGrid::getAllCollidingSelf(const PAtom *atom) const
{
  return findColliding(atom, false, SELF);
}

AtomSet* PGrid::getAllCollidingEither(const PAtom *atom) const
{
  return findColliding(atom, false, EITHER);
}

AtomSet* PGrid::findColliding(const PAtom *atom, bool getOnlyOne, CollisionType type, CollisionMetric inCollision) const
{
	AtomSet *s = new AtomSet();
	Vector3 centerPos = scaleToGrid(atom->getPos());

	for (int xDelta = -m_delta; xDelta <= m_delta; xDelta++)
	{
		for (int yDelta = -m_delta; yDelta <= m_delta; yDelta++)
		{
			for (int zDelta = -m_delta; zDelta <= m_delta; zDelta++)
			{
				Vector3 curPos(centerPos.x + xDelta, centerPos.y + yDelta, centerPos.z + zDelta);
				CollisionMap::const_iterator found = m_collisionGrid.find(curPos); // Loop atoms in this cell
				if (found != m_collisionGrid.end())
				{
					for (AtomSet::const_iterator it = found->second.begin();it != found->second.end(); ++it)
					{
						PAtom *curAtom = *it;
						if (curAtom->isActive())
						{
							if (inCollision(atom, curAtom) &&
									(type == EITHER
									|| (type == SELF && curAtom->getChain()->IsSubChainOf(atom->getChain()))
									|| (type == STATIC && !(curAtom->getChain()->IsSubChainOf(atom->getChain())))))
							{
								s->insert(curAtom);
								if (getOnlyOne)
									return s;
							}
						}
					}
				}
			}
		}
	}
	return s;
}

bool PGrid::covalentCollision(const PAtom *a1, const PAtom *a2)
{
  if (!a1->isActive() || !a2->isActive())
	  return false;

  if (a1 == a2)
	  return false;
  if (atomsBonded(a1, a2))
	  return false;

  Real  radiusSum = a1->getCovalentRadius() + a2->getCovalentRadius(),
    centerDist = a1->getPos().distance(a2->getPos());

  return (radiusSum >= centerDist);
}

bool PGrid::vanDerWaalsCollision(const PAtom *a1, const PAtom *a2)
{
  if ((!a1->isActive()) || (!a2->isActive()))
	  return false;

  if (a1 == a2)
	  return false;
  if (PAtom::shortestBondPath(a1, a2, BOND_THRESHOLD + 1) <= BOND_THRESHOLD)
	  return false;

  Real  radiusSum = a1->getVanDerWaalsRadius() + a2->getVanDerWaalsRadius(),
    centerDist = a1->getPos().distance(a2->getPos());

  return (COLLISION_THRESHOLD * radiusSum >= centerDist);
}

/* Occupancy functions. */

list<PAtom *> PGrid::AtomsNearPoint(const Vector3& point, Real distance) const
{
  list<PAtom *> returnList;
  Vector3 centerPos = scaleToGrid(point);
  int numToSearch = int(distance / Real(m_sideLength)) + 1;

  for(int xDelta = -numToSearch; xDelta <= numToSearch; xDelta++) {
    for(int yDelta = -numToSearch; yDelta <= numToSearch; yDelta++) {
      for(int zDelta = -numToSearch; zDelta <= numToSearch; zDelta++) {
        Vector3 curPos(centerPos.x + xDelta, centerPos.y + yDelta, centerPos.z + zDelta);

        CollisionMap::const_iterator found = m_collisionGrid.find(curPos);  // Loop atoms in this cell
        if (found != m_collisionGrid.end()) {
          for(AtomSet::const_iterator it = found->second.begin(); it != found->second.end(); ++it) {
            PAtom *curAtom = *it;
            if (curAtom->getPos().distance(point) <= distance) {
              returnList.push_back(curAtom);
            }
          }
        }
      }
    }
  }

  return returnList;
}

pair<Vector3, Vector3> PGrid::getCellBoundingBox(const Vector3 &pos) const
{
  Vector3 first = Vector3(pos.x * m_sideLength, pos.y * m_sideLength, pos.z * m_sideLength),
          second = Vector3((pos.x + 1) * m_sideLength, (pos.y + 1) * m_sideLength, (pos.z + 1) * m_sideLength);

  return make_pair<Vector3, Vector3>(first, second);
}

pair<Vector3, Vector3> PGrid::getGlobalBoundingBox() const
{
  Real minX = HUGE_VAL, minY = HUGE_VAL, minZ = HUGE_VAL,
       maxX = -HUGE_VAL, maxY = -HUGE_VAL, maxZ = -HUGE_VAL;
  vector<Vector3> result;

  /* Iterate through all cells in the grid with atoms defined. Update */
  /* the min/max coordinates of the 8 points of the bounding cube.    */
  for(CollisionMap::const_iterator it = m_collisionGrid.begin(); it != m_collisionGrid.end(); ++it) {
    if (it->first.x < minX) minX = it->first.x;
    if (it->first.y < minY) minY = it->first.y;
    if (it->first.z < minZ) minZ = it->first.z;

    if (it->first.x > maxX) maxX = it->first.x;
    if (it->first.y > maxY) maxY = it->first.y;
    if (it->first.z > maxZ) maxZ = it->first.z;
  }

  /* Return the points; order is important because it's assumed by getOccupancy(). */
  return make_pair<Vector3, Vector3>(Vector3(minX, minY, minZ), Vector3(maxX, maxY, maxZ));
}

OccupancyMap PGrid::getOccupancy() const
{
  pair<Vector3, Vector3> bounds = getGlobalBoundingBox(), curBox;
  OccupancyMap result;
  PAtom *curAtom;
  Vector3 curCell, newCell;

  /* Initialize the hash_map to show all cells as unoccupied. */
  for(int x = int(bounds.first.x); x <= bounds.second.x; x++) {
    for(int y = int(bounds.first.y); y <= bounds.second.y; y++) {
      for(int z = int(bounds.first.z); z <= bounds.second.z; z++) {
        result[Vector3(x, y, z)] = false;
      }
    }
  }

  /* For each atom, check occupancy of all its 27 neighboring cells (including its own cell). */
  for(CollisionMap::const_iterator i = m_collisionGrid.begin(); i != m_collisionGrid.end(); ++i) {
    AtomSet curSet = i->second;

    for(AtomSet::const_iterator j = curSet.begin(); j != curSet.end(); ++j) {
      curAtom = *j;
      curCell = scaleToGrid(curAtom->getPos());

      for(int x = -1; x <= 1; x++) {
        for(int y = -1; y <= 1; y++) {
          for(int z = -1; z <= 1; z++) {
            newCell = Vector3(curCell.x + x, curCell.y + y, curCell.z + z);
            curBox = getCellBoundingBox(newCell);

            if (PMath::sphereIntersectsCube(curAtom->getPos(), curAtom->getCovalentRadius(),
               curBox.first, curBox.second)) {
              result[newCell] = true;
            }
          }
        }
      }
    }
  }

  return result;
}

int PGrid:: getOccupancyNeighborNumber(const pair<Vector3, Vector3> &bound, OccupancyMap &occupancy, Vector3 pt){
    int num = -1;
    //checking bound
    if (pt.x<bound.first.x || pt.x > bound.second.x) return num;
    if (pt.y<bound.first.y || pt.y > bound.second.y) return num;
    if (pt.z<bound.first.z || pt.z > bound.second.z) return num;

    num = 0;
    /*
    //checking x-coordinate
    if (pt.x-1 > bound.first.x && occupancy[Vector3(pt.x-1, pt.y, pt.z)]==true) num++;
    if (pt.x+1 < bound.second.x && occupancy[Vector3(pt.x+1, pt.y, pt.z)]==true) num++;    

    //checking y-coordinate
    if (pt.y-1 > bound.first.y && occupancy[Vector3(pt.x, pt.y-1, pt.z)]==true) num++;
    if (pt.y+1 < bound.second.y && occupancy[Vector3(pt.x, pt.y+1, pt.z)]==true) num++;
    
    //checking z-coordinate
    if (pt.z-1 > bound.first.z && occupancy[Vector3(pt.x, pt.y, pt.z-1)]==true) num++;
    if (pt.z+1 < bound.second.z && occupancy[Vector3(pt.x, pt.y, pt.z+1)]==true) num++;
    */

    //checking 8 cells with z = z
    if ( occupancy[Vector3(pt.x-1, pt.y-1, pt.z)]==true) num++;
    if ( occupancy[Vector3(pt.x-1, pt.y, pt.z)]==true) num++;
    if ( occupancy[Vector3(pt.x-1, pt.y+1, pt.z)]==true) num++;
    if ( occupancy[Vector3(pt.x, pt.y-1, pt.z)]==true) num++;
    if ( occupancy[Vector3(pt.x, pt.y+1, pt.z)]==true) num++;
    if ( occupancy[Vector3(pt.x+1, pt.y-1, pt.z)]==true) num++;      
    if ( occupancy[Vector3(pt.x+1, pt.y, pt.z)]==true) num++;
    if ( occupancy[Vector3(pt.x+1, pt.y+1, pt.z)]==true) num++;
    
    //checking 9 cells with z = z-1
    if ( occupancy[Vector3(pt.x-1, pt.y-1, pt.z-1)]==true) num++;
    if ( occupancy[Vector3(pt.x-1, pt.y, pt.z-1)]==true) num++;
    if ( occupancy[Vector3(pt.x-1, pt.y+1, pt.z-1)]==true) num++;
    if ( occupancy[Vector3(pt.x, pt.y-1, pt.z-1)]==true) num++;
    if ( occupancy[Vector3(pt.x, pt.y, pt.z-1)]==true) num++;
    if ( occupancy[Vector3(pt.x, pt.y+1, pt.z-1)]==true) num++;
    if ( occupancy[Vector3(pt.x+1, pt.y-1, pt.z-1)]==true) num++;
    if ( occupancy[Vector3(pt.x+1, pt.y, pt.z-1)]==true) num++;
    if ( occupancy[Vector3(pt.x+1, pt.y+1, pt.z-1)]==true) num++;
    

    //checking 9 cells with z = z+1
    if ( occupancy[Vector3(pt.x-1, pt.y-1, pt.z+1)]==true) num++;
    if ( occupancy[Vector3(pt.x-1, pt.y, pt.z+1)]==true) num++;
    if ( occupancy[Vector3(pt.x-1, pt.y+1, pt.z+1)]==true) num++;
    if ( occupancy[Vector3(pt.x, pt.y-1, pt.z+1)]==true) num++;
    if ( occupancy[Vector3(pt.x, pt.y, pt.z+1)]==true) num++;
    if ( occupancy[Vector3(pt.x, pt.y+1, pt.z+1)]==true) num++;
    if ( occupancy[Vector3(pt.x+1, pt.y-1, pt.z+1)]==true) num++;
    if ( occupancy[Vector3(pt.x+1, pt.y, pt.z+1)]==true) num++;
    if ( occupancy[Vector3(pt.x+1, pt.y+1, pt.z+1)]==true) num++;
    
    return num;
}
