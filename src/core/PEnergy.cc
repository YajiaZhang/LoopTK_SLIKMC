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
#include "PExtension.h"
#include "PResources.h"

class EnergyCalculator: public AtomFunctor {
  public:
  EnergyCalculator(Real threshold, const PSpaceManager *grid, EnergyCalcFn energyFn) {
    m_totalEnergy = 0;
    m_threshold = threshold;
    m_energyFn = energyFn;
    m_grid = grid;

    m_pairCache.clear();
  }

  void operator()(PAtom *atom, PBond *bondFrom) {
    list<PAtom *> neighbors = m_grid->AtomsNearPoint(atom->getPos(), m_threshold);
    PAtom *curAtom;

    for(list<PAtom *>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
      curAtom = *it;

      /* We don't include bonded atoms in the calculation, and we don't
       * want to double-count the same pair of atoms either. */
      int continueflag = 0;
      if (curAtom == atom) continue;
      if (curAtom->isBonded(atom)) continue;
      const vector<PBond *> CloseBonds = *(atom->getBonds());
      vector<PAtom *> CloseAtoms;
      for(int i=0;i<(CloseBonds).size();i++){
        if ((CloseBonds)[i]->getAtom1()!=atom)
        CloseAtoms.push_back((CloseBonds)[i]->getAtom1());
        else CloseAtoms.push_back((CloseBonds)[i]->getAtom2());
      }
      for(int i=0;i<CloseAtoms.size();i++){ if (curAtom->isBonded(CloseAtoms[i])){ continueflag = 1;break;}}
      if (continueflag == 1) continue;
      if (m_pairCache.find(make_pair(atom, curAtom)) != m_pairCache.end()) continue;

      m_totalEnergy += m_energyFn(atom, curAtom);
      m_pairCache.insert(make_pair(atom, curAtom));
    }
  }

  Real getTotalEnergy() {
    return m_totalEnergy;
  }

  private:
  Real m_totalEnergy;
  Real m_threshold;
  EnergyCalcFn m_energyFn;
  AtomCollisions m_pairCache;
  const PSpaceManager *m_grid;
};

Real PEnergy::vanDerWaalsEnergy(const PAtom *a1, const PAtom *a2)
{
  Real e = PResources::GetEpsilonValue(make_pair(a1->getName(), a2->getName()));
  Real s = (a1->getVanDerWaalsRadius() + a2->getVanDerWaalsRadius())/1.122;
  Real r = a1->getPos().distance(a2->getPos());

  return(4 * e * (Pow(s/r, 12) - Pow(s/r, 6)));
}

Real PEnergy::energyOfChain(PChain *chain, EnergyCalcFn energyFn, Real threshold)
{
  Real result = 0;
  EnergyCalculator *calc = new EnergyCalculator(threshold, chain->getSpaceManager(), energyFn);

  chain->traverseFromStart(calc);
  result = calc->getTotalEnergy();

  delete calc;
  return result;
}

Real PEnergy::collisionEnergy(const PAtom *a1, const PAtom *a2) {
  // d0: the two atoms' centers must be apart by COLLISION_THRESHOLD times the sum of their radius
  Real sum = a1->getVanDerWaalsRadius() + a2->getVanDerWaalsRadius();
  Real d0 = sum*COLLISION_THRESHOLD;
  Real d = a1->getPos().distance(a2->getPos());
  Real energy;
  if (d>d0)
    energy = 0;
  else
    energy = 1/Pow(d,2) - 1/Pow(d0,2);
  return energy;
}


