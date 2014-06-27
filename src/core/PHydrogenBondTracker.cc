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

#include "PHydrogenBondTracker.h"


PHydrogenBondTracker *PHydrogenBondTracker::Create(PProtein *toTrack) {
  return new PHydrogenBondTrackerObject(toTrack);
}

PHydrogenBondTrackerObject::PHydrogenBondTrackerObject(PProtein *toTrack) {
  m_trackedProtein = toTrack;
  m_trackedProtein->AddRotateEventHandler(this);
  m_mustRegenerateBonds = true;
  GenerateAtomsToCheck();
}

PHydrogenBondTrackerObject::~PHydrogenBondTrackerObject() {
  m_trackedProtein->RemoveRotateEventHandler(this);

}


struct HydroAtomGetter: AtomFunctor {
public:
  HydroAtomGetter(list<PAtom *> *toUpdate) {
    m_toUpdate = toUpdate;
  }

  void operator()(PAtom *atom, PBond *bondFrom) {
    if (atom->getName()==PID::N||atom->getName()==PID::O)
      m_toUpdate->push_back(atom);
  }

private:
  list<PAtom *> *m_toUpdate;
};

void PHydrogenBondTrackerObject::GenerateAtomsToCheck() {
  HydroAtomGetter old_hag(&m_atomsToCheck);
  m_trackedProtein->traverseFromStart(&old_hag);
}

HydroBondSet& PHydrogenBondTrackerObject::GetHydrogenBonds() {
  if (m_mustRegenerateBonds) FindHydrogenBonds();
  return m_hydroBonds;
}

void PHydrogenBondTrackerObject::FindHydrogenBonds() {
  m_hydroBonds.clear();
  for(list<PAtom *>::iterator it = m_atomsToCheck.begin();it!=m_atomsToCheck.end();++it) {
    PAtom *checkIt = *it;
    string lookingFor = PID::N;
    if (checkIt->getName()==PID::N) lookingFor = PID::O;
    const PSpaceManager *space_manager = checkIt->getSpaceManager();
    list<PAtom *> atomsToCheck = space_manager->AtomsNearPoint(checkIt->getPos(),HYDRO_BOND_CUTOFF);
    for(list<PAtom *>::iterator it = atomsToCheck.begin();it!=atomsToCheck.end();++it) {
      PAtom *compareAtom = *it;
      if (compareAtom->getName()==lookingFor) {
  PHydrogenBond newHydroBond;
  newHydroBond.first = checkIt;
  newHydroBond.second = compareAtom;
  m_hydroBonds.insert(newHydroBond);
      }
    }
  }
  m_mustRegenerateBonds = false;
}

void PHydrogenBondTrackerObject::HandleRotation(PChain *p,ChainMove justExecuted) {
  m_mustRegenerateBonds = true;
}
