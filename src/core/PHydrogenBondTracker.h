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

#ifndef PHYDRO_BOND_TRACKER
#define PHYDRO_BOND_TRACKER

#include "PBasic.h"
#include "PConstants.h"
#include <set>




//should i make this a global float that can be changed by the user through some static function?
#define HYDRO_BOND_CUTOFF 4.5


typedef pair<PAtom *, PAtom *> PHydrogenBond;

struct HydroBondCompare {
  //first atom will be N, second will be O
  PHydrogenBond ConvertHB(PHydrogenBond hb) const {
    PHydrogenBond ret;
    if (hb.first->getName()==PID::N) {
      ret.first = hb.first;
      ret.second = hb.second;
  } else {
      ret.first = hb.second;
      ret.second = hb.first;
    }
    return ret;
  }
  
  bool operator() (const PHydrogenBond &hb1, const PHydrogenBond &hb2) const {
    PHydrogenBond comp1 = ConvertHB(hb1);
    PHydrogenBond comp2 = ConvertHB(hb2);
    if (comp1.first==comp2.first) {
      if (comp1.second==comp2.second) return false;
      else
      return comp1.second<comp2.second;
    } else {
      return comp1.first<comp2.first;
    }  
  }
};

typedef set<PHydrogenBond,HydroBondCompare> HydroBondSet;

class PHydrogenBondTracker {
 public:
  virtual HydroBondSet& GetHydrogenBonds()=0;
  static PHydrogenBondTracker *Create(PProtein *toTrack);
};









  ////IMPLEMENTATION FOLLOWS FROM THIS POINT


class PHydrogenBondTrackerObject:public PRotateEventHandler, public PHydrogenBondTracker {
 public:
  friend class PHydrogenBondTracker;
  ~PHydrogenBondTrackerObject();
  HydroBondSet& GetHydrogenBonds();
  void HandleRotation(PChain *p,ChainMove justExecuted);
 private:
  PHydrogenBondTrackerObject(PProtein *toTrack);
  void GenerateAtomsToCheck();
  void FindHydrogenBonds();
  PProtein *m_trackedProtein;
  HydroBondSet m_hydroBonds;
  list<PAtom *> m_atomsToCheck;
  bool m_mustRegenerateBonds;
};





#endif
