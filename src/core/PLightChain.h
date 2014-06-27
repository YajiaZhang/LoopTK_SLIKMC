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

#ifndef	__P_LIGHT_CHAIN_H
#define __P_LIGHT_CHAIN_H

#include "PAtom.h"
#include "PResidue.h"
//@package Main Infrastructure 
/**
 *
 * A lightweight chain that cannot be manipulated or have subchains. 
 * It can provide position info, collision detection, and other 
 * miscellaneous routines.
 */

class PLightChain {
 public:
  virtual int size() const=0;
  virtual PResidue *getResidue(int index)=0;
  virtual PAtom *getAtomAtRes(const string &atomId, int resNum)=0;
  virtual void setAtomColorAtRes(const string &atomId, int resNum, GLColor color)=0;
  virtual void revertAtomColorAtRes(const string &atomId, int resNum)=0;
  virtual const PSpaceManager* getSpaceManager()=0;

  virtual bool InSelfCollision()=0;
  virtual bool InSelfCollision(int resIndex1, int resIndex2)=0;
  virtual pair<PAtom *, PAtom *> FindSelfCollision()=0;
  virtual pair<PAtom *, PAtom *> FindSelfCollision(int resIndex1, int resIndex2)=0;
  virtual AtomCollisions* getAllCollidingSelf()=0;
  virtual AtomCollisions* getAllCollidingSelf(int resIndex1, int resIndex2)=0;
};

#endif  // __P_LIGHT_CHAIN_H
