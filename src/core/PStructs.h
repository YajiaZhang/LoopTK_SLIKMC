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

#ifndef __P_STRUCTS_H
#define __P_STRUCTS_H

/**
 *
 * Represents a single rotation of a portion of a chain.
 */

struct ChainMove {
  /** The type of block in which the DOF_index indexes into. */
	//NOTE: backbone or sidechain
  string blockType;

  /** Degree of freedom index from which the rotation begins. */
  int DOF_index;

  /** Direction (forward or back) of rotation. */
  BondDirection dir;

  /** Angle of rotation, in degrees. */
  Real degrees;
  inline void print()
  {
	  cout << "blockType:" << blockType << "\tDOF_index:" << DOF_index << "\tBondDirection:" << dir << "\tdegrees:" << this->degrees << endl;
  }
};

/**
 * 
 * Represents DOF specifications used in the optimization
 */
struct CDof {
  /** Direction (forward or back) of rotation. */
  BondDirection dir;
  /** The type of block in which the DOF_index indexes into. */
  string blockType;
  /** Degree of freedom index from which the rotation begins. */
  int DOF_index;
};

/**
 * 
 * Abstract class allowing users of <code>PChain</code>s to
 * specify custom behavior when the chain is rotated.
 */

class PRotateEventHandler {
 public:
  virtual void HandleRotation(PChain *p, ChainMove justExecuted)=0;
};

#endif  // __P_STRUCTS_H
