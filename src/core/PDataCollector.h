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

#ifndef PDATA_COLLECTOR_H
#define PDATA_COLLECTOR_H

#include "PBasic.h"
#include "PConstants.h"
#include "PExtension.h"
#include "PIKAlgorithms.h"

/**
 * A function that returns valid solved conformations starting from the 
 * randomized loop configuration that <code>loop</code> is in.
 */
typedef IKSolutions (*LoopConfFunction)(PProtein *loop);

/**
 * User interface to collect information about loop
 * conformations and output the results to a CS2
 * file.
 */
class PDataCollector {
 public:
  
  /**
   * Creates or overwrites the .cs2 file <code>cs2FileName</code> with the 
   * loop in <code>pdbFileName</code> from the residue indices given as
   * parameters. The .cs2 file will have only one conformation inside, that
   * from the pdb file.
   */
  PDataCollector(const string &cs2FileName, const string &pdbFileName,int loopStartResIndex, int loopEndResIndex);

  /**
   * Creates or overwrites the .cs2 file <code>cs2FileName</code> with 
   * <code>loop</code>. <code>loop</code> should be a sub-loop of a protein.
   */
  PDataCollector(const string &cs2FileName, PProtein *loop); 

  /**
   * Creates a data collector that will append to the given file.
   */
  PDataCollector(const string &cs2FileName);

  /**
   * Clones the data collector and makes it layer off of the new file given.
  */
  void Clone(const string &newActivecs2Filename);

  /**
   * Finds the requested amount of conformations with the collision status
   * necessary for a valid conformation given as a parameter. Behind the 
   * scenes, this function uses Cyclic Coordinate Descent (CCD) followed
   * by an analytical exact IK solver. If the resulting conformation does not
   * match the collision status necessary, it is discarded.
  */

  void FindConformations(int amt, bool shouldBeCollisionFree);
  
  /**
   * Finds the requested amount of conformations using a user-provided function
   * to find a conformation (which is given as a parameter a loop with 
   * randomized degrees of freedom). 
   */

  void FindConformations(int amt, LoopConfFunction FindOneConformation);


  /**
   * Visualizes a conformation and allows the user to manipulate it
   * by hand. This function never returns.
  */

  void VisualizeLoop();

  /**
   * Visualizes the conformation space. This function never returns.
  */

  void VisualizeConformationSpace();


 private:
  //helper functions
  void InitMembers(const string &cs2FileName);

  string m_cs2FileName;





};







#endif
