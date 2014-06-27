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

#ifndef __PPROTEIN_RESIDUE_H
#define __PPROTEIN_RESIDUE_H

#include "PResidue.h"
#include "PBasic.h"

using namespace std;
class PProtein;


/**
 * This class represents a residue within PProtein, so it makes assumptions on 
 * the contents within. More specifically, this residue has a backbone and sidechain 
 * and depends on the configuration data provided with LoopTK. 
 * Do not create PProteinResidue
 * directly - use PProtein::AddResidue instead. 
 * See PResidue, the superclass, for more operations.
 */
class PProteinResidue: public PResidue {
 public:
  PProteinResidue(PProtein *loop, const string &shellName): PResidue((PChain *) loop,shellName) {currRot=0;}

  PProteinResidue(PProtein *loop, const string &shellName, PResidueSpec &spec): PResidue((PChain *) loop,shellName,spec) 
{currRot=0;}

  PProteinResidue(PProtein *loop, const string &shellName, PResidue *toConnect): PResidue((PChain *) loop,shellName,toConnect) 
{currRot=0;}

  PProteinResidue(PProtein *loop, const string &shellName, PResidueSpec &spec, PResidue *toConnect): PResidue((PChain *) 
loop,shellName,spec,toConnect) 
{currRot=0;}

  /**
   * Returns the next residue in the protein chain or NULL if this is the last residue.
   */ 
  PProteinResidue *NextResidue() { return (PProteinResidue *) PResidue::NextResidue(); }

  /**
   * Returns the previous residue in the protein chain or NULL if this is the first residue.
   */ 
  PProteinResidue *PreviousResidue() { return (PProteinResidue *) PResidue::PreviousResidue(); }

  /**
   * Returns the protein which directly owns this residue (the lowest subchain in the protein heirarchy.
   */
  PProtein *getProtein() { return (PProtein *) PResidue::getChain(); } 

  /**
   * Resets the side chain position to saved position. Use SaveSideChain() to save position.
  */

  void ResetSideChain();

  /**
   * Cycles through collision free rotamers.
  */

  bool ApplyRotamer();

  /**
   * Apply a specific rotamer.
  */

  bool ApplyRotamer(int index);

  /**
   * Save the current side chain position. Use ResetSideChain() to apply saved position.
  */
  void SaveSideChain();

  //Yajia added
//  void applyRotamer(const int& index, const vector<double>& rotamer2);
   /**
   * Returns the current Phi angle of the residue.
   */
 
  Real GetPhi();


  /**
   * Returns the current Psi angle of the residue. 
   */

  Real GetPsi();

  /**
   * Sets the indexed Chi angle of the side chain in degrees.
   */

  void SetChi(Real chi, int index);

  /**
   * Returns the indexed Chi angle of the side chain.
   */

  Real GetChi(int index);
  bool TryApplyRotamer(const vector<Real> &rotamer);
 private:

  //given list of chi angles, rotate side chain and check for collision.
//  bool TryApplyRotamer(const vector<Real> &rotamer);

  //use by rotamer
  unsigned currRot;
  vector<Real> originalAngle;

};
















#endif
