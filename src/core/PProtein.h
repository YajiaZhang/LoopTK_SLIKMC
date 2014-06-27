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

#ifndef __P_PROTEIN_H
#define __P_PROTEIN_H

#include "PChain.h"
#include "PProteinResidue.h"


/**
 * A subclass of PChain which makes strong assumptions
 * about the properties the loop should have. Specifically,
 * convenience methods are provided to selectively manipulate
 * the backbone and sidechain or to retrieve anchor/end effector
 * atoms at the amino or carboxy terminus. This class should be used 
 * when using the builtin resources with %LoopTK.
 * 
 * Usage of a PProtein must follow the following sequence of steps:
 * 
 * -# \b Construction - all residues are added to chain. Use the
 *       PProtein::AddResidue() methods to do this.
 * -# \b Finalization - finishes construction of the protein and prepares it
 *       for manipulation. Use PProtein::finalize() to do this.
 * -# \b Analysis - the full range of functionality of PProtein can be used
 *       at this point. However, no more residues can be added.
 * 
 * Creating a PProtein using PDBIO automatically completes the first two steps.
 * 
 */

class PProtein: public PChain {
 public:
  
  /**
   * Constructs a new <code>PProtein</code> with no residues.
   */
  PProtein();

  /**
   * Constructs a new <code>PProtein</code> with one residue,
   * using the default positions for its atoms.
   */
  PProtein(const string &firstResidueName); 

  /**
   * Constructs a new <code>PProtein</code> with one residue,
   * using the specified positions for its atoms.
   */
  PProtein(const string &firstResidueName, PResidueSpec &firstResidueSpec);

  /**
   * Constructs a new <code>PProtein</code> as a subchain of
   * <code>protein</code>, from residue <code>resStartIndex</code>
   * to residue <code>resEndIndex</code>.
   */
  PProtein(PProtein *protein, int resStartIndex, int resEndIndex);

  /**
   * Constructs a new <code>PProtein</code> with
   * <code>numResidues</code> random residues.
   */
  PProtein(int numResidues);
  
  /**
   * Returns the parent <code>PProtein</code> of this <code>PProtein</code>,
   * or <code>NULL</code> if this is the top level protein.
   */
  PProtein *getParent() { return (PProtein *) PChain::getParent(); }
  
  /**
   * Returns the top level protein in the hierarchy this protein belongs to. 
   */
  PProtein *getTopLevelProtein();

  /**
   * Returns the residue at the given index. Residues are indexed from 0 to size()-1 inclusive.
   */
  PProteinResidue *getResidue(int localIndex) { return (PProteinResidue *) PChain::getResidue(localIndex); }

  /**
   * Returns the local index of the residue with PDB index <code>pdb_index</code>.
   */
  int pdbIndexToLocalIndex(int pdb_index);
  
  /**
   * Returns the residue with the specified PDB index, or NULL if no such
   * residue exists in the <code>PProtein</code>.
   */
  PProteinResidue *getResidueByPdbIndex(int pdbIndex);

  /**
   * Returns a subchain of this <code>PProtein</code> in which the first
   * residue has PDB index <code>startPdbIndex</code>, the second residue
   * has PDB index <code>startPdbIndex + 1</code>, ..., and the last
   * residue has PDB index <code>endPdbIndex</code>.
   */
  PProtein *getSubchainByPdbIndices(int startPdbIndex, int endPdbIndex);

  /**
   * Returns the number of degrees of freedom in this
   * <code>PProtein</code>'s backbone.
   */
  int NumBackboneDOFs() const;

  /**
   * Returns the number of degrees of freedom in this
   * <code>PProtein</code>'s side chains.
   */
  int NumSidechainDOFs() const;

  /**
   * Randomizes the DOFs of the side-chain at the given residue index.
   */
  void RandomizeSidechainAtRes(int resIndex);

  /**
   * Randomizes all side-chains in this protein.
   */
  void RandomizeAllSidechains();

  /**
   * Returns the loop closure end effectors on the 
   * specified <code>side</code> of this <code>PProtein</code>
   * (front or back).
   */
  void GetEndEffectors(ProteinSide side, PAtom *&endAtom, PAtom *&endPriorAtom);

  /**
   * Returns the loop closure anchors on the
   * specified side of the protein (front or back).
   */
  void GetAnchors(ProteinSide side, PAtom *&endAnchor, PAtom *&endPriorAnchor);
  
  /**
   * Disables all side chain blocks in the protein.
   */

  void DisableSidechains();

  /**
   * Disables all side chain blocks between the
   * specified start and end residues.
   */
  void DisableSidechains(int resIndex1, int resIndex2);
  
  /**
   * Enables all side chain blocks in the protein.
   */
  void EnableSidechains();

  /**
   * Enables all side chain blocks between the
   * specified start and end residues.
   */
  void EnableSidechains(int resIndex1, int resIndex2);

  /**
   * Rotates the backbone <code>numDegrees</code> degrees,
   * starting from the degree of freedom <code>DOF_index</code>
   * in the direction <code>dir</code>.
   */
  void RotateBackbone(int DOF_index, BondDirection dir, Real numDegrees);

  /**
   * Rotates the sidechain at index <code>DOF_index</code> by
   * <code>numDegrees</code> degrees.
   */
  void RotateSidechain(int DOF_index, float numDegrees);

  /**
   * Creates an exact copy of the protein. If a sub-chain is being copied,
   * it creates a top-level protein as well. So, to free the memory of the
   * sub-chain, Obliterate() must be called.
   */
  PProtein *Clone();

  /**
   * Randomize backbone DOFs, but does not check collision
   */
  void RandomizeBackbone();

  /**
   * Returns the child PProtein to this PProtein at the given index. Use
   * NumChildren() to get the number of child PProtein's.
  */
  PProtein *getChildProtein(unsigned int index) {
    return (PProtein *) getChild(index);
  }

  /**
   * Print all the atom pairs that are in collision.
  **/
  void printCollision();

protected:
  virtual PResidue *CreateResidue(const string &name) {
    return new PProteinResidue(this, name);
  }

  virtual PResidue *CreateResidue(const string &name, PResidueSpec &spec) {
    return new PProteinResidue(this, name, spec); 
  }

  virtual PResidue *CreateResidue(const string &name, PResidue *res) {
    return new PProteinResidue(this, name, res);
  }

  virtual PResidue *CreateResidue(const string &name, PResidueSpec &spec, PResidue *res) {
    return new PProteinResidue(this, name, spec, res);
  }

};



#endif  // __P_PROTEIN_H
