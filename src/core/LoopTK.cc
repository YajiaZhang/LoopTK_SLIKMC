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

#include "PExtension.h"

/*! \mainpage LoopTK: Protein Loop Kinematic Toolkit
 * 
 * \section Summary
 * 
 * %LoopTK is an object-oriented toolkit for modeling, manipulating, and analyzing 
 * proteins. This document describes the design of the toolkit and illustrates how the 
 * various operations that can be done on proteins are made efficient.
 * 
 * The core concept underlying the toolkit is the "chain." A "chain" is a linear 
 * sequence of "residues", where each "residue" can contain an arbitrary number of 
 * atoms connected in an arbitrary way. These definitions of "chain" and "residue" 
 * are extremely general and can handle a broader class of molecules than just 
 * proteins - i.e. molecules that are different from the standard 
 * backbone-sidechain structure of proteins.
 * 
 * An important additional layer of abstraction is added between "residues" and 
 * "atoms" - the "block". A residue is composed of an arbitrary amount of blocks, 
 * while a block defines an arbitrary amount of atoms. This is done for two 
 * reasons. First, there is a lot of overlap between residues in their structure. 
 * For example, nearly every residue has a backbone sequence of atoms.  Abstraction 
 * makes the data structure more efficient. Second, there should be a way to 
 * differentiate between different parts of the protein. Additional layer of 
 * abstraction allows each block to have "type" associated with it - i.e. backbone 
 * and sidechain, which is useful for manipulating specific subsection of the 
 * residue. 
 * 
 * The toolkit requires all resources to be defined before modeling any molecules. 
 * Provided with the toolkit is configuration data for modeling proteins. This data 
 * defines all the residues, blocks, and atoms necessary. %LoopTK uses "hard" sphere 
 * model for atoms and assumes that the nuclear arrangement is fixed.  Most 
 * residues consist of three blocks: the "backbone" block, the "sidechain" block, 
 * and a block consisting of the oxygen atom that connects to the backbone. %LoopTK 
 * represents the coordinate of the atoms using relative coordinates.  When the 
 * protein is initially constructed by connecting block to the next block, 
 * homogeneous coordinate matrix is used to transform the local coordinate of the 
 * block to the appropriate global position.   
 * 
 * %LoopTK can construct molecules that have predefined positions (currently PDB is 
 * supported), or it can be used to construct molecules using only the "typical" 
 * structure for each residue. For example, one could construct a protein by 
 * specifying a sequence of residue names, and the toolkit will automatically 
 * generate positions for all the atoms that correspond to these "typical" 
 * structures.  These "typical" structures supplied with the toolkit are generated 
 * from analyzing the PDB files available on the RSCB website and can be customized
 * by changing the configuration data.
 * 
 * %LoopTK supports variety of collision detection methods through PSpaceManager.  
 * Currently, the "grid" method of Halperin and Overmars (1998) is used, whereby
 * the space is subdivided into cubes of equal size, and each atom is placed in 
 * the cube which contains the atom's center. Each protein has its own grid
 * data structure, implemented by the PGrid class. The toolkit ensures that the grid is always 
 * updated whenever an atom moves so that it is accurate. Since all resources must 
 * be defined before usage, the maximum size of an atom is known before creating 
 * any grids. The size of the side length of each cube is set to be the maximum 
 * diameter of all defined atoms. Therefore, to check if an atom is in collision 
 * with another atom, the algorithm must check the cube the atom is in along with 
 * the 26 neighboring cubes.
 * 
 * %LoopTK provides a set of classes corresponding to the configuration data 
 * provided. These classes provide additional functionality for handling proteins. 
 * 
 * \section Modeling
 * 
 * The five core classes for modeling proteins are PAtom, PBond, PBlock, 
 * PResidue, and PChain. Only PChain is created directly by a user of the 
 * toolkit. A PChain must be fully defined before it can be manipulated and 
 * analyzed. Defining a PChain involves adding residues to the chain with or 
 * without forcing the positions of the atoms. The functions for PChain 
 * construction are: 
 * 
 * \code
void PChain::AddResidue(const string &resName) 
void PChain::AddResidue(const string &resName, const PResidueSpec &resSpec)\endcode
 *
 * 
 * A PResidueSpec defines position information for a residue. By not providing 
 * this, a PChain will automatically generate position information for the atoms 
 * based on the "typical" structure provided in the configuration data.
 * 
 * Once construction is complete, the PChain must be finalized by calling 
 * PChain::finalize(). Finalization signifies that construction is finished and the 
 * chain is ready to be manipulated and analyzed. During finalization, a PChain 
 * caches the data within the molecule in a number of different ways so that later 
 * operations done on the PChain are extremely fast. For example, consider the method:
 * \code
PChain::getAtom(string blockType, int i)\endcode
 * 
 * 
 * This routine selects the atom  at the ith index of the array of atoms belonging to
 * blocks with type  "blockType", where an atom's position in the array is determined
 * by the order in  which the atoms are traversed in a traversal of the chain from
 * start to finish.  Instead of doing this traversal each time the method is called,
 * the traversal is done once during finalization and the results are cached, allowing
 * the method to run in O(1) time. 
 * 
 * One of the strengths of PChain is that it is recursive. Any sequence of residues 
 * within a chain is also a chain. This allows %LoopTK to handle subdivision of the 
 * protein efficiently and is used in prioritized constraint satisfaction approach 
 * (Dhanik et al 2007) of sampling loop conformation, for example.
 * 
 * For example, suppose we have a PChain \c c with 50 residues. To 
 * create a subchain from the residues with indexes 10 to 20, use the following:
 * \code
PChain *subchain = new PChain(c, 10, 20)
\endcode A subchain is an additional 
 * owner of the same residues - operations which affect the residues affect all 
 * owning chains. Subchains have a local view on the residues - in this
 * example, residue 0 in \c subchain is the same as residue 10 in \c c. Because of 
 * the recursion, subchains can be created from subchains. Manipulating a subchain 
 * will only affect residues within the subchain. There is one constraint on 
 * creating subchains: direct children subchains of a chain cannot overlap. For 
 * example, the following two lines would cause an error:
 *
\code
PChain *c2 = new PChain(c, 10, 20);
PChain *c3 = new PChain(c, 15, 17);
\endcode
 *
 * 
 * However, the following commands would be valid:
 *
\code
PChain *c2 = new PChain(c, 10, 20);
PChain *c3 = new PChain(c2, 5, 7);
\endcode
 * 
 * As mentioned before, a PChain makes no assumptions on the data within it - it 
 * can represent any chain of molecules. To deal with proteins, the PProtein class 
 * is provided. PProtein is a subtype of PChain, and provides additional methods 
 * specific to the configuration data. For example, it has a method for getting the 
 * end effectors of a loop, which are the last two atoms of the backbone chain. 
 * Likewise, similar classes are provided for the other pieces of PChain - i.e. 
 * PProteinResidue for PResidue.
 * 
 * \section Configuration Data 
 * 
 * All configuration data is stored in the \c resources directory as XML files.  This
 * directory includes the following files:
 * - \c atoms.xml: Defines all atoms that might be needed by %LoopTK: their name, radius, color, and any other relevant data.
 * - \c blocks.xml: Defines all blocks used by %LoopTK.  A block is defined by specifying all atoms within the block and their relative positions to each other in a "typical" structure. To differentiate between atoms, each atom must be given an ID (i.e. "CA", "CB", or "N"). In addition, it defines which atoms are bonded together and which of those bonds are degrees of freedom.
 * - \c connections.xml: Defines block connections, which show how to connect two previously defined blocks together. Each block connection indicates which atoms between the blocks are bonded together, and whether those bonds are degrees of freedom. It also defines the relative positions between the blocks in a "typical" structure.
 * - \c residues.xml: Defines all residues used by %LoopTK.  A residue is defined by specifying which blocks are within the residue, and which blocks are connected together. Each atom in each block must have a unique ID within the whole residue. Additionally, the "core block" must be specified - the block which is used to connect this residue to other residues. Finally, a residue definition must specify the "start atom" - the atom in which traversals of this residue should start.
 * 
 * 
 * \section Input-Output
 * 
 * %LoopTK provides the module PDBIO for creating proteins from %PDB files and outputting proteins 
 * in memory to to %PDB files. For example, to create a protein with the data from \c 2CRO.pdb, use the
 * following code:
 * 
\code
PProtein *p = PDBIO::readFromFile("2CRO.pdb");
\endcode
 * 
 * Another module for input-output is CS2IO. A <tt>.cs2</tt> file defines a "chain space". 
 * A chain space is a collection of chains in which each chain is a variation of 
 * the others. In particular, each chain closes a particular subchain differently. 
 * All atoms outside of that subchain have the same positions for all chains in the 
 * chain space. <tt>.cs2</tt> files can be used to store all information about a chain 
 * space in a highly memory efficient way. 
 * 
 * \section Inverse Kinematics
 * 
 * A common task to accomplish with %LoopTK is closing a subchain (making it connect 
 * to the rest of the chain). %LoopTK provides two modules to do this, each using a 
 * different algorithm. One module uses Cyclic Coordinate Descent (CCD) (Canutescu, 
 * etl 2003), and the classes for doing this are PCCDSolver (for PChain) and 
 * PProteinCCDSolver (for PProtein). The other module uses exact inverse 
 * kinematics (Coutsias, etl 2004) and is called PExactIKSolver.
 */

bool LoopTK::m_warningsEnabled;

void LoopTK::Initialize(WarningIndicator ind, const string &resourcePath) {
  if (ind==SUPPRESS_WARNINGS) SuppressWarnings();
  else
    EnableWarnings();
  PInit::InitializeResources(resourcePath);
}

void LoopTK::SuppressWarnings() {
  m_warningsEnabled = false;
}


void LoopTK::EnableWarnings() {
  m_warningsEnabled = true;
}

void LoopTK::DisplayWarning(const string &warning)
{
  if (m_warningsEnabled) {
    cerr << "Warning: " << warning << endl;
  }
}
