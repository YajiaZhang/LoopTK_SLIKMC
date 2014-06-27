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
#include "PResources.h"
#include "PExtension.h"

PBlockConnection::PBlockConnection(const string &definedBlock, const string &blockToDefine, 
  const HASH_MAP_STRPAIR_EX(Vector3) &relativePositions, const vector<PBlockBondShell> &blockConnections) {
  m_blockDefined = PResources::GetBlockShell(definedBlock);
  m_blockToDefine = PResources::GetBlockShell(blockToDefine);
  m_blockDefinedRelPositions = relativePositions;
  m_blockConnections = blockConnections;
}

void PBlockConnection::GenerateNewBlock(PBlock *blockDefined, PBlock *incompleteBlock) {
  vector<Vector3> definedRelativeCoordinates;
  vector<Vector3> toDefineLocalCoordinates;
  vector<Vector3> definedGlobalCoordinates;
  vector<string> relativetoDefineAtomID;
  vector<string> noRelativetoDefineAtomID;
  vector<Vector3> definedNewGlobalCoordinates;
  Vector3 newGlobalCoordinates;
  
  PAtomPositionsSpec newPositions;

  HASH_MAP_STR(string) *definedAtoms = m_blockDefined->getAtoms();
  HASH_MAP_STR(string) *toDefineAtoms = m_blockToDefine->getAtoms();  
  PAtomPositionsSpec * localCoordinates = m_blockToDefine->getDefaultPositions();
  
  //build global coord frame
  for(HASH_MAP_STR(string)::iterator itDefined = definedAtoms->begin(); itDefined!=definedAtoms->end(); ++itDefined) {
    definedGlobalCoordinates.push_back((blockDefined->m_atoms[itDefined->first])->getPos());
  }
  
  //store id of todefine -->  expects defined to be backbone .
  for(HASH_MAP_STR(string)::iterator itDefine = toDefineAtoms->begin(); itDefine!=toDefineAtoms->end(); ++itDefine) {
      
      HASH_MAP_STRPAIR_EX(Vector3)::iterator relativeCoordinatesit = m_blockDefinedRelPositions.find(make_pair((definedAtoms->begin())->first,itDefine->first)); 
      if (relativeCoordinatesit != m_blockDefinedRelPositions.end()){
        relativetoDefineAtomID.push_back(itDefine->first);
      } 
      else{
        noRelativetoDefineAtomID.push_back(itDefine->first);
      }       
  }
  
  //setup coord frame
  for(int i=0; i<relativetoDefineAtomID.size(); i++){
    for(HASH_MAP_STR(string)::iterator itDefined = definedAtoms->begin(); itDefined!=definedAtoms->end(); ++itDefined) {
      HASH_MAP_STRPAIR_EX(Vector3)::iterator relativeCoordinatesit = m_blockDefinedRelPositions.find(make_pair(itDefined->first,relativetoDefineAtomID[i]));
      if (relativeCoordinatesit != m_blockDefinedRelPositions.end()){
        definedRelativeCoordinates.push_back(relativeCoordinatesit->second);
      }
    }  
    toDefineLocalCoordinates.push_back(localCoordinates->getAtomPosition(relativetoDefineAtomID[i]));
    
    newGlobalCoordinates = PMath::LocalToGlobalCoord(definedRelativeCoordinates, definedGlobalCoordinates);
    newPositions.addAtom(relativetoDefineAtomID[i],newGlobalCoordinates);

    definedNewGlobalCoordinates.push_back(newGlobalCoordinates);
    definedRelativeCoordinates.clear();  
  }
  
  //define the rest w/ the new coord frame.
  for(int j=0; j<noRelativetoDefineAtomID.size(); j++){
    for(int i=0; i<3; i++){
      definedRelativeCoordinates.push_back(toDefineLocalCoordinates[i] - localCoordinates->getAtomPosition(noRelativetoDefineAtomID[j]));
    }
    newGlobalCoordinates = PMath::LocalToGlobalCoord(definedRelativeCoordinates, definedNewGlobalCoordinates);
    newPositions.addAtom(noRelativetoDefineAtomID[j],newGlobalCoordinates);
    definedRelativeCoordinates.clear();  
  }
  
  //setup block.  
  incompleteBlock->CreateBlock(incompleteBlock->m_residue,newPositions);
  BondBlocksTogether(blockDefined,incompleteBlock);
}

//make sure to update this function if blocks keep track of what they're connected to
void PBlockConnection::BondBlocksTogether(PBlock *blockDefined, PBlock *blockJustDefined) {
  for(int i=0;i<m_blockConnections.size();i++) {
    PBlockBondShell pbbs = m_blockConnections[i];
    new PBond(blockDefined->m_atoms[pbbs.getDefinedAtomId()],blockJustDefined->m_atoms[pbbs.getToDefineAtomId()],pbbs.isDOF());
  }

  blockDefined->m_connectedBlocks.push_back(blockJustDefined);
  blockJustDefined->m_connectedBlocks.push_back(blockDefined);
}
