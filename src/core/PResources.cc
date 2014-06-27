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

#include "PResources.h"
#include "PUtilities.h"
#include "PHashing.h"

BlockConnectionData PResources::m_blockConnections;
AtomShellData PResources::m_atoms;
BlockShellData PResources::m_blocks;
ResidueShellData  PResources::m_residues;

HASH_MAP_STRPAIR_EX(string) PResources::m_atomIDs;
HASH_MAP_STRPAIR_EX(vector<string>) PResources::m_chiIndices;
HASH_MAP_STRPAIR_OR(Real) PResources::m_epsilonVals;
HASH_MAP_STR(vector<vector<Real> >) PResources::m_Rotamers;

vector<string> PResources::m_atomNames, PResources::m_blockNames, PResources::m_residueNames;
vector<StringPair> PResources::m_connectionNames;

void PResources::AddAtomShell(PAtomShell *atomShell) {
  m_atoms[atomShell->getName()] = atomShell;
  m_atomNames.push_back(atomShell->getName());
}

void PResources::AddBlockShell(PBlockShell *blockShell) {
  m_blocks[blockShell->getName()] = blockShell;
  m_blockNames.push_back(blockShell->getName());
}

void PResources::AddBlockConnection(PBlockConnection *blockConnection) {
  m_blockConnections[blockConnection->getName()] = blockConnection;
  m_connectionNames.push_back(blockConnection->getName());
}

void PResources::AddResidueShell(PResidueShell *resShell) {
  m_residues[resShell->getName()] = resShell;
  m_residueNames.push_back(resShell->getName());
}

void PResources::AddAtomIDMapping(const string &resName, const string &fromID, const string &toID)
{
  m_atomIDs[make_pair(resName, fromID)] = toID;
}

void PResources::AddChiIndex(const string &resName, int chiIndex, const vector<string> &bondedAtoms)
{
  m_chiIndices[make_pair(resName, PUtilities::toStr(chiIndex))] = bondedAtoms;
}

void PResources::AddEpsilonValue(const StringPair &atomTypes, Real val)
{
  m_epsilonVals[atomTypes] = val;
}

void PResources::AddRotamer(const string &resName, const vector<string> &chiDegrees)
{
  vector<Real> parsedChiDegrees;

  for(unsigned i = 0; i < chiDegrees.size(); i++) {
    parsedChiDegrees.push_back(atof(chiDegrees[i].c_str()));
  }
  m_Rotamers[resName].push_back(parsedChiDegrees);
}

PAtomShell* PResources::GetAtomShell(const string &name) {
  PAtomShell *ret = m_atoms[name];
  if (ret==NULL) ResourceError(name);
  return ret;
}

PBlockShell* PResources::GetBlockShell(const string &name) {
  PBlockShell *ret = m_blocks[name];
  if (ret==NULL) ResourceError(name);
  return ret;
}

PBlockConnection *PResources::GetBlockConnection(const string &definedName, const string &toDefine) {
  StringPair sp;
  sp.first = definedName;
  sp.second = toDefine;
  PBlockConnection *ret = m_blockConnections[sp];
  if (ret==NULL) ResourceError("BlockConnection(" + definedName + "," + toDefine + ")");
  return ret; 
}

PResidueShell *PResources::GetResidueShell(const string &name) {
  PResidueShell *ret = m_residues[name];
  if (ret==NULL) ResourceError(name);
  return ret;
}

Real PResources::GetEpsilonValue(const StringPair &elemPair) {
  HASH_MAP_STRPAIR_OR(Real)::const_iterator it = m_epsilonVals.find(elemPair);
  if (it == m_epsilonVals.end()) {
    ResourceError("Epsilon(" + elemPair.first + ", " + elemPair.second + ")");
  } else {
    return it->second;
  }
  return 0;
}

string PResources::GetAtomIDMapping(const string &id, const string &resName)
{
  HASH_MAP_STRPAIR_EX(string)::const_iterator it = m_atomIDs.find(make_pair(resName, id));
  if (it == m_atomIDs.end()) {
    ResourceError("atom ID map(" + resName + ", " + id + ")");
  } else {
    return it->second;
  }
}


vector<string> PResources::GetChiIndex(const string &resName, int chiIndex)
{
  HASH_MAP_STRPAIR_EX(vector<string>)::const_iterator it;

  it = m_chiIndices.find(make_pair(resName, PUtilities::toStr(chiIndex)));
  if (it == m_chiIndices.end()) {
    ResourceError("chi index(" + resName + ", " + PUtilities::toStr(chiIndex) + ")");
  } else {
    return it->second;
  }
}

vector<vector<Real> > PResources::GetRotamer(const string &resName)
{
  HASH_MAP_STR(vector<vector<Real> >)::const_iterator it = m_Rotamers.find(resName);

  if (it == m_Rotamers.end()) {
    ResourceError("rotamer(" + resName + ")");
  } else {
    return it->second;
  }
}

int PResources::GetRotamerSize(const string &resName){
  if (ContainsRotamer(resName)){
    vector<vector<Real> > rotamerAngles = GetRotamer(resName);
    return rotamerAngles.size();
  }
  else{
    return 0;
  }
}


int PResources::numChiIndices(const string &resName)
{
  int total = 0;

  for(int i = 1; ; i++) {
    if (ContainsChiIndex(resName, i)) {
      total++;
    } else break;
  }

  return total;
}

void PResources::FreeResources()
{
  FreeResMap(m_atoms);
  FreeResMap(m_blockConnections);
  FreeResMap(m_blocks);
  FreeResMap(m_residues);

  m_atomNames.clear();
  m_blockNames.clear();
  m_connectionNames.clear();
  m_residueNames.clear();
}

void PResources::ResourceError(const string &id) {
  PUtilities::AbortProgram("Resource \"" + id + "\" could not be found.");
}

