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
#include "PLibraries.h"
#include "PResources.h"
#include "PUtilities.h"

#include <libxml/parser.h>
#include <libxml/tree.h>

HASH_MAP_STRPAIR_EX(vector<PBlockBondShell>)    PInit::blockBonds;
HASH_MAP_STRPAIR_EX(HASH_MAP_STRPAIR_EX(Vector3))  PInit::relPos;


const string atomsFileName = "atoms.xml";
const string blocksFileName = "blocks.xml";
const string chiIndexFileName = "chi.xml";
const string connectionsFileName = "connections.xml";
const string epsilonFileName = "epsilon.xml";
const string mapsFileName = "maps.xml";
const string residuesFileName = "residues.xml";
const string rotamersFileName = "rotamers.xml";

void PInit::InitializeResources(const string &resourceDir)  {
  srand(unsigned(time(NULL)));
  string dir = resourceDir;
  if(resourceDir[resourceDir.size()-1]!='/') {
    dir+='/';
  }

  setupBlockShells(dir + blocksFileName);
  setupBlockConnections(dir + connectionsFileName);
  setupResidueShells(dir + residuesFileName);

  setupAtomShells(dir + atomsFileName);
  setupChiIndices(dir + chiIndexFileName);
  setupEpsilonVals(dir + epsilonFileName);
  setupIDMaps(dir + mapsFileName);
  setupRotamers(dir + rotamersFileName);
}

// XML helper functions.
string PInit::GetName(xmlNodePtr node) {
  return reinterpret_cast<const char *>(node->name);
}

string PInit::GetAttribute(xmlNodePtr node, const string &attr) {
  return reinterpret_cast<const char *>(xmlGetProp(node,
    reinterpret_cast<const xmlChar *>(attr.c_str())));
}

void PInit::setupAtomShells(const string &atomsFile) {
  xmlDocPtr doc = xmlReadFile(atomsFile.c_str(), NULL, 0);
  if (doc == NULL) {
    PUtilities::AbortProgram("Could not read " + atomsFile);
  }

  xmlNodePtr start_node = doc->children->children;
  for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
    if (GetName(cur_node) == "atom") {
      string atom_id = GetAttribute(cur_node, "id");
             
      Real covalent_radius = atof(GetAttribute(cur_node, "covalent_radius").c_str()),
           vanderwaals_radius = atof(GetAttribute(cur_node, "vanderwaals_radius").c_str());

      Real color_r = atof(GetAttribute(cur_node, "color_r").c_str()),
           color_g = atof(GetAttribute(cur_node, "color_g").c_str()),
           color_b = atof(GetAttribute(cur_node, "color_b").c_str()),
           color_a = atof(GetAttribute(cur_node, "color_a").c_str());

      GLColor::GLColor cur_color = GLColor::GLColor(color_r, color_g, color_b, color_a);

      PResources::AddAtomShell(new PAtomShell(covalent_radius, vanderwaals_radius, atom_id, cur_color));
    }
  }

  xmlFreeDoc(doc);
}

void PInit::setupBlockShells(const string &blocksFile) {
  xmlDocPtr doc = xmlReadFile(blocksFile.c_str(), NULL, 0);
  if (doc == NULL) {
    PUtilities::AbortProgram("Could not read " + blocksFile);
  }

  xmlNodePtr start_node = doc->children->children;
  for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
    if (GetName(cur_node) == "block") {
      // Initialize a new PBlockShell with this block name and type.
      string block_name = GetAttribute(cur_node, "name"),
             block_type = GetAttribute(cur_node, "type");
      PBlockShell *shell = new PBlockShell(block_name, block_type);

      for (xmlNodePtr cur_child_node = cur_node->children;
                      cur_child_node != NULL;
                      cur_child_node = cur_child_node->next) {
        if (GetName(cur_child_node) == "atom") {
          shell->addAtom(GetAttribute(cur_child_node, "name"),
                         GetAttribute(cur_child_node, "element"),
                         Vector3(atof(GetAttribute(cur_child_node, "x").c_str()),
                                 atof(GetAttribute(cur_child_node, "y").c_str()),
                                 atof(GetAttribute(cur_child_node, "z").c_str())));
        } else if (GetName(cur_child_node) == "bond") {
          shell->addBond(GetAttribute(cur_child_node, "first"),
                         GetAttribute(cur_child_node, "second"),
                         atoi(GetAttribute(cur_child_node, "dof").c_str()));
        }
      }

      PResources::AddBlockShell(shell);
    }
  }

  xmlFreeDoc(doc);
}

void PInit::setupBlockConnections(const string &connectionsFile) {
  StringPair connectionPair;
  HASH_MAP_STRPAIR_EX(Vector3) relativePos;

  xmlDocPtr doc = xmlReadFile(connectionsFile.c_str(), NULL, 0);
  if (doc == NULL) {
    PUtilities::AbortProgram("Could not read " + connectionsFile);
  }

  xmlNodePtr start_node = doc->children->children;
  for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
    if (GetName(cur_node) == "connection") {
      string first_block = GetAttribute(cur_node, "first"),
             second_block = GetAttribute(cur_node, "second");
      connectionPair = make_pair(first_block, second_block);

      for (xmlNodePtr cur_child_node = cur_node->children;
                      cur_child_node != NULL;
                      cur_child_node = cur_child_node->next) {
        if (GetName(cur_child_node) == "bond") {
          PBlockBondShell curBondShell = PBlockBondShell(GetAttribute(cur_child_node, "first"),
                                                         GetAttribute(cur_child_node, "second"),
                                                         atoi(GetAttribute(cur_child_node, "dof").c_str()));
          if (blockBonds.find(connectionPair) == blockBonds.end()) {
            blockBonds[connectionPair] = vector<PBlockBondShell>(1, curBondShell);
          } else {
            blockBonds[connectionPair].push_back(curBondShell);
          }
        } else if (GetName(cur_child_node) == "relative_pos") {
          StringPair atomPair = make_pair(GetAttribute(cur_child_node, "from"),
                                          GetAttribute(cur_child_node, "to"));
          relativePos[atomPair] = Vector3(atof(GetAttribute(cur_child_node, "x").c_str()),
                                          atof(GetAttribute(cur_child_node, "y").c_str()),
                                          atof(GetAttribute(cur_child_node, "z").c_str()));
        }
      }

      relPos[connectionPair] = relativePos;
      if (blockBonds.find(connectionPair) != blockBonds.end()) {
        PResources::AddBlockConnection(new PBlockConnection(connectionPair.first, connectionPair.second,
          relativePos, blockBonds[connectionPair]));
      } else {
        LoopTK::DisplayWarning(connectionPair.first + " -> " + connectionPair.second + " has no bonds defined.");
        PResources::AddBlockConnection(new PBlockConnection(connectionPair.first, connectionPair.second,
          relativePos, vector<PBlockBondShell>()));
      }
    }
  }

  xmlFreeDoc(doc);
}

void PInit::setupResidueShells(const string &residuesFile)
{
  xmlDocPtr doc = xmlReadFile(residuesFile.c_str(), NULL, 0);
  if (doc == NULL) {
    PUtilities::AbortProgram("Could not read " + residuesFile);
  }

  xmlNodePtr start_node = doc->children->children;
  for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
    if (GetName(cur_node) == "residue") {
      // Initialize a new PBlockShell with this block name and type.
      string residue_name = GetAttribute(cur_node, "name"),
             start_atom = GetAttribute(cur_node, "start_atom");
      PResidueShell *shell;

      for (xmlNodePtr cur_child_node = cur_node->children;
                      cur_child_node != NULL;
                      cur_child_node = cur_child_node->next) {
        if (GetName(cur_child_node) == "block") {
          if (GetAttribute(cur_child_node, "core") == "true") {
            shell = new PResidueShell(residue_name, GetAttribute(cur_child_node, "id"),
                                      GetAttribute(cur_child_node, "name"), start_atom);
          } else {
            shell->addBlock(GetAttribute(cur_child_node, "id"),
                            GetAttribute(cur_child_node, "name"));
          }
        } else if (GetName(cur_child_node) == "connection") {
          shell->addBlockConnection(GetAttribute(cur_child_node, "first"),
                                    GetAttribute(cur_child_node, "second"));
        }
      }

      PResources::AddResidueShell(shell);
    }
  }

  xmlFreeDoc(doc);
}

void PInit::setupChiIndices(const string &chiIndexFile)
{
  xmlDocPtr doc = xmlReadFile(chiIndexFile.c_str(), NULL, 0);
  if (doc == NULL) {
    PUtilities::AbortProgram("Could not read " + chiIndexFile);
  }

  xmlNodePtr start_node = doc->children->children;
  for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
    if (GetName(cur_node) == "residue") {
      string residue_name = GetAttribute(cur_node, "name");

      for (xmlNodePtr cur_child_node = cur_node->children;
                      cur_child_node != NULL;
                      cur_child_node = cur_child_node->next) {
        if (GetName(cur_child_node) == "chi") {
          int index = atoi(GetAttribute(cur_child_node, "index").c_str());

          vector<string> chi_atoms;

          chi_atoms.push_back(GetAttribute(cur_child_node, "first"));
          chi_atoms.push_back(GetAttribute(cur_child_node, "second"));
          chi_atoms.push_back(GetAttribute(cur_child_node, "third"));
          chi_atoms.push_back(GetAttribute(cur_child_node, "fourth"));

          PResources::AddChiIndex(residue_name, index, chi_atoms);
        }
      }
    }
  }

  xmlFreeDoc(doc);
}

void PInit::setupEpsilonVals(const string &epsilonFile) {
  xmlDocPtr doc = xmlReadFile(epsilonFile.c_str(), NULL, 0);
  if (doc == NULL) {
    PUtilities::AbortProgram("Could not read " + epsilonFile);
  }

  xmlNodePtr start_node = doc->children->children;
  for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
    if (GetName(cur_node) == "atom_pair") {
      string first_elem = "", second_elem = "", epsilon_val = "";
      
      for (xmlNodePtr cur_child_node = cur_node->children;
                      cur_child_node != NULL;
                      cur_child_node = cur_child_node->next) {

        if (GetName(cur_child_node) == "first") {
          first_elem = reinterpret_cast<const char *>(cur_child_node->children->content);
        } else if (GetName(cur_child_node) == "second") {
          second_elem = reinterpret_cast<const char *>(cur_child_node->children->content);
        } else if (GetName(cur_child_node) == "epsilon") {
          epsilon_val = reinterpret_cast<const char *>(cur_child_node->children->content);
        }
      }

      PResources::AddEpsilonValue(make_pair(first_elem, second_elem),
                                  atof(epsilon_val.c_str()));
    }
  }

  xmlFreeDoc(doc);
}

void PInit::setupIDMaps(const string &mapsFile)
{
  xmlDocPtr doc = xmlReadFile(mapsFile.c_str(), NULL, 0);
  if (doc == NULL) {
    PUtilities::AbortProgram("Could not read " + mapsFile);
  }

  xmlNodePtr start_node = doc->children->children;
  for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
    if (GetName(cur_node) == "mapping") {
      string residue = "", from_name = "", to_name = "";
      
      for (xmlNodePtr cur_child_node = cur_node->children;
                      cur_child_node != NULL;
                      cur_child_node = cur_child_node->next) {

        if (GetName(cur_child_node) == "residue") {
          residue = reinterpret_cast<const char *>(cur_child_node->children->content);
        } else if (GetName(cur_child_node) == "from_name") {
          from_name = reinterpret_cast<const char *>(cur_child_node->children->content);
        } else if (GetName(cur_child_node) == "to_name") {
          to_name = reinterpret_cast<const char *>(cur_child_node->children->content);
        }
      }

      PResources::AddAtomIDMapping(residue, from_name, to_name);
    }
  }

  xmlFreeDoc(doc);
}

void PInit::setupRotamers(const string &rotamersFile)
{
  xmlDocPtr doc = xmlReadFile(rotamersFile.c_str(), NULL, 0);
  if (doc == NULL) {
    PUtilities::AbortProgram("Could not read " + rotamersFile);
  }

  xmlNodePtr start_node = doc->children->children;
  for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
    if (GetName(cur_node) == "rotamer") {
      string residue = "";
      vector<string> chi_degrees;
      
      for (xmlNodePtr cur_child_node = cur_node->children;
                      cur_child_node != NULL;
                      cur_child_node = cur_child_node->next) {

        if (GetName(cur_child_node) == "residue") {
          residue = reinterpret_cast<const char *>(cur_child_node->children->content);
        } else if (GetName(cur_child_node) == "chi") {
          chi_degrees.push_back(reinterpret_cast<const char *>(cur_child_node->children->content));
        }
      }

      PResources::AddRotamer(residue, chi_degrees);
    }
  }
}
