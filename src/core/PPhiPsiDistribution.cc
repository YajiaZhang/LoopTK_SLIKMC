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
#include "PPhiPsiDistribution.h"
using namespace std;

const string &ramachandranFile = "resources/ramachandran.xml";

PPhiPsiDistribution::PPhiPsiDistribution(const vector<vector<double> > &v, const string &name) {
  AA_Name = name;

  distribution = v;
  PhiIntervalNum = v.size();
  if (PhiIntervalNum != 0) {
    PsiIntervalNum = v[0].size();
    PhiIntervalSize = 360 / PhiIntervalNum;
    PsiIntervalSize = 360 / PsiIntervalNum;
  }
  else {
    PsiIntervalNum = 0;
    PhiIntervalSize = 360;
    PsiIntervalSize = 360;
  }

  // Compute unconditional phi distribution
  for (int i = 0; i < PhiIntervalNum; i++) {
    double sum = 0;
    for (int j = 0; j < PsiIntervalNum; j++) {
      sum += distribution[i][j];
    }
    PhiDistribution.push_back(sum);
  }

  // Compute unconditional psi distribution
  for (int i = 0; i < PsiIntervalNum; i++) {
    double sum = 0;
    for (int j = 0; j < PhiIntervalNum; j++) { 
      sum += distribution[j][i];
    }
    PsiDistribution.push_back(sum);
  }
}

void PPhiPsiDistribution::print(ostream &out) const {
  out << "Distribution for " << AA_Name << endl;
  for (int i=0; i<PhiIntervalNum; ++i) {
    out << i*PhiIntervalSize-180 << "~" << (i+1)*PhiIntervalSize-180-1 << "\t";
    for (int j=0; j<PsiIntervalNum; ++j) 
      out << distribution[i][j] << "\t";
    out << endl;
  }
  out << "\t";
  for (int i=0; i<PsiIntervalNum; ++i)
    out << i*PsiIntervalSize-180 << "~" << (i+1)*PsiIntervalSize-180-1 << "\t";
  out << endl;

  out << "PhiIntervalNum = " << PhiIntervalNum << "\t";
  out << "PsiIntervalNum = " << PsiIntervalNum << endl;
  out << "PhiIntervalSize = " << PhiIntervalSize << "\t";
  out << "PsiIntervalSize = " << PsiIntervalSize << endl;
}

double PPhiPsiDistribution::getProbPhiPsi(double phi, double psi) const {
  int phiIntervalIdx = getPhiIntervalIdx(phi);
  int psiIntervalIdx = getPsiIntervalIdx(psi);
  return distribution[phiIntervalIdx][psiIntervalIdx];
}

vector<double> PPhiPsiDistribution::getPhiDistribution(double psiValue) const {
  int psiIntervalIdx = getPsiIntervalIdx(psiValue);
  vector<double> phiDist;
  for (int i=0; i<PhiIntervalNum; ++i) {
    phiDist.push_back(distribution[i][psiIntervalIdx]);
  }
  PMath::normalize(phiDist);
  return phiDist;
}

vector<double> PPhiPsiDistribution::getPsiDistribution(double phiValue) const {
  int phiIntervalIdx = getPhiIntervalIdx(phiValue);
  vector<double> psiDist;
  for (int i=0; i<PsiIntervalNum; ++i) {
    psiDist.push_back(distribution[phiIntervalIdx][i]);
  }
  PMath::normalize(psiDist);
  return psiDist;
}

int PPhiPsiDistribution::getPhiIntervalIdx (double phiValue) const {
  return (int)(floor( (phiValue+180)/PhiIntervalSize ));
}

int PPhiPsiDistribution::getPsiIntervalIdx (double psiValue) const {
  return (int)(floor( (psiValue+180)/PsiIntervalSize ));
}

map<string,PPhiPsiDistribution> PPhiPsiDistribution::generateRamachandran () {
  map<string, PPhiPsiDistribution> ramachandran;

  xmlDocPtr doc = xmlReadFile(ramachandranFile.c_str(), NULL, 0);
  if (doc == NULL) {
    PUtilities::AbortProgram("Could not read " + ramachandranFile);
  }

  xmlNodePtr start_node = doc->children->children;
  for (xmlNodePtr cur_node = start_node; cur_node != NULL; cur_node = cur_node->next) {
    if (PInit::GetName(cur_node) == "distribution") {
      string residue = "";
      vector<vector<double> > distribution;
      
      for (xmlNodePtr cur_child_node = cur_node->children;
                      cur_child_node != NULL;
                      cur_child_node = cur_child_node->next) {

        if (PInit::GetName(cur_child_node) == "residue") {
          residue = reinterpret_cast<const char *>(cur_child_node->children->content);
        } else if (PInit::GetName(cur_child_node) == "phi_intervals") {
          string phi_intervals = reinterpret_cast<const char *>(cur_child_node->children->content);
          distribution.resize(atoi(phi_intervals.c_str()));
        } else if (PInit::GetName(cur_child_node) == "psi_intervals") {
          string psi_intervals = reinterpret_cast<const char *>(cur_child_node->children->content);
          for (unsigned i = 0; i < distribution.size(); i++) {
            distribution[i].resize(atoi(psi_intervals.c_str()));
          }
        } else if (PInit::GetName(cur_child_node) == "probability") {
          int phi = atoi(PInit::GetAttribute(cur_child_node, "phi_interval").c_str()),
              psi = atoi(PInit::GetAttribute(cur_child_node, "psi_interval").c_str());
          string value = reinterpret_cast<const char *>(cur_child_node->children->content);
          distribution[phi - 1][psi - 1] = atof(value.c_str());
        }
      }

      ramachandran[residue] = PPhiPsiDistribution(distribution, residue);
    }
  }

  xmlFreeDoc(doc);
  return ramachandran;
}
