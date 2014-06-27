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
#include "PExtension.h"

void CS2IO::appendConformation(const string &cs2File, PProtein *loop)
{
  int numConfsInFile = getNumConformations(cs2File);
  pair<int, int> topIndices, sanityCheck;

  if (loop == NULL) {
    PUtilities::AbortProgram("Error: Cannot append null loop conformation to file " + cs2File + "!");
  }

  if (numConfsInFile == 0) {
    PUtilities::AbortProgram("File with at least 1 conformation must exist first.");
  } else {
    /* Make sure this is the same loop as is already in the file! */
    topIndices = loop->getTopLevelIndices();
    sanityCheck = parseLoopLine(cs2File, PUtilities::getLinesStarting(cs2File, "loop ", "", true)[0]);

    if (!(topIndices == sanityCheck)) {
      PUtilities::AbortProgram("Error: Loop to be appended to file " + cs2File +
        " does not start at residue " + PUtilities::toStr(sanityCheck.first) +
        " and/or end at residue " + PUtilities::toStr(sanityCheck.second) +
        " in the top-level chain!");
    }

    /* Output loop as the next conformation number. */
    outputHelper(cs2File, loop, numConfsInFile + 1, topIndices.first, topIndices.second);
  }
}

PProtein* CS2IO::readConformation(const string &cs2File, int confNum)
{
  vector<string> atomLines, loopLines, confLines;
  parseFile(cs2File, confNum, atomLines, loopLines, confLines);

  PProtein *protein = PDBIO::readFromLines(atomLines);
  pair<int, int> sanityCheck = parseLoopLine(cs2File, loopLines[0]);

  /* If a valid protein was read, apply the appropriate changes. */
  if (protein != NULL) {
    for(unsigned i = 0; i < confLines.size(); i++) {
      csConfLine curConfLine = parseConfLine(cs2File, confLines[i], sanityCheck);
      applyConfLine(protein, curConfLine);
    }
  } else {
    PUtilities::AbortProgram("Error: Couldn't read conformation " + PUtilities::toStr(confNum) +
        " in file " + cs2File + " (check file integrity?)");
  }

  return new PProtein(protein, sanityCheck.first, sanityCheck.second);
}

pair<PLightChain *, PLightChain *> CS2IO::getNonLoopChains(const string &cs2File)
{
  vector<string> atomLines, loopLines;
  parseFile(cs2File, atomLines, loopLines);

  return getNonLoopChains(atomLines, loopLines, cs2File);
}

pair<PLightChain *, PLightChain *> CS2IO::getNonLoopChains(const vector<string> &atomLines, const vector<string> &loopLines,
              const string &fileName)
{
  vector<PDBIO::pdbAtomLine> parsedLines;
  PProtein *head = new PProtein(), *tail = new PProtein();

  pair<int, int> loopIndices = parseLoopLine(fileName, loopLines[0]);

  for(unsigned i = 0; i < atomLines.size(); i++) {
    PDBIO::pdbAtomLine curLine = PDBIO::parseAtomLine(atomLines[i]);
    if (!parsedLines.empty() && curLine.resNum != parsedLines[0].resNum) {
      PResidueSpec curSpec = PDBIO::getSpec(parsedLines);

      if (parsedLines[0].resNum < loopIndices.first + 1) {
        head->AddResidue(parsedLines[0].resName, curSpec);
      } else if (parsedLines[0].resNum > loopIndices.second + 1) {
        tail->AddResidue(parsedLines[0].resName, curSpec);
      }

      parsedLines.clear();
    }
    parsedLines.push_back(curLine);
  }
  PResidueSpec prs = PDBIO::getSpec(parsedLines);
  tail->AddResidue(parsedLines[0].resName, prs);

  return make_pair(head, tail);
}

PLightChain* CS2IO::readConformationLoopOnly(const string &cs2File, int confNum)
{
  vector<string> atomLines, loopLines, confLines;
  parseFile(cs2File, confNum, atomLines, loopLines, confLines);

  return readConformationLoopOnly(atomLines, loopLines, confLines, cs2File);
}

PLightChain* CS2IO::readConformationLoopOnly(const vector<string> &atomLines, const vector<string> &loopLines,
            const vector<string> &confLines, const string &fileName)
{
  vector<PDBIO::pdbAtomLine> parsedLines;
  PProtein *loop = new PProtein();
  pair<int, int> loopIndices = parseLoopLine(fileName, loopLines[0]);

  for(unsigned i = 0; i < atomLines.size(); i++) {
    PDBIO::pdbAtomLine curLine = PDBIO::parseAtomLine(atomLines[i]);
    if (!parsedLines.empty() && curLine.resNum != parsedLines[0].resNum) {
      PResidueSpec curSpec = PDBIO::getSpec(parsedLines);

      if (parsedLines[0].resNum >= loopIndices.first + 1 &&
         parsedLines[0].resNum <= loopIndices.second + 1) {
        loop->AddResidue(parsedLines[0].resName, curSpec);
      }

      parsedLines.clear();
    }
    parsedLines.push_back(curLine);
  }

  for(unsigned i = 0; i < confLines.size(); i++) {
    csConfLine curConfLine = parseConfLine(fileName, confLines[i], loopIndices);
    applyConfLine(loop, curConfLine, loopIndices.first);
  }

  return loop;
}

PConformationSpace* CS2IO::readConformationSpace(const string &cs2File)
{
  vector<vector<string> > confMap;
  vector<string> atomLines, loopLines;

  parseFile(cs2File, atomLines, loopLines, confMap);

  PConformationSpace *result = new PConformationSpace(getNonLoopChains(atomLines, loopLines, cs2File));

  for(unsigned i = 0; i < confMap.size(); i++) {
    result->push_back(readConformationLoopOnly(atomLines, loopLines, confMap[i], cs2File));
  }

  return result;
}

void CS2IO::writeToFile(const string &cs2File, PProtein *loop)
{
  PChain *topLevel = loop->getTopLevelChain();
  pair<int, int> topIndices = loop->getTopLevelIndices();

  /* Output the top-level protein, i.e. the "parent" of loop. */
  PDBIO::writeToFile(topLevel, cs2File, "pdb ");

  /* Output the "loop" line indicating loop's residue indices. */
  PUtilities::appendToFile(cs2File, "loop " + PUtilities::toStr(topIndices.first) +
         " " + PUtilities::toStr(topIndices.second));

  /* Output loop as the first conformation in the file. */
  outputHelper(cs2File, loop, 1, topIndices.first, topIndices.second);
}

void CS2IO::createPDBDirectory(const string &cs2File, const string &directoryName)
{
  int num = getNumConformations(cs2File);
  PProtein *curProtein;

  PUtilities::makeDirectory(directoryName);
  for(int i = 1; i <= num; i++) {
    curProtein = readConformation(cs2File, i);
    PDBIO::writeToFile(curProtein->getTopLevelProtein(), directoryName + 
           "/conf" + PUtilities::toStr(i) + ".pdb");
    delete curProtein;
  }
}

/* Private methods */

csConfLine CS2IO::parseConfLine(const string &cs2File, const string &cs2Line, const pair<int, int> &sanityCheck)
{
  vector<string> lineTokens = PUtilities::getTokens(cs2Line);
  if (lineTokens.size() != 5) {
    PUtilities::AbortProgram("Error: Malformed conformation line: [" + cs2Line +
      "] in CS2 file " + cs2File + "!");
  }

  int resNum = atoi(lineTokens[0].c_str());
  if (resNum < sanityCheck.first || resNum > sanityCheck.second) {
    PUtilities::AbortProgram("Error: Residue index in line: [" + cs2Line +
      "] in CS2 file " + cs2File + " does not match loop indices!");
  }

  string atomID = lineTokens[1];
  Vector3 pos;

  pos.x = atof(lineTokens[2].c_str());
  pos.y = atof(lineTokens[3].c_str());
  pos.z = atof(lineTokens[4].c_str());

  return make_pair(make_pair(resNum, atomID), pos);
}

pair<int, int> CS2IO::parseLoopLine(const string &cs2File, const string &loopLine)
{
  vector<string> lineTokens = PUtilities::getTokens(loopLine);
  int resStart, resEnd;

  if (lineTokens.size() != 2) {
    PUtilities::AbortProgram("Error: Malformed \"loop\" line: [" + loopLine +
      "] in CS2 file " + cs2File + "!");
  }

  /* Parse the start and end residue indices. */
  resStart = atoi(lineTokens[0].c_str());
  resEnd = atoi(lineTokens[1].c_str());

  if (resStart < 0 || resEnd < 0 || resStart >= resEnd) {
    PUtilities::AbortProgram("Error: \"loop\" line [" + loopLine + "] in CS2 file " +
      cs2File + " contains invalid residue indices!");
  }

  return make_pair(resStart, resEnd);
}

void CS2IO::parseFile(const string &cs2File, vector<string> &atomLines, vector<string> &loopLines)
{
  vector<string> fileLines = PUtilities::getLines(cs2File);

  atomLines = extractAtomLines(fileLines, cs2File);
  loopLines = extractLoopLines(fileLines, cs2File);
}

void CS2IO::parseFile(const string &cs2File, vector<string> &atomLines,
                      vector<string> &loopLines, vector<vector<string> > &allConfLines)
{
  vector<string> fileLines = PUtilities::getLines(cs2File);

  atomLines = extractAtomLines(fileLines, cs2File);
  loopLines = extractLoopLines(fileLines, cs2File);
  allConfLines = extractConfLines(fileLines, cs2File);
}

void CS2IO::parseFile(const string &cs2File, int confNum, vector<string> &atomLines,
      vector<string> &loopLines, vector<string> &confLines)
{
  vector<vector<string> > confMap;

  parseFile(cs2File, atomLines, loopLines, confMap);
  confLines = confMap[confNum - 1];
}

int CS2IO::getNumConformations(const string &cs2File)
{
  vector<string> confLines = PUtilities::getLinesStarting(cs2File, "cn", "", true), curTokens;
  int maxConfNum = 0;

  for(unsigned i = 0; i < confLines.size(); i++) {
    curTokens = PUtilities::getTokens(confLines[i]);

    if (atoi(curTokens[0].c_str()) > maxConfNum) {
      maxConfNum = atoi(curTokens[0].c_str());
    }
  }

  return maxConfNum;
}

void CS2IO::applyConfLine(PProtein* &loop, const csConfLine &conf, int offset)
{
  if (loop != NULL) {
    pair<int, string> resAtomPair = conf.first;
    Vector3 newPos = conf.second;

    PAtom *targetAtom = loop->getAtomAtRes(resAtomPair.second, resAtomPair.first - offset);
    targetAtom->changePosition(newPos);
  }
}

void CS2IO::outputHelper(const string &cs2File, PProtein *loop, int confNum,
       int startIndex, int endIndex)
{
  string curLine;
  HASH_MAP_STR(PAtom *) *atomMap;

  for(int i = 0; i < loop->size(); i++) {
    PResidue *curResidue = loop->getResidue(i);
    atomMap = curResidue->getAtomMap();

    for(HASH_MAP_STR(PAtom *)::const_iterator it = atomMap->begin(); it != atomMap->end(); ++it) {
      /* First token in output consists of the conformation number. */
      curLine = "cn" + PUtilities::toStr(confNum) + " ";
      
      /* Second token is the residue number; third is atom ID. */
      curLine += PUtilities::toStr(startIndex + i) + " ";
      curLine += it->first + " ";

      /* Remaining tokens are the (x, y, z) coordinates of this atom. */
      curLine += PUtilities::toStr(it->second->getPos().x) + " " +
           PUtilities::toStr(it->second->getPos().y) + " " +
           PUtilities::toStr(it->second->getPos().z) + " ";

      /* Write the line to file. */
      PUtilities::appendToFile(cs2File, curLine);
    }
  }
}

vector<string> CS2IO::extractAtomLines(const vector<string> &cs2Lines, const string &fileName)
{
  vector<string> pdbLines  = PUtilities::getLinesStarting(cs2Lines, "pdb ", "", true),
                 atomLines = PUtilities::getLinesStarting(pdbLines, "ATOM  ", "", false);

  if (atomLines.size() == 0) {
    PUtilities::AbortProgram("Error: CS2 file " + fileName + " does not contain PDB data!");
  }

  return atomLines;
}

vector<string> CS2IO::extractLoopLines(const vector<string> &cs2Lines, const string &fileName)
{
  vector<string> loopLines = PUtilities::getLinesStarting(cs2Lines, "loop ", "", true);

  if (loopLines.size() != 1) {
    PUtilities::AbortProgram("Error: CS2 file " + fileName + " does not contain exactly one \"loop\" line!");
  }

  return loopLines;
}

vector<vector<string> > CS2IO::extractConfLines(const vector<string> &cs2Lines, const string &fileName)
{
  vector<vector<string> > result;
  vector<string> confLines = PUtilities::getLinesStarting(cs2Lines, "cn", "", true), curTokens;

  int curConfNum;
  string curConfLine;

  if (confLines.size() == 0) {
    PUtilities::AbortProgram("Error: CS2 file " + fileName + " contains no conformation data!");
  }

  for(unsigned i = 0; i < confLines.size(); i++) {
    curTokens = PUtilities::getTokens(confLines[i]);
    if (curTokens.size() != 6) {
      PUtilities::AbortProgram("Error: Malformed conformation line \"" + confLines[i] +
             "\" in CS2 file " + fileName);
    }

    curConfNum = atoi(curTokens[0].c_str());
    curConfLine = curTokens[1] + " " + curTokens[2] + " " + curTokens[3] + " " +
                  curTokens[4] + " " + curTokens[5] + " ";

    if (result.size() < curConfNum) {
      result.resize(result.size() + 1);
    }
    result[curConfNum - 1].push_back(curConfLine);
  }

  return result;
}
