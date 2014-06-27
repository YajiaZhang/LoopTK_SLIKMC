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

#include "PUtilities.h"

void PUtilities::copyFile(const string &newFileName, const string &fileToCopy)
{
  char *buffer;
  long size;

  ifstream infile(fileToCopy.c_str(), ifstream::binary);
  ofstream outfile(newFileName.c_str(), ofstream::binary);

  infile.seekg(0,ifstream::end);
  size = infile.tellg();    /* Get size of input file. */
  infile.seekg(0);

  buffer = new char [size];  /* Allocate input buffer   */
  infile.read(buffer,size);  /* Read input to buffer.   */

  outfile.write(buffer,size);  /* Write output to buffer. */
  
  delete[] buffer;    /* Free dynamic memory.    */

  outfile.close();
  infile.close();
}

void PUtilities::makeDirectory(const string &dirName)
{
  if (mkdir(dirName.c_str(), S_IRWXU) != 0) {
    AbortProgram("Error: Could not create directory " + dirName + ".");
  }
}

/*
 * String utilities.
 */

string PUtilities::padLeft(const string &s, unsigned int len, char padChar)
{
  string sCopy(s), pad(1, padChar);
  while(sCopy.size() < len) {
    sCopy = pad + sCopy;
  }

  return sCopy;
}

string PUtilities::padRight(const string &s, unsigned int len, char padChar)
{
  string sCopy(s), pad(1, padChar);
  while(sCopy.size() < len) {
    sCopy += pad;
  }

  return sCopy;
}

void PUtilities::trimSpaces(string& str)
{
  string::size_type pos = str.find_last_not_of(' ');

  if (pos != string::npos) {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if (pos != string::npos) str.erase(0, pos);
  } else {
    str.erase(str.begin(), str.end());
  }
}

vector<string> PUtilities::getLinesStarting(const vector<string> &lines, const string &prefix,
        const string &chainId, bool trimPrefix, const string &terminationCode)
{
  vector<string> v;
  string curString;

  for(unsigned i = 0; i < lines.size(); i++) {
    curString = lines[i];

    if (chainId!="" && curString.substr(21,1)!=chainId)
      continue;

    if (terminationCode != "" && curString.find(terminationCode) == 0) {
      break;
    }

    if (curString.find(prefix) == 0) {
      if (trimPrefix) {
        curString.replace(0, prefix.length(), "");
      }

                  v.push_back(curString);  // curString starts with prefix
    }
  }

  return v;    
}

vector<string> PUtilities::getTokens(const string &curLine)
{
        vector<string> v;
        stringstream ss(curLine);
        string curString;
 
        while(1) {
                ss >> curString;
                if (ss.fail()) break;
 
                v.push_back(curString);
        }
 
        return v;
}

/*
 * Text I/O utilities.
 */

vector<string> PUtilities::getLinesStarting(const string &fileName, const string &prefix,
        const string &chainId, bool trimPrefix, const string &terminationCode)
{
        vector<string> fileLines;
        ifstream inputFile(fileName.c_str());
        string curString;
 
        while(1) {
        	getline(inputFile, curString);
        	if (inputFile.fail())
        		break;
 
        	fileLines.push_back(curString);
        }
 
        return getLinesStarting(fileLines, prefix, chainId, trimPrefix, terminationCode);
}

vector<string> PUtilities::getLines(const string &fileName, const string &chainId)
{
	return getLinesStarting(fileName, "", chainId);
}
 
void PUtilities::appendToFile(const string &fileName, const string &line)
{
  ofstream outputFile;

  outputFile.open(fileName.c_str(), ofstream::out | ofstream::app);
  outputFile << line << endl;
  outputFile.close();
}

void PUtilities::AbortProgram(const string &message)
{
   cerr << message << endl;
//  exit(1);
   abort();
}

const void* PUtilities::PointerThatIsNot(const void *p1, const void *p2, const void *have)
{
  if (p1 == have) return p2;
  else return p1;
}
