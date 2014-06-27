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

#ifndef PUTILITIES_H
#define PUTILITIES_H

#include "PLibraries.h"

#define		isNAN(x)		!((x) <= 0 || (x) > 0)
#define		ends_with(str, x)	(string(str).find(string(x)) != string::npos && \
					 string(str).find(string(x)) == string(str).length() - \
					 string(x).length())
//@package Utilities
/**
 * Contains miscellaneous utility methods
 * for string manipulation, text file I/O,
 * system interaction, and so on.  These
 * methods are made public in case they are
 * valuable to users of <code>LoopTK</code>
 * but are generally intended for internal
 * use.
 */

class PUtilities {
  public:

// File/directory manipulation functions -- might not be portable
// to all operating systems (Windows?)

    /**
     * Copies <code>fileToCopy</code> to the output file
     * <code>newFileName</code>.
     */

    static void copyFile(const string &newFileName, const string &fileToCopy);

    /**
     * Creates the directory specified by <code>dirName</code>.
     * Aborts the program if the directory cannot be created.
     */

    static void makeDirectory(const string &dirName);

// String manipulation functions

    /**
     * Removes leading and trailing spaces
     * from <code>str</code>.
     */

    static void trimSpaces(string &str);

    /**
     * Returns a string version of <code>x</code> if possible,
     * i.e. if <code>stringstream::operator<< </code> supports
     * objects of type <code>T</code>.
     */

     template <typename T>
     static string toStr(const T &x)
     {
       stringstream ss;
       ss << x;
       return ss.str();
     }

    /**
     * Left-pads <code>s</code> with <code>padChar</code>
     * until its length is <code>len</code>.
     */

    static string padLeft(const string &s, unsigned int len, char padChar = ' ');

    /**
     * Right-pads <code>s</code> with <code>padChar</code>
     * until its length is <code>len</code>.
     */

    static string padRight(const string &s, unsigned int len, char padChar = ' ');

    /**
     * Returns a <code>vector</code> of all lines in the
     * specified <code>vector</code> beginning with the
     * specified <code>prefix</code>.
     *
     * Optionally, if <code>trimPrefix == true</code>, the lines
     * in the resulting vector will have the prefix removed. Also,
     * if <code>terminationCode</code> is not the empty string,
     * parsing will stop immediately if a line beginning with
     * <code>terminationCode</code> is found. 
     */

    static vector<string> getLinesStarting(const vector<string> &lines,
				const string &prefix, 
				const string &chainId = "",
				bool trimPrefix = false,
				const string &terminationCode = "");

    /**
     * Returns a vector of all the tokens
     * in the specified string.
     */

    static vector<string> getTokens(const string &curLine);

// Text I/O functions

    /**
     * Returns a <code>vector</code> of all lines in the
     * specified <code>fileName</code> beginning with the
     * specified <code>prefix</code>.
     *
     * Optionally, if <code>trimPrefix == true</code>, the lines
     * in the resulting vector will have the prefix removed. Also,
     * if <code>terminationCode</code> is not the empty string,
     * parsing will stop immediately if a line beginning with
     * <code>terminationCode</code> is found. 
     */

    static vector<string> getLinesStarting(const string &fileName,
				const string &prefix, 
				const string &chainId="",
				bool trimPrefix = false,
				const string &terminationCode = "");

    /**
     * Returns a vector of all lines in the
     * specified file, regardless of prefix.
     */

    static vector<string> getLines(const string &fileName, const string &chainId="");

    /**
     * Appends <code>line</code> to the file <code>fileName</code>.
     */

    static void appendToFile(const string &fileName, const string &line);

// Miscellaneous functions

    /**
     * Outputs <code>message</code> to <code>stderr</code>
     * and immediately terminates execution of the program.
     */

    static void AbortProgram(const string &message);

    /**
     * Returns <code>p2</code> if <code>p1</code> and <code>have</code>
     * point to the same address; otherwise, returns <code>p1</code>.
     */

    static const void *PointerThatIsNot(const void *p1, const void *p2, const void *have);

    /**
     * Creates a newly-allocated 2D table of <code>T</code>'s with
     * dimension <code>size</code> by <code>size</code>.
     */

    template <typename T>
    static T** New2DArray(unsigned size) {
      T** table = new T*[size];

      for(unsigned i = 0; i < size; i++) {
        table[i] = new T[size];
      }

      return table;
    }

    /**
     * Deletes the 2D array <code>table</code>, freeing all
     * dynamically-allocated memory associated with it.
     */

    template <typename T>
    static void Delete2DArray(T **table, unsigned size) {
      for(unsigned i = 0; i < size; i++) {
        delete[] table[i];
      }

      delete[] table;
    }

};

#endif
