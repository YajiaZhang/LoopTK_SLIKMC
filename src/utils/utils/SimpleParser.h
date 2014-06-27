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

#ifndef UTILS_SIMPLE_PARSER_H
#define UTILS_SIMPLE_PARSER_H

#include <iostream>
#include <assert.h>
using namespace std;

//simple class to parse files (actually, to lex them).
//comment lines: whitespace, followed by comment char, up to endline
//rest of lines are broken up into tokens and punctuation, possibly separated by whitespace
class SimpleParser
{
public:
  enum Result { Continue, Stop, Error };
  SimpleParser(istream& in);
  virtual ~SimpleParser() {}
  virtual bool IsSpace(char c) const { return isspace(c); }
  virtual bool IsComment(char c) const { return c=='#'; }
  virtual bool IsToken(char c) const { return isalnum(c); }
  virtual bool IsPunct(char c) const { return !IsSpace(c) && !IsComment(c) && !IsToken(c); }
  virtual Result InputToken(const string& word)=0;
  virtual Result InputPunct(const string& punct)=0;
  virtual Result InputEndLine()=0;

  bool Read();

  //called by sub classes on error, stop
  bool EatSpace();
  bool ReadLine(string& str);

  istream& in;
  int lineno;
};

#endif
