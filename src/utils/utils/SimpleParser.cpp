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

#include "SimpleParser.h"

SimpleParser::SimpleParser(istream& _in)
  :in(_in),lineno(1)
{}

bool SimpleParser::EatSpace()
{
  int c;
  while((c=in.peek()) != EOF) {
    if(!isspace(c))
      return true;
    c=in.get();
  }
  return true;
}

bool SimpleParser::ReadLine(string& str)
{
  //read until comment characters or eof
  str.erase();
  int c;
  while((c=in.peek()) != EOF) {
    if(!in) {
      cerr<<"Error while reading line!"<<endl;
      return false;
    }
    if(c == '\\') { //read literal (or skip endline)
      c=in.get();
      c=in.peek();
      if(c == '\r') {  //carriage return, this may be a DOS-style file
	c = in.get();
	c = in.peek();
      }
      if(c == EOF) {
	cerr<<"literal character \\ before end of file"<<endl;
	return false;
      }
      else if(c == '\n') { } //skip
      else str += c;
    }
    else if(c == '\n' || IsComment(c)) {
      return true;
    }
    else str += c;
    c=in.get();
  }
  cout<<"Reached end of file"<<endl;
  return true;
}

bool SimpleParser::Read()
{
  int c;
  int mode=0;  //0 whitespace/idle, 1 comment, 2 token, 3 punct
  string str;
  while((c=in.peek()) != EOF) {
    if(!in) {
      cerr<<"Error while reading characters!"<<endl;
      return false;
    }
    switch(mode) {
    case 0:
      if(IsSpace(c)) {  //continue
      }
      else if(IsComment(c)) mode=1;
      else if(IsToken(c)) { str+=c; mode=2;  }
      else if(IsPunct(c)) { str+=c; mode=3; }
      else {
	cerr<<"Illegal character "<<(char)c<<endl;
	return false;
      }
      break;
    case 1:
      if(c=='\n')
	mode=0;
      break;
    case 2:
      if(IsToken(c)) {str+=c;}
      else {
	Result res=InputToken(str);
	if(res == Stop) return true;
	else if(res == Error) {
	  cerr<<"Error on token "<<str<<endl;
	  return false;
	}
	str.erase();
	if(IsSpace(c)) mode=0;
	else if(IsComment(c)) mode=1;
	else if(IsPunct(c)) { str+=c; mode=3;	}
	else {
	  cerr<<"Illegal character "<<(char)c<<endl;
	  return false;
	}
      }
      break;
    case 3:
      if(IsPunct(c)) {str+=c;}
      else {
	Result res=InputPunct(str);
	if(res == Stop) return true;
	else if(res == Error) {
	  cerr<<"Error on token "<<str<<endl;
	  return false;
	}
	str.erase();
	if(IsSpace(c)) mode=0;
	else if(IsComment(c)) mode=1;
	else if(IsToken(c)) { str+=c; mode=2; }
	else {
	  cerr<<"Illegal character "<<(char)c<<endl;
	  return false;
	}
      }
      break;
    }
    if(c == '\n') lineno++;
    c = in.get();
  }
  in.get();
  assert(!in);
  return true;
}

