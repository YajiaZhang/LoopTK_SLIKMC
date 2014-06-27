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

#include "ioutils.h"

//if c is preceded by a \, returns the translated ascii character
int TranslateEscape(int c)
{
  switch(c) {
  case 'a':  return '\a';
  case 'b':  return '\b';
  case 'n':  return '\n';
  case 'r':  return '\r';
  case 't':  return '\t';
  case 'v':  return '\v';
  default:   return c;
  }

}

bool InputQuotedString(std::istream& in, char* out, int n)
{
  int c;
  int state=0;  //0 = not yet read begin ", 1 is reading chars until end "
  int i=0;
  while((c=in.peek())!=EOF) {
    switch (state) {
    case 0:
      if(c=='\"')
	state = 1;
      else if(!isspace(c)) return false;
    case 1:
      if(c=='\"') {
	in.get();
	out[i]='\0';
	return true;
      }
      else if(c=='\\') {
	c=in.get();
	c=in.peek();
	out[i] = c;
      }
      else {
	if(i>=n) return false;
	out[i] = c;
	i++;
      }
      break;
    }
    c=in.get();
  }
  return false;
}

bool InputQuotedString(std::istream& in, std::string& out)
{
  int c;
  int state=0;
  while((c=in.peek())!=EOF) {
    switch (state) {
      case 0:
	if(c=='\"')
	  state = 1;
	else if(!isspace(c))
	  return false;
	break;
    case 1:
      if(c=='\"') { in.get(); return true; }
      else if(c=='\\') {
	c=in.get();
	c=in.peek();
	out+= c;
      }
      else
	out+= c;
      break;
    }
    c = in.get();
  }
  return false;
}

void OutputQuotedString(std::ostream& out, const char* str)
{
  if(StringContainsQuote(str)) {
    //output quotes using \"
    out<<'\"';
    const char* c=str;
    while(*c) {
      if(*c == '\"') out<<"\\\"";
      else out<<*c;
    }
    out<<'\"';
  }
  else
    out<<'\"'<<str<<'\"';
}

void OutputQuotedString(std::ostream& out, const std::string& str)
{
  OutputQuotedString(out,str.c_str());
}

bool StringContainsQuote(const char* str)
{
	return strchr(str,'\"')!=NULL;
}

bool StringContainsQuote(const std::string& str)
{
	return str.rfind('\"')!=str.npos;
}
