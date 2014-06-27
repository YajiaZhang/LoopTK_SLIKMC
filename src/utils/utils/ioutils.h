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

#ifndef IO_UTILITIES_H
#define IO_UTILITIES_H

#include <iostream>
#include <string>
#include <vector>
#include <myfile.h>

/** @file utils/ioutils.h
 * @ingroup Utils
 * @brief Utilities for I/O.
 */

/** @addtogroup Utils */
/*@{*/

///Gets a quoted string from the istream into a char buffer (or string)
bool InputQuotedString(std::istream& in, char* str, int n);
bool InputQuotedString(std::istream& in, std::string&);
///Outputs "str".  Outputs \" for quote characters in str.
void OutputQuotedString(std::ostream& out, const char* str);
void OutputQuotedString(std::ostream& out, const std::string&);
///Returns true if str has a quote character
bool StringContainsQuote(char* str);
bool StringContainsQuote(const std::string& str);

///If c is preceded by a \, returns the translated ascii character.
///e.g. n->\n, t->\t, etc.
int TranslateEscape(int c);


///ReadFile() for STL vectors.  See myfile.h
template <class type>
bool ReadFile(File& f, std::vector<type>& v)
{
	size_t size;
	if(!ReadFile(f,size)) return false;
	v.resize(size);
	for(size_t i=0;i<size;i++)
		if(!ReadFile(f,v[i])) return false;
	return true;
}

///WriteFile() for STL vectors.  See myfile.h
template <class type>
bool WriteFile(File& f, const std::vector<type>& v)
{
	size_t size=v.size();
	if(!WriteFile(f,size)) return false;
	for(size_t i=0;i<size;i++)
		if(!WriteFile(f,v[i])) return false;
	return true;
}

///Inputs a vector from an iostream using the format 
///size  v[0] v[1] ... v[size-1]
template <class type>
bool InputVector(std::istream& in, std::vector<type>& v)
{
	size_t size;
	in>>size;
	if(in.bad()) return false;
	v.resize(size);
	for(size_t i=0;i<size;i++) {
	  in>>v[i];
	  if(in.bad()) return false;
	}
	return true;
}

///Outputs a vector to an iostream using the format 
///size  v[0] v[1] ... v[size-1]
template <class type>
bool OutputVector(std::ostream& out, const std::vector<type>& v)
{
  out<<v.size()<<'\t';
  for(size_t i=0;i<v.size();i++) {
    out<<v[i]<<" ";
  }
  return true;
}

/*@}*/

#endif
