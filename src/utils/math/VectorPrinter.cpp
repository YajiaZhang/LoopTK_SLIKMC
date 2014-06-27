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

#include "VectorPrinter.h"
#include "complex.h"
#include <utils/stringutils.h>
#include "ASCIIShade.h"
using namespace std;

namespace Math {

VectorPrinter::VectorPrinter(const fVector& v)
  :fv(&v),dv(NULL),cv(NULL),delim(' '),bracket('['),mode(Normal)
{}

VectorPrinter::VectorPrinter(const dVector& v)
  :fv(NULL),dv(&v),cv(NULL),delim(' '),bracket('['),mode(Normal)
{}

VectorPrinter::VectorPrinter(const cVector& v)
  :fv(NULL),dv(NULL),cv(&v),delim(' '),bracket('['),mode(Normal)
{}

template<class T>
void PrintVector(const VectorTemplate<T>& x,ostream& out,char delim,char bracket)
{
  char closebracket = CloseBracket(bracket);
  if(bracket) out<<bracket;
  VectorIterator<T> v=x.begin();
  for(int i=0; i<x.n; i++,v++)
    out<<*v<<delim;
  if(bracket) out<<closebracket;
}

template <class T>
void OutputPlusMinus(ostream& out,const VectorTemplate<T>& x,T zeroTolerance=Epsilon)
{
  for(int i=0;i<x.n;i++) {
    if(x(i) < -zeroTolerance) out<<'-';
    else if(x(i) > zeroTolerance) out<<'+';
    else out<<'0';
  }
}

void VectorPrinter::Print(ostream& out) const
{
  switch(mode) {
  case Normal:
    if(fv) PrintVector(*fv,out,delim,bracket);
    else if(dv) PrintVector(*dv,out,delim,bracket);
    else if(cv) PrintVector(*cv,out,delim,bracket);
    break;
  case AsciiShade:
    if(fv) OutputASCIIShade(out,*fv);
    else if(dv) OutputASCIIShade(out,*dv);
    else if(cv) { cerr<<"Unable to output an ASCII-shaded complex matrix"<<endl; }
    break;
  case PlusMinus:
    if(fv) OutputPlusMinus(out,*fv);
    else if(dv) OutputPlusMinus(out,*dv);
    else if(cv) { cerr<<"Unable to output an +/- shaded complex matrix"<<endl; }
    break;
  }
}

ostream& operator << (ostream& out,const VectorPrinter& vp)
{
  vp.Print(out);
  return out;
}

} //namespace Math
