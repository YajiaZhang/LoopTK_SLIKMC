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

#include "ASCIIShade.h"
using namespace std;

namespace Math {

const static int kNumAsciiShades = 15;
const static char kAsciiShades[15] = {
  'W',
  'w',
  '#',
  '%',
  '&',
  '*',
  '+',
  ' ',  //0 is element 7
  '.',
  ':',
  'o',
  '0',
  '8',
  'O',
  '@',
};
  

char ASCIIShade(double x)
{
  if(IsNaN(x)) return 'E';
  if(IsInf(x)==1) return 'I';
  else if(IsInf(x)==-1) return 'i';
  int index = (int)Trunc(x*8) + 7;
  if(index < 0) index=0;
  if(index >= kNumAsciiShades) index=kNumAsciiShades-1;
  if(index == 7) {
    if(x > 0) return kAsciiShades[8];
    else if(x < 0) return kAsciiShades[6];
    else return kAsciiShades[7];
  }
  return kAsciiShades[index];
}

void OutputASCIIShade(ostream& out,double x)
{
  out<<ASCIIShade(x);
}

void OutputASCIIShade(ostream& out,const fVector& x,float scale)
{
  if(scale == 0) scale = x.maxAbsElement();
  out<<scale<<" x ";
  out<<'[';
  for(int i=0;i<x.n;i++)
    out<<ASCIIShade(x(i)/scale);
  out<<']';
}

void OutputASCIIShade(ostream& out,const dVector& x,double scale)
{
  if(scale == 0) scale = x.maxAbsElement();
  out<<scale<<" x ";
  out<<'[';
  for(int i=0;i<x.n;i++)
    out<<ASCIIShade(x(i)/scale);
  out<<']';
}

void OutputASCIIShade(ostream& out,const fMatrix& A,float scale)
{
  if(scale == 0) scale = A.maxAbsElement();
  out<<scale<<" x"<<endl;
  for(int i=0;i<A.m;i++) {
    out<<'[';
    for(int j=0;j<A.n;j++) 
      out<<ASCIIShade(A(i,j)/scale);
    out<<']';
    if(i+1 < A.m)
      out<<endl;
  }
}


void OutputASCIIShade(ostream& out,const dMatrix& A,double scale)
{
  if(scale == 0) scale = A.maxAbsElement();
  out<<scale<<" x"<<endl;
  for(int i=0;i<A.m;i++) {
    out<<'[';
    for(int j=0;j<A.n;j++) 
      out<<ASCIIShade(A(i,j)/scale);
    out<<']';
    if(i+1 < A.m)
      out<<endl;
  }
}


} //namespace Math
