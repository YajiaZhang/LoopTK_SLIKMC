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

#include "ProgressPrinter.h"
#include "../math/math.h"
#include <math.h>
using namespace std;

ProgressPrinter::ProgressPrinter(ostream& _out,int _max,int _increments)
  :out(_out),max(_max),increments(_increments),iter(0)
{
}

ProgressPrinter::ProgressPrinter(int _max,int _increments)
  :out(cout),max(_max),increments(_increments),iter(0)
{
}

void ProgressPrinter::Update()
{
  Update(iter+1);
}

void ProgressPrinter::Update(int _iter)
{
  iter=_iter;
  if((iter*increments) / max != ((iter-1)*increments) / max) {
    Print(float(iter)/float(max));
  }
}

void ProgressPrinter::Done()
{
  iter=0;
  out<<"100%"<<endl;
}

void ProgressPrinter::Print(float fraction)
{
  int n=int(floorf(fraction*100.0f+0.5f));
  out<<n<<"%..."; out.flush();
}

