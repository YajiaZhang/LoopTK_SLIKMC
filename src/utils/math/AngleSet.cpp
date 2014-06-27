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

#include "AngleSet.h"
using namespace Math;
using namespace std;

void AngleSet::Intersect(const vector<AngleInterval>& s)
{
  AngleSet newIntervals;
  AngleInterval temp;
  for(size_t i=0; i<size(); i++) {
    for(size_t j=0;j<s.size();j++) {
      temp.setIntersection(operator[](i),s[j]);
      if(!temp.isEmpty())
	newIntervals.push_back(temp);
    }
  }
  ::swap(*this,newIntervals);
}

//void AngleSet::Union(const AngleInterval& r);

void AngleSet::Intersect(const AngleInterval& r)
{
  AngleSet newIntervals;
  AngleInterval temp;
  for(size_t i=0; i<size(); i++) {
    temp.setIntersection(operator[](i),r);
    if(!temp.isEmpty())
      newIntervals.push_back(temp);
  }
  ::swap(*this,newIntervals);
}

//void AngleSet::Subtract(const AngleInterval& r);

