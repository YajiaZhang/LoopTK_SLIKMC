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

#include "IntervalSet.h"
#include <algorithm>
#include <iostream>
#include <errors.h>
using namespace Math;
using namespace std;

struct LeftPointLess
{
  bool operator() (const Interval& a,const Interval& b) const
  { return a.a < b.a; }
};

void OpenIntervalSet::Union(const BaseT& set)
{
  //sort all the segments
  BaseT x=*this,y=set;
  LeftPointLess cmp;
  sort(x.begin(),x.end(),cmp);
  sort(y.begin(),y.end(),cmp);
  for(size_t i=0;i+1<x.size();i++)
    Assert(x[i].b <= x[i+1].a);
  for(size_t i=0;i+1<y.size();i++)
    Assert(y[i].b <= y[i+1].a);

  size_t i=0,j=0;
  OpenInterval temp;
  resize(0);
  while(i != x.size() && j != y.size()) {
    if(x[i].a < y[j].a) {    //shift y intervals
      while(j != y.size() && y[j].b <= x[i].b) j++;
      temp.a = x[i].a;
      if(j != y.size() && y[j].a < x[i].b) {
	temp.b = y[j].b;
	j++;
      }
      else {
	temp.b=x[i].b;
      }
      i++;
    }
    else {   //shift x intervals
      while(i != x.size() && x[i].b <= y[j].b) i++;
      temp.a = y[j].a;
      if(i != x.size() && x[i].a < y[j].b) {
	temp.b = x[i].b;
	i++;
      }
      else {
	temp.b=y[j].b;
      }
      j++;
    }
    push_back(temp);
  }
  //one of the two (or both) are at the end
  for(;j<y.size();j++)
    push_back(y[j]);
  for(;i<x.size();i++)
    push_back(x[i]);
}

/*
  void Intersect(const BaseT&);
//  void Subtract(const ClosedBaseT&);
  void Union(const OpenInterval&);
*/

void OpenIntervalSet::Intersect(const OpenInterval& interval)
{
  for(size_t i=0;i<size();i++) {
    (*this)[i].setIntersection((*this)[i],interval);
    if((*this)[i].isEmpty()) {
      erase(begin()+i);
      i--;
    }
  }
}

//  void Subtract(const ClosedInterval&);


void ClosedIntervalSet::Union(const BaseT& set)
{
  //sort all the segments
  BaseT x=*this,y=set;
  LeftPointLess cmp;
  sort(x.begin(),x.end(),cmp);
  sort(y.begin(),y.end(),cmp);
  for(size_t i=0;i+1<x.size();i++)
    Assert(x[i].b < x[i+1].a);
  for(size_t i=0;i+1<y.size();i++)
    Assert(y[i].b < y[i+1].a);

  size_t i=0,j=0;
  ClosedInterval temp;
  resize(0);
  while(i != x.size() && j != y.size()) {
    if(x[i].a < y[j].a) {  //shift y intervals
      while(j != y.size() && y[j].b <= x[i].b) j++;
      temp.a = x[i].a;
      if(j != y.size() && y[j].a <= x[i].b) {
	temp.b = y[j].b;
	j++;
      }
      else {
	temp.b=x[i].b;
      }
      i++;
    }
    else {   //shift x intervals
      while(i != x.size() && x[i].b <= y[j].b) i++;
      temp.a = y[j].a;
      if(i != x.size() && x[i].a <= y[j].b) {
	temp.b = x[i].b;
	i++;
      }
      else {
	temp.b=y[j].b;
      }
      j++;
    }
    push_back(temp);
  }
  //one of the two (or both) are at the end
  for(;j<y.size();j++)
    push_back(y[j]);
  for(;i<x.size();i++)
    push_back(x[i]);
}

void ClosedIntervalSet::Intersect(const BaseT& set)
{
  //sort all the segments
  BaseT x=*this,y=set;
  LeftPointLess cmp;
  sort(x.begin(),x.end(),cmp);
  sort(y.begin(),y.end(),cmp);
  for(size_t i=0;i+1<x.size();i++)
    Assert(x[i].b <= x[i+1].a);
  for(size_t i=0;i+1<y.size();i++)
    Assert(y[i].b <= y[i+1].a);

  /*
  cout<<"Intersecting intervals: "<<endl;
  for(size_t i=0;i<x.size();i++)
    cout<<"["<<x[i].a<<","<<x[i].b<<"] ";
  cout<<endl;
  for(size_t i=0;i<y.size();i++)
    cout<<"["<<y[i].a<<","<<y[i].b<<"] ";
  cout<<endl;
  */

  size_t i=0,j=0;
  ClosedInterval temp;
  resize(0);
  while(i != x.size() && j != y.size()) {
    if(x[i].b < y[j].a) {  //shift x
      i++;
    }
    else if(y[j].b < x[i].a) {  //shift y
      j++;
    }
    else {
      //intersection
      temp.setIntersection(x[i],y[j]);
      Assert(!temp.isEmpty());
      push_back(temp);
      if(x[i].a < y[j].a) i++;
      else j++;
    }
  }

  /*
  cout<<"Res: "<<endl;
  for(size_t i=0;i<size();i++)
    cout<<"["<<(*this)[i].a<<","<<(*this)[i].b<<"] ";
  cout<<endl;
  */
}
/*
  void Subtract(const OpenBaseT&);
  void Union(const ClosedInterval&);
  void Intersect(const ClosedInterval&);
  void Subtract(const OpenInterval&);

*/
