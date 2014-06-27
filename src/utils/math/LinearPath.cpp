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

#include "LinearPath.h"
#include "metric.h"
#include <myfile.h>
#include <algorithm>
#include <errors.h>
using namespace Math;
using namespace std;

bool LessTime(const LinearPath::ControlPoint& a, Real b)
{
  return a.t < b;
}

void LinearPath::PreEval(Real t)
{
  //TODO: find the segment for t so we don't have to later
}

void LinearPath::Eval(Real t,Vector& x)
{
  vector<ControlPoint>::iterator prev=GetSegment(t);
  if(prev == --points.begin()) {
    x = points.front().x;
    return;    
  }
  Assert(t >= prev->t);
  vector<ControlPoint>::iterator next=prev; next++;
  if(next == points.end()) {
    x = prev->x;
    return;
  }
  Real dt = next->t-prev->t;
  if(dt <= Zero) {
    cout<<"LinearPath: Invalid range of times: ["<<prev->t<<","<<next->t<<"]"<<endl;
    cout<<"Size of path: "<<points.size()<<endl;
    cout<<"Time range: "<<BeginTime()<<" "<<EndTime()<<endl;
    cout<<"Query time: "<<t<<endl;
    Abort();
  }
  Real u = (t-prev->t)/dt;
  Interpolate(prev->x,next->x,u,x);
}

void LinearPath::Interpolate(const Vector& a,const Vector& b,Real u,Vector& x) const
{
  x.mul(a,(One-u));
  x.madd(b,u);
}

void LinearPath::Deriv(Real t,Vector& dx)
{
  vector<ControlPoint>::iterator prev=GetSegment(t);
  if(prev == --points.begin()) {
    dx.resize(points.front().x.n,Zero);
    return;    
  }
  Assert(t >= prev->t);
  vector<ControlPoint>::iterator next=prev; next++;
  if(next == points.end()) {
    dx.resize(points.front().x.n,Zero);
    return;
  }
  Real dt = next->t-prev->t;
  if(dt <= Zero) {
    cout<<"LinearPath: Invalid range of times: ["<<prev->t<<","<<next->t<<"]"<<endl;
    Abort();
  }
  Difference(next->x,prev->x,dx);
  dx /= dt;
}

void LinearPath::Difference(const Vector& a,const Vector& b,Vector& dx) const
{
  dx.sub(a,b);
}

Real LinearPath::Distance(const Vector& a,const Vector& b) const
{
  return Distance_L2(a,b);
}

vector<LinearPath::ControlPoint>::iterator LinearPath::GetSegment(Real t)
{
  return --std::lower_bound(points.begin(),points.end(),t,LessTime);
}

void LinearPath::SetSequence(const VectorSequence& s)
{
  points.resize(s.size());
  for(size_t i=0;i<s.size();i++) {
    points[i].t = (Real)i;
    points[i].x = s[i];
  }
}

void LinearPath::ArcLengthParameterize()
{
  if(points.empty()) return;
  points[0].t = Zero;
  for(size_t i=1;i<points.size();i++) {
    points[i].t = points[i-1].t + Distance(points[i].x,points[i-1].x);
  }
}

Real LinearPath::Length() const
{
  if(points.empty()) return Zero;
  Real len=Zero;
  for(size_t i=0;i+1<points.size();i++) {
    len += Distance(points[i].x,points[i+1].x);
  }
  return len;
}

void LinearPath::ScaleTime(Real s)
{
  for(size_t i=0;i<points.size();i++) points[i].t*=s;
}

void LinearPath::OffsetTime(Real off)
{
  for(size_t i=0;i<points.size();i++) points[i].t+=off;
}

void LinearPath::Concat(const LinearPath& p)
{
  size_t offset = points.size();
  Real toffset = (offset==0?Zero:EndTime());
  cout<<"Concat, offset index "<<offset<<endl;
  cout<<"     offset time  "<<toffset<<endl;
  cout<<"     new points "<<p.points.size()<<endl;
  points.resize(points.size()+p.points.size());
  for(size_t i=0;i<p.points.size();i++) {
    points[i+offset].t = p.points[i].t+toffset;
    points[i+offset].x = p.points[i].x;
  }
}

void LinearPath::Append(const Vector& x,Real dt)
{
  ControlPoint cp;
  cp.x=x;
  if(!points.empty()) cp.t=EndTime()+dt;
  else cp.t=dt;
  points.push_back(cp);
}

bool LinearPath::Read(File& f)
{
  int numPoints=0;
  if(!ReadFile(f,numPoints)) return false;
  points.resize(numPoints);
  for(int i=0;i<numPoints;i++) {
    if(!ReadFile(f,points[i].t)) return false;
    if(!points[i].x.Read(f)) return false;
  } 
  return true;
}

bool LinearPath::Write(File& f) const
{
  int numPoints=(int)points.size();
  if(!WriteFile(f,numPoints)) return false;
  for(int i=0;i<numPoints;i++) {
    if(!WriteFile(f,points[i].t)) return false;
    if(!points[i].x.Write(f)) return false;
  } 
  return true;
}


