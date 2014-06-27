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

#ifndef UTILS_TRACE_H
#define UTILS_TRACE_H

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <myfile.h>
#include "Timer.h"

struct TraceFunctionCall;

struct TraceFunction
{
  std::string name;  
  int numCalls;
  double time;
};

struct TraceItem
{
  TraceItem();
  TraceFunctionCall* call;
  std::string text;
};

struct TraceFunctionCall
{
  TraceFunctionCall();
  TraceFunctionCall(TraceFunction* t,TraceFunctionCall* p);
  ~TraceFunctionCall();
  void ClearChildren();

  TraceFunction* type;
  TraceFunctionCall* parent;
  std::string args;
  std::vector<TraceItem> children;
  std::string ret;
  double calltime,endtime;
};

class Trace
{
public:
  Trace();
  ~Trace();
  void Clear();
  void ResetLoop();
  bool Load(const char* fn);
  bool Save(const char* fn);
  void Dump(std::ostream& out=std::cout) const;

  void Call(const char* function,const char* args=NULL);
  void CallFmt(const char* function,const char* fmt,...);
  void EndCall(const char* function,const char* ret=NULL);
  void EndCallFmt(const char* function,const char* fmt,...);
  void Log(const char* txt);

private:
  bool LoadIter(File& f,TraceItem& item,TraceFunctionCall* parent);
  bool SaveIter(File& f,const TraceItem& item);
  void DumpIter(std::ostream& out,const TraceFunctionCall* call,int depth) const;
  TraceFunction* FindFunction(const char* func);
  TraceFunctionCall* FindParentIter(TraceFunctionCall* call,TraceFunction* func);
  void BeginCall(TraceFunctionCall* call);
  void EndCall(TraceFunctionCall* call);

public:
  TraceFunctionCall root;
  TraceFunctionCall* cur;
  std::list<TraceFunction> funcs;
  Timer timer;
};

#endif
