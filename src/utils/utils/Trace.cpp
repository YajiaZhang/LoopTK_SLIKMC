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

#include "Trace.h"
#include <iostream>
#include <utils.h>
#include <stdarg.h>
#include <errors.h>
using namespace std;

#define MAXBUF 4096

void Indent(std::ostream& out,int indent)
{
  for(int i=0;i<indent;i++) out<<" ";
}


TraceItem::TraceItem()
:call(NULL)
{}

TraceFunctionCall::TraceFunctionCall()
:type(NULL),parent(NULL),calltime(0),endtime(0)
{}

TraceFunctionCall::TraceFunctionCall(TraceFunction* t,TraceFunctionCall* p)
:type(t),parent(p),calltime(0),endtime(0)
{}

TraceFunctionCall::~TraceFunctionCall()
{
  ClearChildren();
}

void TraceFunctionCall::ClearChildren()
{
  for(size_t i=0;i<children.size();i++)
    SafeDelete(children[i].call);
  children.clear();
}

Trace::Trace()
{
  cur = &root;
}

Trace::~Trace()
{
  Clear();
}

void Trace::Clear()
{
  timer.Reset();
  root.ClearChildren();
  funcs.clear();
  cur = &root;
}

void Trace::ResetLoop()
{
  if(cur != &root) {
    cout<<"Trace::ResetLoop: Warning, there looks like an unended call in the trace log"<<endl;
    abort();
  }
  timer.Reset();
  root.ClearChildren();
  cur = &root;
}

bool Trace::Load(const char* fn)
{
  Clear();
  File f;
  if(!f.Open(fn,FILEREAD)) return false;
  size_t n;
  char buf[MAXBUF];

  if(!ReadFile(f,n)) return false;
  TraceFunction tf;
  for(size_t i=0;i<n;i++) {
    if(!f.ReadString(buf,MAXBUF)) return false;
    tf.name = buf;
    if(!ReadFile(f,tf.numCalls)) return false;
    if(!ReadFile(f,tf.time)) return false;
    funcs.push_back(tf);
  }

  if(!ReadFile(f,n)) return false;
  root.children.resize(n);
  for(size_t i=0;i<root.children.size();i++)
    if(!LoadIter(f,root.children[i],&root)) return false;
  if(!ReadFile(f,root.calltime)) return false;
  if(!ReadFile(f,root.endtime)) return false;
  return true; 
}

bool Trace::LoadIter(File& f,TraceItem& item,TraceFunctionCall* parent)
{
  char type;
  char buf[MAXBUF];
  if(!ReadFile(f,type)) return false;
  if(type == 'c') {
    if(!f.ReadString(buf,MAXBUF)) return false;
    TraceFunction* func = FindFunction(buf);
    item.call = new TraceFunctionCall(func,parent);
    TraceFunctionCall* call = item.call;
    if(!f.ReadString(buf,MAXBUF)) return false;
    call->args = buf;
    size_t n;
    if(!ReadFile(f,n)) return false;
    call->children.resize(n);
    for(size_t i=0;i<n;i++)
      if(!LoadIter(f,call->children[i],call)) return false;
    if(!f.ReadString(buf,MAXBUF)) return false;
    call->ret = buf;
    if(!ReadFile(f,call->calltime)) return false;
    if(!ReadFile(f,call->endtime)) return false;
    return true;
  }
  else if(type == 'x') {
    if(!f.ReadString(buf,MAXBUF)) return false;
    item.text = buf;
    return true;
  }
  cout<<"Trace::Load(): Wrong type, got "<<type<<endl;
  return false;
}

bool Trace::Save(const char* fn) 
{
  File f;
  if(!f.Open(fn,FILEWRITE)) return false;
  size_t n=funcs.size();
  if(!WriteFile(f,n)) return false;
  for(list<TraceFunction>::iterator i=funcs.begin();i!=funcs.end();i++) {
    TraceFunction* func = &(*i);
    if(!f.WriteString(func->name.c_str())) return false;
    if(!WriteFile(f,func->numCalls)) return false;
    if(!WriteFile(f,func->time)) return false;
  }

  n=root.children.size();
  if(!WriteFile(f,n)) return false;
  for(size_t i=0;i<root.children.size();i++)
    if(!SaveIter(f,root.children[i])) return false;
  if(!WriteFile(f,root.calltime)) return false;
  if(!WriteFile(f,root.endtime)) return false;
  return true; 
}

bool Trace::SaveIter(File& f,const TraceItem& item)
{
  if(item.call) {
    WriteFile(f,'c');
    TraceFunctionCall* call = item.call;
    Assert(call->type != NULL);
    if(!f.WriteString(call->type->name.c_str())) return false;
    if(!f.WriteString(call->args.c_str())) return false;
    size_t n=call->children.size();
    if(!WriteFile(f,n)) return false;
    for(size_t i=0;i<call->children.size();i++)
      if(!SaveIter(f,call->children[i])) return false;
    if(!f.WriteString(call->ret.c_str())) return false;
    if(!WriteFile(f,call->calltime)) return false;
    if(!WriteFile(f,call->endtime)) return false;
    return true;
  }
  else {
    WriteFile(f,'x');
    return f.WriteString(item.text.c_str());
  }
}

void Trace::Dump(ostream& out) const
{
  DumpIter(out,&root,-2);
}

void Trace::DumpIter(ostream& out,const TraceFunctionCall* call,int indent) const
{
  if(call->type) {
    Indent(out,indent); out<<call->type->name<<"("<<call->args<<") : "<<call->calltime<<endl;
  }
  for(size_t i=0;i<call->children.size();i++) {
    const TraceItem& item = call->children[i];
    if(item.call)
      DumpIter(out,item.call,indent+2);
    else {
      Indent(out,indent+2); out<<item.text<<endl;
    }
  }
  if(call->type) {
    Indent(out,indent); 
    if(!call->ret.empty())
      out<<"~"<<call->type->name<<" ("<<call->ret<<") : "<<call->endtime<<endl;
  }
  if(call == cur) {
    Indent(out,indent); out<<"*** current position ***"<<endl;
  }
}

void Trace::CallFmt(const char* function,const char* fmt,...)
{
  char buf [MAXBUF];
  va_list args;
	va_start(args, fmt);
#ifdef WIN32
	_vsnprintf(buf, MAXBUF, fmt, args);
#else
  vsnprintf(buf, MAXBUF, fmt, args);
#endif
  Call(function,buf);
}

void Trace::Call(const char* function,const char* args)
{
  TraceFunction* func = FindFunction(function);
  if(!func) {
    funcs.resize(funcs.size()+1);
    func = &funcs.back();
    func->name = function;
    func->numCalls = 0;
    func->time = 0;
  }
  TraceItem item;
  item.call = new TraceFunctionCall(func,cur);
  cur->children.push_back(item);
  BeginCall(item.call);
  cur = item.call;
}

void Trace::EndCallFmt(const char* function,const char* fmt,...)
{
  char buf [MAXBUF];
  va_list args;
	va_start(args, fmt);
#ifdef WIN32
	_vsnprintf(buf, MAXBUF, fmt, args);
#else
  vsnprintf(buf, MAXBUF, fmt, args);
#endif
  EndCall(function,buf);
}

void Trace::EndCall(const char* function,const char* ret)
{
  TraceFunction* func = FindFunction(function);
  if(!func) {
    FatalError("Trace::EndCall(): attempted to end a nonexistent function");
    //Dump(cout);
  }
  TraceFunctionCall* pcall=FindParentIter(cur,func);
  if(!pcall) {
    FatalError("Trace::EndCall(): Fatal error, attempted to end a non-parent function");
    //Dump(cout);
  }
  while(cur != pcall) {
    cur->ret = "implicit";
    EndCall(cur);
    cur = cur->parent;
  }
  if(ret) cur->ret = ret;
  else cur->ret = "void";
  EndCall(cur);
  cur = cur->parent;
}

void Trace::Log(const char* txt)
{
  TraceItem item; item.text = txt;
  cur->children.push_back(item);
}

TraceFunction* Trace::FindFunction(const char* func)
{
  //really stupid, linear time algorithm for now
  for(list<TraceFunction>::iterator i=funcs.begin();i!=funcs.end();i++)
    if(0 == strcmp(i->name.c_str(),func)) return &(*i);
  return NULL;
}

TraceFunctionCall* Trace::FindParentIter(TraceFunctionCall* call,TraceFunction* func)
{
  if(!call) return NULL;
  if(call->type == func) return call;
  return FindParentIter(call->parent,func);
}

void Trace::BeginCall(TraceFunctionCall* call)
{
  call->calltime = timer.ElapsedTime();
  call->type->numCalls++;
}

void Trace::EndCall(TraceFunctionCall* call)
{
  call->endtime = timer.ElapsedTime();
  call->type->time += (call->endtime-call->calltime);
}
