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

#ifndef UNIQUE_POINTER_H
#define UNIQUE_POINTER_H

#include "Pointers.h"
#include "PointerWrapper.h"

namespace Pointers {

/***********************************************
 * Pointers::Unique
 *
 * A class that maintains unique "ownership" of
 * dynamically allocated heap data.
 * Acts like a normal pointer, except for the 
 * following:
 * ~(destructor),clear() deletes owned objects.
 * attach(P) clears and attaches P to A.
 * release() returns the owned object, and releases
 *   it from ownership.
 * steal(A) clears, and takes ownership of A's
 *   objects.
 * swap(A) swaps the owned objects.
 *
 * Operator semantics:
 * =(T*) is equivalent to attach().
 * =(Unique&) is equivalent to steal().
 * The rest are the same as regular pointers.
 ***********************************************/
template <class T>
class Unique : public Wrapper<T>
{
public:
  Unique() {}
  Unique(Unique& p) { steal(p); } 
  Unique(T* p):Wrapper<T>(p) {}
  ~Unique() { clear(); }
  inline void clear() { if(this->ptr) destroy(this->ptr); this->ptr=NULL;; }
  inline void attach(T* p) { clear(); this->ptr=p; }
  inline T* release() { T* p=this->ptr; this->ptr=NULL; return p; }
  inline void steal(Unique<T>& p) { clear(); std::swap(this->ptr,p.ptr); }
  inline void swap(Unique<T>& p) { std::swap(this->ptr,p.ptr); }
  inline void operator=(T* p) { attach(p); } 
  inline void operator=(Unique<T>& p) { steal(p); } 
};

} //namespace Pointers

namespace std {
  template <class T>
  inline void swap(Pointers::Unique<T>& a,Pointers::Unique<T>& b) {
    a.swap(b);
  }
} //namespace std

#endif
