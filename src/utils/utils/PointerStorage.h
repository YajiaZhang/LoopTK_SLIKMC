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

#ifndef POINTER_STORAGE_H
#define POINTER_STORAGE_H

#include "PointerWrapper.h"
#include "Pointers.h"

namespace Pointers {

/****************************************************************
 * Pointers::Storage
 *
 * Semantics:
 * Acts like a pointer, but:
 * - Deletes the object on destruction
 * - Copies the object on copy constructor, or = operator
 * - Takes ownership of pointers given through constructor or =
 * - Whenever assigning to a new pointer, deletes owned pointers.
 ****************************************************************/

template <class T>
class Storage : public Wrapper<T>
{
public:
  inline Storage() {}
  inline Storage(T* p):Wrapper<T>(p) {}
  inline Storage(const Storage<T>& p) { if(p.ptr) this->ptr=copy(p.ptr); }
  virtual ~Storage() { clear(); }
  void clear() { if(this->ptr) destroy(this->ptr); this->ptr=NULL; }
  const Storage<T>& operator = (const Storage<T>& p) { clear(); if(p.ptr) this->ptr=copy(p.ptr); return *this; }
  const Storage<T>& operator = (T* p) { clear(); this->ptr=p; return *this; }
};

} //namespace Pointers

#endif
