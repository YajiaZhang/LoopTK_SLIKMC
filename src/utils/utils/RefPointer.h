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

#ifndef REF_POINTER_H
#define REF_POINTER_H

#include <stdlib.h>

class RefObjectBase
{
public:
	RefObjectBase();
	void Ref();
	void Unref();
	void UnrefNoDelete();
	inline int NumRefs() const { return numRefs; }

protected:
	//protected so you can't delete this automatically
	virtual ~RefObjectBase();

private:
	int numRefs;
};

//obj must define the methods Ref and Unref
template <class obj>
class RefPointer
{
public:
	typedef RefPointer<obj> Pointer;

	RefPointer()
	:object(NULL)
	{}

	RefPointer(obj* obj_in)
	:object(obj_in)
	{ if(object != NULL) object->Ref();	}

	RefPointer(const RefPointer& rhs)
	:object(rhs.object)
	{ if(object != NULL) object->Ref(); }

	~RefPointer()
	{ if(object != NULL) object->Unref(); }

	inline operator obj* () const { return object; }
	inline obj* operator ->() const { return object; }
	inline const Pointer& operator = (const Pointer& rhs) { return operator=(rhs.object); }
	inline const Pointer& operator = (obj* obj_in) {
		if(object) object->Unref();
		object=obj_in;
		if(object) object->Ref();
		return *this;
	}
	inline bool operator == (const Pointer& rhs) const { return (object==rhs.object); }
	inline bool operator != (const Pointer& rhs) const { return (object!=rhs.object); }
	inline bool isNull() const { return (object==NULL); }

private:
	obj* object;
};

//RefStorage
//stores a RefPointer that keeps an object persistent
//obj must define the methods Ref and Unref
template <class obj>
class RefStorage
{
public:
	typedef RefStorage<obj> Storage;
	typedef RefPointer<obj> Pointer;

	RefStorage()
	:pointer(new obj)
	{}

	RefStorage(const Pointer& ptr)
	:pointer(ptr)
	{}

	RefStorage(const Storage& rhs)
	:pointer(rhs.pointer)
	{}

	Pointer GetPointer() const { return pointer; }

	const Storage& operator = (const Storage& storage)
	{ pointer=storage.pointer; return *this; }
	const Storage& operator = (const Pointer& ptr)
	{ pointer=ptr; return *this; }

	operator const obj& () const { return *pointer; }
	void Copy(const Storage& storage) { pointer->Copy(*storage.pointer); }
	void Copy(const obj& storage) { pointer->Copy(storage); }

private:
	Pointer pointer;
};



#endif
