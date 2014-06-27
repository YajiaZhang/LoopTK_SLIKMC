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

#ifndef BASIC_MYHASHTABLE_H
#define BASIC_MYHASHTABLE_H

#include <debug.h>
#include <utils.h>
#include <stdlib.h>
#include <memory.h>

template <class type>
class myhashtable
{
public:
	myhashtable();
	myhashtable(unsigned int size);
	~myhashtable();

	void resize(unsigned int size);
	void cleanup();

	void insert(const type&);
	bool remove(type&);
	bool contains(type&);

	unsigned int (*HashProc) (const type&);

private:
	struct hash_entry
	{
		type data;
		hash_entry* next;
	};

	hash_entry** entries;
	unsigned int num_entries;
};

template <class type>
myhashtable<type>::myhashtable(unsigned int size)
:entries(NULL), num_entries(0)
{
	resize(size);
}

template <class type>
myhashtable<type>::~myhashtable()
{
	cleanup();
}

template <class type>
void myhashtable<type>::resize(unsigned int size)
{
	cleanup();
	num_entries=size;
	entries = new hash_entry* [num_entries];
	memset(entries, 0, num_entries*sizeof(hash_entry*));
}

template <class type>
void myhashtable<type>::cleanup()
{
	hash_entry *temp, *prev;
	for(unsigned int i=0; i<num_entries; i++)
	{
		temp=entries[i];
		while(temp)
		{
			prev=temp->next;
			delete(temp);
			temp = prev;
		}
	}
	SafeArrayDelete(entries);
	num_entries = 0;
}

template <class type>
bool myhashtable<type>::contains(type& data)
{
	hash_entry* e = entries[(HashProc(data)%num_entries)];
	while(e)
	{
		if(e->data == data)
		{
			data = e->data;
			return true;
		}
		e = e->next;
	}

	return false;
}

template <class type>
void myhashtable<type>::insert(const type& data)
{
	unsigned int i = (HashProc(data)%num_entries);

	hash_entry* e = new hash_entry;
	e->data = data;

	if(entries[i])
		temp->next = entries[i];

	entries[i]=temp;
}

template <class type>
bool myhashtable<type>::remove(type& data)
{
	int hashslot = (HashProc(data)%num_entries);
	hash_entry* e = entries[hashslot];
	if(e->data == data)
	{
		data = e->data;
		entries[hashslot] = e->next;
		delete e;
		return true;
	}
	hash_entry* prev = e;
	while(e)
	{
		if(e->data == data)
		{
			data = e->data;
			prev->next = e->next;
			delete e;
			return true;
		}
		prev = e;
		e = e->next;
	}

	return false;
}

#endif
