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

#ifndef BASIC_MYQUEUE_H
#define BASIC_MYQUEUE_H

#include "mylist.h"

template<class type>
class myqueue : public mylist<type>
{
public:
	virtual void enqueue(const type& t)
	{
		push_back(t);
	}
	virtual void dequeue()
	{
		pop_front();
	}
	virtual void dequeue(type& t)
	{
		pop_front(t);
	}
};

//assumes a member variable 'priority' is defined in 'type'
template<class type>
class mypriorityqueue : public myqueue<type>
{
public:
	void enqueue(const type& t)
	{
		type_node* temp = new type_node;
		temp->t = t;
		temp->prev = NULL;
		temp->next = NULL;
		num_nodes++;
		insert(temp);
	}

	void refresh(type* t)
	{
		type_node* temp = (type_node*)t;
		detach(temp);
		insert(temp);
	}

	void refresh_all()
	{
		type_node *frontsave = _front, *backsave = _back;
		_front = _back = NULL;

		type_node *temp = frontsave, *next;
		while(temp)
		{
			next = temp->next;
			insert(temp);
			temp = next;
		}
	}

private:
	void insert(type_node* t)
	{
		type_node* temp = _front;
		while(temp)
		{
			if(temp->t.priority <= t->t.priority)
			{
				t->prev = temp->prev;
				t->next = temp;
				if(temp->prev)
					temp->prev->next = t;
				else
					_front = t;
				temp->prev = t;
				return;
			}
			temp = temp->next;
		}
		//insert in back
		t->prev = _back;
		t->next = NULL;
		if(_front == NULL)
			_front = t;
		if(_back == NULL)
			_back = t;
		else
			_back->next = t;
		_back = t;
	}
	void detach(type_node* t)
	{
		if(t->prev)
			t->prev->next = t->next;
		else
			_front = t->next;
		if(t->next)
			t->next->prev = t->prev;
		else
			_back = t->prev;
	}
};

#endif
